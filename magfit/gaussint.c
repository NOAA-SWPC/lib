/*
 * gaussint.c
 *
 * Fit internal potential field based on Gauss coefficients
 * to magnetic vector data
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <satdata/satdata.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_blas.h>

#include <common/common.h>
#include <common/interp.h>
#include <common/oct.h>

#include "green.h"
#include "track.h"

#include "magfit.h"

typedef struct
{
  size_t n;         /* total number of measurements in system */
  size_t p;         /* number of coefficients */
  size_t nmax;      /* maximum number of measurements in LS system */

  size_t lmax;      /* maximum spherical harmonic degree */
  size_t mmax;      /* maximum spherical harmonic order */

  gsl_matrix *X;    /* LS matrix */
  gsl_vector *rhs;  /* rhs vector */
  gsl_vector *wts;  /* weight vector */
  gsl_matrix *cov;  /* covariance matrix */
  gsl_vector *c;    /* solution vector */

  gsl_multifit_linear_workspace *multifit_p;
  gsl_multilarge_linear_workspace *multilarge_p;
  green_workspace *green_workspace_p;
} gaussint_state_t;

static void *gaussint_alloc(const void * params);
static void gaussint_free(void * vstate);
static int gaussint_reset(void * vstate);
static size_t gaussint_ncoeff(void * vstate);
static int gaussint_add_datum(const double t, const double r, const double theta, const double phi,
                              const double qdlat, const double B[3], void * vstate);
static int gaussint_fit(double * rnorm, double * snorm, void * vstate);
static int gaussint_eval_B(const double t, const double r, const double theta, const double phi,
                           double B[3], void * vstate);
static int gaussint_eval_J(const double r, const double theta, const double phi,
                           double J[3], void * vstate);
static double gaussint_eval_chi(const double theta, const double phi, void * vstate);

static int build_matrix_row(const double r, const double theta, const double phi,
                            gsl_vector *X, gsl_vector *Y, gsl_vector *Z,
                            gaussint_state_t *state);

/*
gaussint_alloc()
  Allocate gaussint workspace

Inputs: params - parameters

Return: pointer to workspace
*/

static void *
gaussint_alloc(const void * params)
{
  gaussint_state_t *state;

  (void) params;

  state = calloc(1, sizeof(gaussint_state_t));
  if (!state)
    return 0;

  state->lmax = 60;
  state->mmax = 30;
  state->nmax = 200000;
  state->n = 0;

  state->green_workspace_p = green_alloc(state->lmax, state->mmax, R_EARTH_KM);

  state->p = green_nnm(state->green_workspace_p);

  state->X = gsl_matrix_alloc(state->nmax, state->p);
  state->c = gsl_vector_alloc(state->p);
  state->rhs = gsl_vector_alloc(state->nmax);
  state->wts = gsl_vector_alloc(state->nmax);
  state->cov = gsl_matrix_alloc(state->p, state->p);
  state->multifit_p = gsl_multifit_linear_alloc(state->nmax, state->p);
  state->multilarge_p = gsl_multilarge_linear_alloc(gsl_multilarge_linear_normal, state->p);

  return state;
}

static void
gaussint_free(void * vstate)
{
  gaussint_state_t *state = (gaussint_state_t *) vstate;

  if (state->X)
    gsl_matrix_free(state->X);

  if (state->c)
    gsl_vector_free(state->c);

  if (state->rhs)
    gsl_vector_free(state->rhs);

  if (state->wts)
    gsl_vector_free(state->wts);

  if (state->cov)
    gsl_matrix_free(state->cov);

  if (state->multifit_p)
    gsl_multifit_linear_free(state->multifit_p);

  if (state->multilarge_p)
    gsl_multilarge_linear_free(state->multilarge_p);

  if (state->green_workspace_p)
    green_free(state->green_workspace_p);

  free(state);
}

/* reset workspace to work on new data set */
static int
gaussint_reset(void * vstate)
{
  gaussint_state_t *state = (gaussint_state_t *) vstate;
  state->n = 0;
  return 0;
}

static size_t
gaussint_ncoeff(void * vstate)
{
  gaussint_state_t *state = (gaussint_state_t *) vstate;
  return state->p;
}

/*
gaussint_add_datum()
  Add single vector measurement to LS system

Inputs: t      - timestamp (CDF_EPOCH)
        r      - radius (km)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        qdlat  - QD latitude (degrees)
        B      - magnetic field vector NEC (nT)
        vstate - state

Return: success/error

Notes:
1) state->n is updated with the number of total data added
*/

static int
gaussint_add_datum(const double t, const double r, const double theta, const double phi,
                   const double qdlat, const double B[3], void * vstate)
{
  gaussint_state_t *state = (gaussint_state_t *) vstate;
  size_t rowidx = state->n;
  double wi = 1.0;
  gsl_vector_view vx = gsl_matrix_row(state->X, rowidx);
  gsl_vector_view vy = gsl_matrix_row(state->X, rowidx + 1);
  gsl_vector_view vz = gsl_matrix_row(state->X, rowidx + 2);

  (void) t;
  (void) qdlat;

  /* set rhs vector */
  gsl_vector_set(state->rhs, rowidx, B[0]);
  gsl_vector_set(state->rhs, rowidx + 1, B[1]);
  gsl_vector_set(state->rhs, rowidx + 2, B[2]);

  /* set weight vector */
  gsl_vector_set(state->wts, rowidx, wi);
  gsl_vector_set(state->wts, rowidx + 1, wi);
  gsl_vector_set(state->wts, rowidx + 2, wi);

  /* build 3 rows of the LS matrix */
  build_matrix_row(r, theta, phi, &vx.vector, &vy.vector, &vz.vector, state);
  rowidx += 3;

  state->n = rowidx;

  return GSL_SUCCESS;
}

/*
gaussint_fit()
  Fit model to previously added tracks

Inputs: rnorm  - residual norm || y - A x ||
        snorm  - solution norm || x ||
        vstate - state

Return: success/error

Notes:
1) Data must be added to workspace via gaussint_add_datum()
*/

static int
gaussint_fit(double * rnorm, double * snorm, void * vstate)
{
  gaussint_state_t *state = (gaussint_state_t *) vstate;
  gsl_matrix_view A = gsl_matrix_submatrix(state->X, 0, 0, state->n, state->p);
  gsl_vector_view b = gsl_vector_subvector(state->rhs, 0, state->n);
  gsl_vector_view wts = gsl_vector_subvector(state->wts, 0, state->n);
  double rcond;

  if (state->n < state->p)
    return -1;

  fprintf(stderr, "\n");
  fprintf(stderr, "\t n = %zu\n", state->n);
  fprintf(stderr, "\t p = %zu\n", state->p);

#if 0

  {
    double chisq;

    /* solve system */
    gsl_multifit_wlinear(&A.matrix, &wts.vector, &b.vector, state->c, state->cov, &chisq, state->multifit_p);

    rcond = gsl_multifit_linear_rcond(state->multifit_p);
    *rnorm = sqrt(chisq);
    *snorm = gsl_blas_dnrm2(state->c);
  }

#else

  /* apply weights */
  gsl_multifit_linear_applyW(&A.matrix, &wts.vector, &b.vector, &A.matrix, &b.vector);

  /* fold (A,b) in to LS system */
  gsl_multilarge_linear_accumulate(&A.matrix, &b.vector, state->multilarge_p);

  /* solve system */
  gsl_multilarge_linear_solve(0.0, state->c, rnorm, snorm, state->multilarge_p);

  /* compute condition number */
  gsl_multilarge_linear_rcond(&rcond, state->multilarge_p);

#endif

  fprintf(stderr, "\t rnorm = %g\n", *rnorm);
  fprintf(stderr, "\t snorm = %g\n", *snorm);
  fprintf(stderr, "\t cond(X) = %g\n", 1.0 / rcond);

  return 0;
}

/*
gaussint_eval_B()
  Evaluate magnetic field at a given (r,theta,phi) using
previously computed coefficients

Inputs: t      - timestamp (CDF_EPOCH)
        r      - radius (km)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        B      - (output) magnetic field vector (nT)
        vstate - state

Notes:
1) state->c must contain fit coefficients
*/

static int
gaussint_eval_B(const double t, const double r, const double theta, const double phi,
                double B[3], void * vstate)
{
  int status = GSL_SUCCESS;
  gaussint_state_t *state = (gaussint_state_t *) vstate;
  gsl_vector_view vx = gsl_matrix_row(state->X, 0);
  gsl_vector_view vy = gsl_matrix_row(state->X, 1);
  gsl_vector_view vz = gsl_matrix_row(state->X, 2);

  (void) t;

  build_matrix_row(r, theta, phi, &vx.vector, &vy.vector, &vz.vector, state);

  gsl_blas_ddot(&vx.vector, state->c, &B[0]);
  gsl_blas_ddot(&vy.vector, state->c, &B[1]);
  gsl_blas_ddot(&vz.vector, state->c, &B[2]);

  return status;
}

/*
gaussint_eval_J()
  Evaluate current density at a given (r,theta,phi) using
previously computed coefficients

Inputs: r      - radius (km)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        J      - (output) current density vector [A/km]
        vstate - workspace

Notes:
1) state->c must contain coefficients
*/

static int
gaussint_eval_J(const double r, const double theta, const double phi,
                double J[3], void * vstate)
{
  int status;
  gaussint_state_t *state = (gaussint_state_t *) vstate;

  (void) r; /* unused parameter */

  status = green_eval_sheet_int(R_EARTH_KM + 110.0, theta, phi, state->c, J,
                                state->green_workspace_p);

  return status;
}

/*
gaussint_eval_chi()
  Evaluate current stream function at a given (theta,phi) using
previously computed coefficients

Inputs: theta  - colatitude (radians)
        phi    - longitude (radians)
        vstate - workspace

Return: current stream function chi in kA/nT

Notes:
1) state->c must contain coefficients
*/

static double
gaussint_eval_chi(const double theta, const double phi, void * vstate)
{
  gaussint_state_t *state = (gaussint_state_t *) vstate;
  double chi;

  chi = green_eval_chi_int(R_EARTH_KM + 110.0, theta, phi, state->c, state->green_workspace_p);

  return chi;
}

static int
build_matrix_row(const double r, const double theta, const double phi,
                 gsl_vector *X, gsl_vector *Y, gsl_vector *Z,
                 gaussint_state_t *state)
{
  int s = green_calc_int(r, theta, phi, X->data, Y->data, Z->data, state->green_workspace_p);
  return s;
}

static const magfit_type gaussint_type =
{
  "gaussint",
  gaussint_alloc,
  gaussint_reset,
  gaussint_ncoeff,
  gaussint_add_datum,
  gaussint_fit,
  gaussint_eval_B,
  gaussint_eval_J,
  gaussint_eval_chi,
  gaussint_free
};

const magfit_type *magfit_gaussint = &gaussint_type;
