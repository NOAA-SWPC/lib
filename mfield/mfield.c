/*
 * mfield.c
 *
 * This module contains code to perform magnetic main field modeling
 *
 * Calling sequences:
 * 1. mfield_data_alloc     - allocate mfield_data_workspace structure
 * 2. mfield_data_copy      - copy all satellite data to internal
 *                            structure, ignoring flagged data
 * 3. mfield_data_map       - print out map of spatial coverage of data
 * 4. mfield_alloc          - allocate mfield_workspace
 * 5. mfield_init           - initialize various internal parameters
 *                            (time scaling, etc)
 * 6. mfield_set_damping    - set coefficent damping parameters (requires
 *                            time scaling parameters)
 * 7. mfield_calc_linear    - build linear LS matrix and solve system
 * 8. mfield_calc_nonlinear - solve nonlinear LS system; the initial
 *                            guess for the nonlinear procedure is the
 *                            linearized solution
 * 9. mfield_free
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <errno.h>
#include <string.h>

#include <satdata/satdata.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sort_vector.h>

#include "mfield_green.h"

#include "common.h"
#include "euler.h"
#include "lls.h"
#include "mfield.h"
#include "oct.h"
#include "track_weight.h"

/*
 * scale time to dimensionless units for SV/SA terms - this
 * helps improve the condition of the least squares A^T A matrix
 */
#define MFIELD_SCALE_TIME            1

#define MFIELD_REGULARIZE            1

/* weighting factors for data (multiplies spatial weights) */
#define MFIELD_WEIGHT_F              (1.0)
#define MFIELD_WEIGHT_X              (1.0)
#define MFIELD_WEIGHT_Y              (1.0)
#define MFIELD_WEIGHT_Z              (1.0)

static int mfield_green(const double r, const double theta, const double phi,
                        mfield_workspace *w);
static int mfield_eval_g(const double t, const double r, const double theta, const double phi,
                         const gsl_vector *c, double B[4], mfield_workspace *w);
static int mfield_print_residuals(const int header, FILE *fp,
                                  const gsl_vector *r,
                                  const gsl_vector *w_rob,
                                  const gsl_vector *w_spatial);
static int mfield_update_histogram(const gsl_vector *r, gsl_histogram *h);
static int mfield_print_histogram(FILE *fp, gsl_histogram *h);

#include "mfield_euler.c"
#include "mfield_nonlinear.c"
#include "lapack_inverse.c"

/*
mfield_alloc()
  Allocate a mfield workspace

Inputs: params - model parameters
                 epoch   - model epoch t0 (years)
                 R       - reference radius (km)
                 nmax_mf - maximum spherical harmonic degree for MF
                 nmax_sv - maximum spherical harmonic degree for SV
                 nmax_sa - maximum spherical harmonic degree for SA
                 nsat    - number of different satellites
                 flags   - flags
*/

mfield_workspace *
mfield_alloc(const mfield_parameters *params)
{
  mfield_workspace *w;
  const size_t plm_size = gsl_sf_legendre_array_n(params->nmax_mf);
  const size_t ntheta = 100;
  const size_t nphi = 100;

  w = calloc(1, sizeof(mfield_workspace));
  if (!w)
    return 0;

  w->nsat = params->nsat;
  w->epoch = params->epoch;
  w->R = params->R;
  w->nmax_mf = params->nmax_mf;
  w->nmax_sv = params->nmax_sv;
  w->nmax_sa = params->nmax_sa;
  w->data_workspace_p = params->mfield_data_p;

  w->params = *params;

  w->weight_workspace_p = track_weight_alloc(ntheta, nphi);

  w->nbins_euler = calloc(1, w->nsat * sizeof(size_t));
  w->offset_euler = calloc(1, w->nsat * sizeof(size_t));
  if (!w->nbins_euler || !w->offset_euler)
    {
      mfield_free(w);
      return 0;
    }

  /*
   * Add up all the contributions to the coefficient vector, which will
   * be partitioned as:
   *
   * c = [ MF | SV | SA | Euler | external ]
   */

  /* subtract 1 to exclude the (0,0) coefficient */
  w->nnm_mf = (w->nmax_mf + 1) * (w->nmax_mf + 1) - 1;
  w->nnm_sv = (w->nmax_sv + 1) * (w->nmax_sv + 1) - 1;
  w->nnm_sa = (w->nmax_sa + 1) * (w->nmax_sa + 1) - 1;

#if !MFIELD_FIT_SECVAR
  w->nnm_sv = 0;
#endif

#if !MFIELD_FIT_SECACC
  w->nnm_sa = 0;
#endif

  /* total (internal) model coefficients */
  w->p_int = w->nnm_mf + w->nnm_sv + w->nnm_sa;
  w->p = w->p_int;

#if MFIELD_FIT_EULER
  if (w->data_workspace_p)
    {
      size_t i;
      size_t sum = 0;

      /* compute total number of Euler bins for each satellite */
      for (i = 0; i < w->nsat; ++i)
        {
          magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
          double t0, t1, dt;

          magdata_t(&t0, &t1, mptr);
          dt = (t1 - t0) / 86400000.0; /* convert to days */

          if (params->euler_period <= 0.0)
            w->nbins_euler[i] = 1;
          else
            w->nbins_euler[i] = (size_t) (dt / params->euler_period) + 1;

          w->offset_euler[i] = 3 * sum;
          sum += w->nbins_euler[i];
        }

      w->neuler = 3 * sum;
      /*w->neuler = 3 * w->nsat;*/
    }
#else
  w->neuler = 0;
#endif

#if MFIELD_FIT_EXTFIELD
  /* allow for 3 years of data plus 1 month of spill */
  w->next = 3 * 366 + 30;
#else
  w->next = 0;
#endif

  w->p += w->neuler + w->next;

  w->sv_offset = w->nnm_mf;
  w->sa_offset = w->sv_offset + w->nnm_sv;
  w->euler_offset = w->sa_offset + w->nnm_sa;
  w->ext_offset = w->euler_offset + w->neuler;

  w->cosmphi = malloc((w->nmax_mf + 1) * sizeof(double));
  w->sinmphi = malloc((w->nmax_mf + 1) * sizeof(double));

  w->Plm = malloc(plm_size * sizeof(double));
  w->dPlm = malloc(plm_size * sizeof(double));
  if (!w->Plm || !w->dPlm)
    {
      mfield_free(w);
      return 0;
    }

  w->c = gsl_vector_calloc(w->p);
  w->c_copy = gsl_vector_alloc(w->p);

  /* covariance matrix of the internal field coefficients */
  w->covar = gsl_matrix_alloc(w->p_int, w->p_int);

  w->dX = malloc(w->nnm_mf * sizeof(double));
  w->dY = malloc(w->nnm_mf * sizeof(double));
  w->dZ = malloc(w->nnm_mf * sizeof(double));

  w->green_workspace_p = mfield_green_alloc(w->nmax_mf, w->R);

  w->diag = gsl_vector_alloc(w->p_int);

  w->hf = gsl_histogram_alloc(500);
  w->hz = gsl_histogram_alloc(500);
  gsl_histogram_set_ranges_uniform(w->hf, -40.0, 40.0);
  gsl_histogram_set_ranges_uniform(w->hz, -40.0, 40.0);

  w->nobs_cnt = 0;

  /* these are computed later in mfield_init() */
  w->t_mu = 0.0;
  w->t_sigma = 1.0;
  w->t0_data = 0.0;

  w->lambda_mf = 0.0;
  w->lambda_sv = 0.0;
  w->lambda_sa = 0.0;

  w->niter = 0;

  /* maximum observations to accumulate at once in LS system */
  w->data_block = MFIELD_BLOCK_SIZE;
  
  /* add factor 4 for (X,Y,Z,F) */
  w->block_J = gsl_matrix_alloc(4 * w->data_block, w->p);
  w->block_f = gsl_vector_alloc(4 * w->data_block);
  w->wts = gsl_vector_alloc(4 * w->data_block);

  w->JTJ_vec = gsl_matrix_alloc(w->p_int, w->p_int);

  w->eigen_workspace_p = gsl_eigen_symm_alloc(w->p);

  w->lambda_diag = gsl_vector_calloc(w->p);
  w->LTL = gsl_vector_calloc(w->p);

  w->block_dX = gsl_matrix_alloc(w->data_block, w->nnm_mf);
  w->block_dY = gsl_matrix_alloc(w->data_block, w->nnm_mf);
  w->block_dZ = gsl_matrix_alloc(w->data_block, w->nnm_mf);
  w->fp_dX = fopen("mat/dX.dat", "w+");
  w->fp_dY = fopen("mat/dY.dat", "w+");
  w->fp_dZ = fopen("mat/dZ.dat", "w+");

  return w;
} /* mfield_alloc() */

void
mfield_free(mfield_workspace *w)
{
  if (w->cosmphi)
    free(w->cosmphi);

  if (w->sinmphi)
    free(w->sinmphi);

  if (w->Plm)
    free(w->Plm);

  if (w->dPlm)
    free(w->dPlm);

  if (w->c)
    gsl_vector_free(w->c);

  if (w->c_copy)
    gsl_vector_free(w->c_copy);

  if (w->covar)
    gsl_matrix_free(w->covar);

  if (w->dX)
    free(w->dX);

  if (w->dY)
    free(w->dY);

  if (w->dZ)
    free(w->dZ);

  if (w->diag)
    gsl_vector_free(w->diag);

  if (w->hf)
    gsl_histogram_free(w->hf);

  if (w->hz)
    gsl_histogram_free(w->hz);

  if (w->mat_dX)
    gsl_matrix_free(w->mat_dX);

  if (w->mat_dY)
    gsl_matrix_free(w->mat_dY);

  if (w->mat_dZ)
    gsl_matrix_free(w->mat_dZ);

  if (w->lambda_diag)
    gsl_vector_free(w->lambda_diag);

  if (w->LTL)
    gsl_vector_free(w->LTL);

  if (w->wts_spatial)
    gsl_vector_free(w->wts_spatial);

  if (w->wts_final)
    gsl_vector_free(w->wts_final);

  if (w->multifit_nlinear_p)
    gsl_multifit_nlinear_free(w->multifit_nlinear_p);

  if (w->nbins_euler)
    free(w->nbins_euler);

  if (w->offset_euler)
    free(w->offset_euler);

  if (w->fvec)
    gsl_vector_free(w->fvec);

  if (w->wfvec)
    gsl_vector_free(w->wfvec);

  if (w->robust_workspace_p)
    gsl_multifit_robust_free(w->robust_workspace_p);

  if (w->green_workspace_p)
    mfield_green_free(w->green_workspace_p);

  if (w->weight_workspace_p)
    track_weight_free(w->weight_workspace_p);

  if (w->block_J)
    gsl_matrix_free(w->block_J);

  if (w->block_f)
    gsl_vector_free(w->block_f);

  if (w->wts)
    gsl_vector_free(w->wts);

  if (w->JTJ_vec)
    gsl_matrix_free(w->JTJ_vec);

  if (w->block_dX)
    gsl_matrix_free(w->block_dX);

  if (w->block_dY)
    gsl_matrix_free(w->block_dY);

  if (w->block_dZ)
    gsl_matrix_free(w->block_dZ);

  if (w->fp_dX)
    fclose(w->fp_dX);

  if (w->fp_dY)
    fclose(w->fp_dY);

  if (w->fp_dZ)
    fclose(w->fp_dZ);

  if (w->nlinear_workspace_p)
    gsl_multilarge_regnlinear_free(w->nlinear_workspace_p);

  if (w->eigen_workspace_p)
    gsl_eigen_symm_free(w->eigen_workspace_p);

  free(w);
} /* mfield_free() */

/*
mfield_copy()
  Make a copy of an mfield workspace that can be used
interchangeably with mfield_eval(). This is useful for
OpenMP applications, since mfield_eval() uses internal
arrays (dX,dY,dZ) so we need separate workspaces for each
thread for evaluation
*/

mfield_workspace *
mfield_copy(const mfield_workspace *w)
{
  mfield_workspace *w_copy;

  w_copy = mfield_alloc(&(w->params));
  if (!w_copy)
    return 0;

  /* copy coefficients */
  gsl_vector_memcpy(w_copy->c, w->c);

  /* copy time scaling parameters */
  w_copy->t_mu = w->t_mu;
  w_copy->t_sigma = w->t_sigma;

  return w_copy;
} /* mfield_copy() */

/*
mfield_reset()
  Reset mfield workspace to be ready for a new
LS iteration
*/

int
mfield_reset(mfield_workspace *w)
{
  int s = 0;

#if 0
  /* reset weight histogram */
  s += weight_reset(w->weight_workspace_p);
#endif

  w->nobs_cnt = 0;

  return s;
}

/*
mfield_init()
  Initialize model

Notes:
1) On input, w->data_workspace_p must be initialized with satellite data

2) On output, w->t_mu and w->t_sigma are initialized if time
scaling is desired

3) On output, w->t0_data is initialized to the time of the first available
data point (CDF_EPOCH)

4) Area/density weights are computed from previously constructed
histogram
*/

int
mfield_init(mfield_workspace *w)
{
  int s = 0;
  size_t i, j;

  /* initialize t_mu and t_sigma */
  mfield_data_init(w->data_workspace_p);

#if MFIELD_SCALE_TIME
  w->t_mu = w->data_workspace_p->t_mu;
  w->t_sigma = w->data_workspace_p->t_sigma;
#endif

  fprintf(stderr, "mfield_init: t_mu    = %g [years]\n", w->t_mu);
  fprintf(stderr, "mfield_init: t_sigma = %g [years]\n", w->t_sigma);

  /* find time of first available data in CDF_EPOCH */
  w->t0_data = w->data_workspace_p->t0_data;

  /* initialize spatial weighting histogram and time scaling */
  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

      for (j = 0; j < mptr->n; ++j)
        {
          const double u = satdata_epoch2year(mptr->t[j]) - w->epoch;

          /* center and scale time */
          mptr->ts[j] = (u - w->t_mu) / w->t_sigma;

          if (mptr->flags[j] & MAGDATA_FLG_DISCARD)
            continue;

#if 0
          if (mptr->flags[j] & MAGDATA_FLG_X)
            track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);
          if (mptr->flags[j] & MAGDATA_FLG_Y)
            track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);
          if (mptr->flags[j] & MAGDATA_FLG_Z)
            track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);
          if (mptr->flags[j] & MAGDATA_FLG_F)
            track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);
#else
          track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);
#endif
        }
    }

  /* compute data weights with histogram */
  track_weight_calc(w->weight_workspace_p);

  mfield_init_nonlinear(w);

  return s;
} /* mfield_init() */

/*
mfield_set_damping()
  Set damping parameters for SV, SA coefficients

Inputs: lambda_sv - damping parameter for SV coefficients (years)
        lambda_sa - damping parameter for SA coefficients (years^2)

Notes:
1) The above units ensure the Tikhonov term ||Lg||^2 has units of nT^2,
similar to the residual vector term ||f||^2
*/

int
mfield_set_damping(const double lambda_sv, const double lambda_sa,
                   mfield_workspace *w)
{
  /* convert to dimensionless units with time scale */
  w->lambda_mf = 0.0;
  w->lambda_sv = lambda_sv / w->t_sigma;
  w->lambda_sa = lambda_sa / (w->t_sigma * w->t_sigma);

  return GSL_SUCCESS;
} /* mfield_set_damping() */

/*
mfield_calc_uncertainties()
  Compute uncertainties in coefficients using the variance-covariance
matrix:

C = [J^T W J]^{-1}

and

W = diag(w_i) = diag(1 / sigma^2)  with

sigma^2 = chi^2 / (n-p+1) for all i

Notes:
1) On output, w->covar contains the covariance matrix
   C = sigma^2 (J^T J)^{-1}
*/

int
mfield_calc_uncertainties(mfield_workspace *w)
{
  int s = GSL_SUCCESS;
  /*
   * compute uncertainties only for internal field coefficients; the
   * external field coefficients could be zero for many days (if no
   * data) leading to a singular Jacobian
   */
  const size_t p = w->p_int;
  gsl_vector * f = gsl_multilarge_regnlinear_residual(w->nlinear_workspace_p);
  gsl_matrix * JTJ = gsl_multilarge_regnlinear_JTJ(w->nlinear_workspace_p);
  gsl_matrix_view m = gsl_matrix_submatrix(JTJ, 0, 0, p, p);
  double dof = (double)f->size - (double)w->p;
  double chisq, sigmasq;

  fprintf(stderr, "\n");

  /* compute chi^2 = f^T f */
  fprintf(stderr, "\t computing final chisq...");
  gsl_blas_ddot(f, f, &chisq);
  fprintf(stderr, "done (chisq = %f)\n", chisq);

  /* use total number of coefficients for this calculation of dof */
  fprintf(stderr, "\t computing final sigmasq...");
  sigmasq = chisq / dof;
  fprintf(stderr, "done (sigmasq = %f)\n", sigmasq);

  /* copy lower triangle of JTJ to covar */
  gsl_matrix_tricpy('L', 1, w->covar, &m.matrix);

  /* compute (J^T J)^{-1} using Cholesky decomposition */
  fprintf(stderr, "\t inverting J^T J matrix...");
  gsl_linalg_cholesky_decomp(w->covar);
  gsl_linalg_cholesky_invert(w->covar);
  fprintf(stderr, "done\n");

  /* scale by sigma^2 */
  gsl_matrix_scale(w->covar, sigmasq);

  return s;
} /* mfield_calc_uncertainties() */

int
mfield_calc_evals(gsl_vector *evals, mfield_workspace *w)
{
  int s;
  gsl_matrix * JTJ = gsl_multilarge_regnlinear_JTJ(w->nlinear_workspace_p);
  gsl_matrix * A = gsl_matrix_alloc(w->p, w->p);

  /* copy lower triangle of JTJ to A */
  gsl_matrix_tricpy('L', 1, A, JTJ);

  s = gsl_eigen_symm(A, evals, w->eigen_workspace_p);
  gsl_sort_vector(evals);

  gsl_matrix_free(A);

  return s;
}

/*
mfield_coeffs()
  After the least squares system has been solved, convert
coefficients to physical units and store in a given vector

Inputs: dir  - direction of transformation
                1: from dimensionless to physical
               -1: from physical to dimensionless
        gin  - input coefficients
        gout - (output) where to store coefficients
               dir = 1:
                 MF coefficients: nT
                 SV coefficients: nT/year
                 SA coefficients: nT/year^2
               dir = -1:
                 MF coefficients: nT
                 SV coefficients: nT/dimensionless_time
                 SA coefficients: nT/dimensionless_time^2
        w    - workspace

Notes:
1) It is allowed for gin = gout
*/

int
mfield_coeffs(const int dir, const gsl_vector *gin, gsl_vector *gout,
              const mfield_workspace *w)
{
  int s = 0;
  size_t n;
  int m;
  double ratio = w->t_mu / w->t_sigma;
  double ratio_sq = ratio * ratio;
  double sigma_inv = 1.0 / w->t_sigma;
  double sigma_inv_sq = sigma_inv * sigma_inv;
  double sigma = w->t_sigma;
  double sigma_sq = sigma * sigma;
  double mu = w->t_mu;
  double mu_sq = mu * mu;

  /* this will copy any external coefficients over */
  gsl_vector_memcpy(gout, gin);

  for (n = 1; n <= w->nmax_mf; ++n)
    {
      int ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);
          double gnm = mfield_get_mf(gin, cidx, w);
          double dgnm = mfield_get_sv(gin, cidx, w);
          double ddgnm = mfield_get_sa(gin, cidx, w);

          if (dir == 1)
            {
              /*
               * From dimensionless to physical units;
               * undo the scaling transformation t~ = (t - t0 - mu) / sigma
               *
               * gnm_orig = gnm - (mu/sigma)*dgnm + 1/2*(mu/sigma)^2*ddgnm
               * dgnm_orig = (dgnm - (mu/sigma)*ddgnm) / sigma
               * ddgnm_orig = ddgnm / sigma^2
               */
              mfield_set_mf(gout, cidx, gnm - ratio * dgnm + 0.5 * ratio_sq * ddgnm, w);
              mfield_set_sv(gout, cidx, (dgnm - ratio * ddgnm) * sigma_inv, w);
              mfield_set_sa(gout, cidx, ddgnm * sigma_inv_sq, w);
            }
          else
            {
              /*
               * From physical to dimensionless units
               *
               * gnm_new = gnm + mu*dgnm + 1/2*mu^2*ddgnm
               * dgnm_new = sigma*(dgnm + mu*ddgnm)
               * ddgnm_new = sigma^2 ddgnm
               */
              mfield_set_mf(gout, cidx, gnm + mu * dgnm + 0.5 * mu_sq * ddgnm, w);
              mfield_set_sv(gout, cidx, sigma * (dgnm + mu * ddgnm), w);
              mfield_set_sa(gout, cidx, ddgnm * sigma_sq, w);
            }
        }
    }

  return s;
} /* mfield_coeffs() */

/*
mfield_eval()
  Evaluate magnetic field model at given point

Inputs: t     - timestamp (CDF_EPOCH)
        r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - geocentric longitude (radians)
        B     - (output) magnetic field
                B[0] = B_x
                B[1] = B_y
                B[2] = B_z
                B[3] = |B|
        w     - workspace
*/

int
mfield_eval(const double t, const double r, const double theta,
            const double phi, double B[4], mfield_workspace *w)
{
  int s = 0;
  gsl_vector *c = w->c_copy;

  /* convert coefficients to physical time units */
  mfield_coeffs(1, w->c, c, w);

  s = mfield_eval_g(t, r, theta, phi, c, B, w);

  return s;
} /* mfield_eval() */

/*
mfield_eval_dBdt()
  Evaluate magnetic field model time derivative at given point

Inputs: t     - timestamp (CDF_EPOCH)
        r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - geocentric longitude (radians)
        dBdt  - (output) magnetic field
                dBdt[0] = d/dt B_x
                dBdt[1] = d/dt B_y
                dBdt[2] = d/dt B_z
                dBdt[3] = |dBdt|
        w     - workspace
*/

int
mfield_eval_dBdt(const double t, const double r, const double theta,
                 const double phi, double dBdt[4], mfield_workspace *w)
{
  int s = 0;
  gsl_vector *c = w->c_copy;
#if MFIELD_FIT_SECVAR
  gsl_vector_view dg = gsl_vector_subvector(c, w->sv_offset, w->nnm_sv);
#else
  gsl_vector_view dg;
#endif
#if MFIELD_FIT_SECACC
  gsl_vector_view ddg = gsl_vector_subvector(c, w->sa_offset, w->nnm_sa);
#else
  gsl_vector_view ddg;
#endif

  /* convert coefficients to physical time units */
  mfield_coeffs(1, w->c, c, w);

  s = mfield_eval_dgdt(t, r, theta, phi, &dg.vector, &ddg.vector, dBdt, w);

  return s;
} /* mfield_eval_dBdt() */

/*
mfield_eval_g()
  Evaluate magnetic field model for given coefficients

Inputs: t     - timestamp (CDF_EPOCH)
        r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        c     - field coefficients (nT,nT/year,nT/year^2)
        B     - (output) magnetic field
                B[0] = B_x
                B[1] = B_y
                B[2] = B_z
                B[3] = |B|
        w     - workspace
*/

static int
mfield_eval_g(const double t, const double r, const double theta, const double phi,
              const gsl_vector *c, double B[4], mfield_workspace *w)
{
  int s = 0;
  size_t n;
  int m;

  /* convert to years and subtract epoch */
  const double t1 = satdata_epoch2year(t) - w->epoch; /* SV term (years) */
  const double t2 = 0.5 * t1 * t1;                    /* SA term (years^2) */

  s += mfield_green(r, theta, phi, w);

  B[0] = B[1] = B[2] = 0.0;

  for (n = 1; n <= w->nmax_mf; ++n)
    {
      int ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);
          double g = mfield_get_mf(c, cidx, w);
          double dg = mfield_get_sv(c, cidx, w);
          double ddg = mfield_get_sa(c, cidx, w);
          double gnm = g + dg * t1 + ddg * t2;

          B[0] += gnm * w->dX[cidx];
          B[1] += gnm * w->dY[cidx];
          B[2] += gnm * w->dZ[cidx];
        }
    }

  B[3] = gsl_hypot3(B[0], B[1], B[2]);

  return s;
} /* mfield_eval_g() */

/*
mfield_eval_dgdt()
  Evaluate magnetic field model time derivative for given coefficients

Inputs: t     - timestamp (CDF_EPOCH)
        r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        dg    - secular variation coefficients (nT/year)
        ddg   - secular acceleration coefficients (nT/year^2)
        dBdt  - (output) magnetic field
                dBdt[0] = d/dt B_x
                dBdt[1] = d/dt B_y
                dBdt[2] = d/dt B_z
                dBdt[3] = |dBdt|
        w     - workspace
*/

int
mfield_eval_dgdt(const double t, const double r, const double theta,
                 const double phi, const gsl_vector *dg,
                 const gsl_vector *ddg, double dBdt[4],
                 mfield_workspace *w)
{
  int s = 0;
  size_t n;
  int m;

  /* convert to years and subtract epoch */
  /*const double t1 = satdata_epoch2year(t) - w->epoch;*/
  const double t1 = t - w->epoch;

  s += mfield_green(r, theta, phi, w);

  dBdt[0] = dBdt[1] = dBdt[2] = 0.0;

  for (n = 1; n <= w->nmax_mf; ++n)
    {
      int ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);
          double dg0 = gsl_vector_get(dg, cidx);
          double ddg0 = gsl_vector_get(ddg, cidx);
          double gnm = dg0 + ddg0 * t1;

          dBdt[0] += gnm * w->dX[cidx];
          dBdt[1] += gnm * w->dY[cidx];
          dBdt[2] += gnm * w->dZ[cidx];
        }
    }

  dBdt[3] = gsl_hypot3(dBdt[0], dBdt[1], dBdt[2]);

  return s;
} /* mfield_eval_dgdt() */

/*
mfield_eval_ext()
  Evaluate external magnetic field model at given point

Inputs: t     - timestamp (CDF_EPOCH)
        r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - geocentric longitude (radians)
        B     - (output) magnetic field
                B[0] = B_x
                B[1] = B_y
                B[2] = B_z
                B[3] = |B|
        w     - workspace
*/

int
mfield_eval_ext(const double t, const double r, const double theta, const double phi,
                double B[4], mfield_workspace *w)
{
  int s = 0;

#if MFIELD_FIT_EXTFIELD

  size_t extidx = mfield_extidx(t, w);
  double extcoeff = gsl_vector_get(w->c, extidx);
  size_t i;

  s = mfield_nonlinear_model_ext(r, theta, phi, w->c, B, w);

  /*
   * if there was not enough data on a particular day, the computed
   * coefficient could be wildly wrong
   */
  if (fabs(extcoeff) > 50.0)
    extcoeff = 0.0;

  for (i = 0; i < 3; ++i)
    B[i] *= extcoeff;

#else
  
  B[0] = B[1] = B[2] = 0.0;

#endif

  B[3] = gsl_hypot3(B[0], B[1], B[2]);

  return s;
} /* mfield_eval_ext() */

/*
mfield_eval_ext_coeff()
  Evaluate external magnetic field model at given point

Inputs: r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - geocentric longitude (radians)
        extcoeff - external coefficient (nT)
        B     - (output) magnetic field
                B[0] = B_x
                B[1] = B_y
                B[2] = B_z
                B[3] = |B|
        w     - workspace
*/

int
mfield_eval_ext_coeff(const double r, const double theta, const double phi,
                      const double extcoeff, double B[4], mfield_workspace *w)
{
  int s = 0;
  size_t i;

  s = mfield_nonlinear_model_ext(r, theta, phi, w->c, B, w);

  for (i = 0; i < 3; ++i)
    B[i] *= extcoeff;

  B[3] = gsl_hypot3(B[0], B[1], B[2]);

  return s;
} /* mfield_eval_ext_coeff() */

/*
mfield_eval_g_ext()
  Evaluate external magnetic field model at given point

Inputs: t     - timestamp (CDF_EPOCH)
        r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - geocentric longitude (radians)
        B     - (output) magnetic field
                B[0] = B_x
                B[1] = B_y
                B[2] = B_z
                B[3] = |B|
        w     - workspace
*/

int
mfield_eval_g_ext(const double t, const double r, const double theta, const double phi,
                  const double E_st, const double I_st,
                  const gsl_vector *g, const gsl_vector *dg,
                  double B[4], mfield_workspace *w)
{
  int s = 0;
  int m;
  double sint = sin(theta);
  double rterm;
  double dt = satdata_epoch2year(t) - 2012.0; /* HDGM 2013 epoch */

  /* no radial term for n = 1 external, only internal */
  rterm = pow(w->R / r, 3.0);

  /* compute associated legendres */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, 1, cos(theta), w->Plm, w->dPlm);

  B[0] = B[1] = B[2] = 0.0;

  for (m = 0; m <= 1; ++m)
    {
      size_t cidx = mfield_coeff_nmidx(1, m);
      size_t pidx = gsl_sf_legendre_array_index(1, m);
      double g1m = gsl_vector_get(g, cidx) + dt * gsl_vector_get(dg, cidx);
      double h1m = 0.0;

      if (m != 0)
        {
          cidx = mfield_coeff_nmidx(1, -m);
          h1m = gsl_vector_get(g, cidx) + dt * gsl_vector_get(dg, cidx);
        }

      /* external contribution */
      B[0] += E_st * (g1m * cos(m * phi) + h1m * sin(m * phi)) * w->dPlm[pidx];
      B[1] += E_st * m / sint * (g1m * sin(m * phi) - h1m * cos(m * phi)) * w->Plm[pidx];
      B[2] += E_st * (g1m * cos(m * phi) + h1m * sin(m * phi)) * w->Plm[pidx];

      /* internal contribution */
      B[0] += I_st * rterm * (g1m * cos(m * phi) + h1m * sin(m * phi)) * w->dPlm[pidx];
      B[1] += I_st * rterm * m / sint * (g1m * sin(m * phi) - h1m * cos(m * phi)) * w->Plm[pidx];
      B[2] -= I_st * 2.0 * rterm * (g1m * cos(m * phi) + h1m * sin(m * phi)) * w->Plm[pidx];
    }

  B[3] = gsl_hypot3(B[0], B[1], B[2]);

  return s;
} /* mfield_eval_g_ext() */

/*
mfield_eval_static()
  Evaluate magnetic field model

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        g     - spherical harmonic coefficients in units of
                nT, nT/year, nT/year^2
        B     - (output) magnetic field in nT
                B[0] = B_x
                B[1] = B_y
                B[2] = B_z
                B[3] = |B|
        w     - workspace
*/

int
mfield_eval_static(const double r, const double theta, const double phi,
                   const gsl_vector *g, double B[4], mfield_workspace *w)
{
  int s = 0;
  size_t n;
  int m;

  s += mfield_green(r, theta, phi, w);

  B[0] = B[1] = B[2] = 0.0;

  for (n = 1; n <= w->nmax_mf; ++n)
    {
      int ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);
          double gnm = gsl_vector_get(g, cidx);

          B[0] += gnm * w->dX[cidx];
          B[1] += gnm * w->dY[cidx];
          B[2] += gnm * w->dZ[cidx];
        }
    }

  B[3] = gsl_hypot3(B[0], B[1], B[2]);

  return s;
} /* mfield_eval_static() */

double
mfield_spectrum(const size_t n, const mfield_workspace *w)
{
  int m, ni = (int) n;
  double sum = 0.0;
  gsl_vector *c = w->c_copy;

  /* undo time scaling of coefficients */
  mfield_coeffs(1, w->c, c, w);

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = mfield_coeff_nmidx(n, m);
      double gnm = gsl_vector_get(c, cidx);

      sum += gnm * gnm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
} /* mfield_spectrum() */

/* SV spectrum */
double
mfield_spectrum_sv(const size_t n, const mfield_workspace *w)
{
  int m, ni = (int) n;
  double sum = 0.0;
  gsl_vector *c = w->c_copy;

  /* undo time scaling of coefficients */
  mfield_coeffs(1, w->c, c, w);

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = mfield_coeff_nmidx(n, m);
      double dgnm = mfield_get_sv(c, cidx, w);

      sum += dgnm * dgnm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
} /* mfield_spectrum_sv() */

/* SA spectrum */
double
mfield_spectrum_sa(const size_t n, const mfield_workspace *w)
{
  int m, ni = (int) n;
  double sum = 0.0;
  gsl_vector *c = w->c_copy;

  /* undo time scaling of coefficients */
  mfield_coeffs(1, w->c, c, w);

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = mfield_coeff_nmidx(n, m);
      double ddgnm = mfield_get_sa(c, cidx, w);

      sum += ddgnm * ddgnm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
} /* mfield_spectrum_sa() */

int
mfield_write(const char *filename, mfield_workspace *w)
{
  int s = 0;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "mfield_write: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  fwrite(&(w->params), sizeof(mfield_parameters), 1, fp);
  fwrite(&(w->t_mu), sizeof(double), 1, fp);
  fwrite(&(w->t_sigma), sizeof(double), 1, fp);
  fwrite(&(w->t0_data), sizeof(double), 1, fp);

  /*
   * only write internal coefficients since when we later read
   * the file we won't be able to recalculate w->neuler
   */
  fwrite(w->c->data, sizeof(double), w->p_int, fp);

  gsl_matrix_fwrite(fp, w->covar);

  fclose(fp);

  return s;
} /* mfield_write() */

mfield_workspace *
mfield_read(const char *filename)
{
  mfield_workspace *w;
  mfield_parameters params;
  FILE *fp;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "mfield_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  fread(&params, sizeof(mfield_parameters), 1, fp);

  params.mfield_data_p = NULL;

  w = mfield_alloc(&params);

  fread(&(w->t_mu), sizeof(double), 1, fp);
  fread(&(w->t_sigma), sizeof(double), 1, fp);
  fread(&(w->t0_data), sizeof(double), 1, fp);
  fread(w->c->data, sizeof(double), w->p_int, fp);
  gsl_matrix_fread(fp, w->covar);

  fclose(fp);

  return w;
} /* mfield_read() */

/*
mfield_write_ascii()
  Write ascii coefficient file
*/

int
mfield_write_ascii(const char *filename, const double epoch,
                   const int write_delta, mfield_workspace *w)
{
  int s = 0;
  FILE *fp;
  size_t n, i;
  gsl_vector_view v = gsl_matrix_diagonal(w->covar);
  gsl_vector_view cov = gsl_vector_subvector(w->c_copy, 0, w->p_int);

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "mfield_write_ascii: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  /* convert coefficients to physical time units */
  mfield_coeffs(1, w->c, w->c, w);

  gsl_vector_memcpy(&cov.vector, &v.vector);
  for (i = 0; i < w->p_int; ++i)
    {
      double ci = gsl_vector_get(&cov.vector, i);
      gsl_vector_set(&cov.vector, i, sqrt(ci));
    }

  /* convert covariance vector to physical units */
  mfield_coeffs(1, w->c_copy, w->c_copy, w);

  /* print header information */
  fprintf(fp, "%% Magnetic field model coefficients\n");
  fprintf(fp, "%% nmax:  %zu\n", w->nmax_mf);
  fprintf(fp, "%% epoch: %.4f\n", epoch);
  fprintf(fp, "%% radius: %.1f\n", w->R);

  if (write_delta)
    {
      fprintf(fp, "%% %3s %5s %20s %20s %20s %20s %20s %20s\n",
              "n",
              "m",
              "MF gnm (nT)",
              "MF dgnm (nT)",
              "SV gnm (nT/year)",
              "SV dgnm (nT/year)",
              "SA gnm (nT/year^2)",
              "SA dgnm (nT/year^2)");
    }
  else
    {
      fprintf(fp, "%% %3s %5s %20s %20s %20s\n",
              "n",
              "m",
              "MF gnm (nT)",
              "SV gnm (nT/year)",
              "SA gnm (nT/year^2)");
    }

  for (n = 1; n <= w->nmax_mf; ++n)
    {
      int m, ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);
          double gnm = mfield_get_mf(w->c, cidx, w); 
          double dgnm = mfield_get_sv(w->c, cidx, w);
          double ddgnm = mfield_get_sa(w->c, cidx, w); 

          if (write_delta)
            {
              double gnm_delta = mfield_get_mf(&cov.vector, cidx, w); 
              double dgnm_delta = mfield_get_sv(&cov.vector, cidx, w); 
              double ddgnm_delta = mfield_get_sa(&cov.vector, cidx, w); 

              fprintf(fp, "%5zu %5d %20.4f %20.4f %20.4f %20.4f %20.4f %20.4f\n",
                      n,
                      m,
                      gnm,
                      gnm_delta,
                      dgnm,
                      dgnm_delta,
                      ddgnm,
                      ddgnm_delta);
            }
          else
            {
              fprintf(fp, "%5zu %5d %20.4f %20.4f %20.4f\n",
                      n,
                      m,
                      gnm,
                      dgnm,
                      ddgnm);
            }
        }
    }

  fclose(fp);

  /* convert coefficients to dimensionless time units */
  mfield_coeffs(-1, w->c, w->c, w);

  return s;
} /* mfield_write_ascii() */

/*
mfield_new_epoch()
  Extrapolate model coefficients to new epoch. The
old Taylor series is:

gnm(t) = gnm + (t - t0) dgnm + 1/2 (t - t0)^2 ddgnm

We want a taylor series around a new epoch:

gnm(t) = new_gnm + (t - new_t0) new_dgnm + 1/2 (t - new_t0)^2 new_ddgnm

The new taylor coefficients are given as:

new_gnm   = gnm + a dgnm + 1/2 a^2 ddgnm
new_dgnm  = dgnm + a ddgnm
new_ddgnm = ddgnm

with: a = new_t0 - t0
*/

int
mfield_new_epoch(const double new_epoch, mfield_workspace *w)
{
  int s = 0;
  size_t n;
  gsl_vector *c = w->c_copy;
  double t1 = new_epoch - w->epoch;
  double t2 = 0.5 * t1 * t1;

  /* convert coefficients to physical time units */
  mfield_coeffs(1, w->c, c, w);

  for (n = 1; n <= w->nmax_mf; ++n)
    {
      int m, ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);
          double gnm = mfield_get_mf(c, cidx, w); 
          double dgnm = mfield_get_sv(c, cidx, w);
          double ddgnm = mfield_get_sa(c, cidx, w);

          mfield_set_mf(c, cidx, gnm + t1 * dgnm + t2 * ddgnm, w);
          mfield_set_sv(c, cidx, dgnm + t1 * ddgnm, w);
          mfield_set_sa(c, cidx, ddgnm, w);
        }
    }

  /* convert back to dimensionless units */
  mfield_coeffs(-1, c, w->c, w);

  return s;
} /* mfield_new_epoch() */

/*******************************************************
 *      INTERNAL ROUTINES                              *
 *******************************************************/

/*
mfield_green()
  Compute Green's functions for X,Y,Z spherical harmonic expansion. These
are simply the basis functions multiplying the g_{nm} and h_{nm} coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Plm
w->dPlm
w->sinmphi
w->cosmphi

2) The output Green's functions are stored in w->dX, w->dY, w->dZ
*/

static int
mfield_green(const double r, const double theta, const double phi, mfield_workspace *w)
{
  int s = 0;
  size_t n;
  int m;
  const double sint = sin(theta);
  const double cost = cos(theta);
  double ratio = w->R / r;
  double term = ratio * ratio;     /* (a/r)^{n+2} */

  /* precompute cos(m phi) and sin(m phi) */
  for (n = 0; n <= w->nmax_mf; ++n)
    {
      w->cosmphi[n] = cos(n * phi);
      w->sinmphi[n] = sin(n * phi);
    }

  /* compute associated legendres */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, w->nmax_mf, cost,
                                  w->Plm, w->dPlm);

  for (n = 1; n <= w->nmax_mf; ++n)
    {
      int ni = (int) n;

      /* (a/r)^{n+2} */
      term *= ratio;

      for (m = -ni; m <= ni; ++m)
        {
          int mabs = abs(m);
          size_t cidx = mfield_coeff_nmidx(n, m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);

          if (m < 0)
            {
              /* h_{nm} */
              w->dX[cidx] = term * w->sinmphi[mabs] * w->dPlm[pidx];
              w->dY[cidx] = -term / sint * mabs * w->cosmphi[mabs] * w->Plm[pidx];
              w->dZ[cidx] = -(n + 1.0) * term * w->sinmphi[mabs] * w->Plm[pidx];
            }
          else
            {
              /* g_{nm} */
              w->dX[cidx] = term * w->cosmphi[mabs] * w->dPlm[pidx];
              w->dY[cidx] = term / sint * mabs * w->sinmphi[mabs] * w->Plm[pidx];
              w->dZ[cidx] = -(n + 1.0) * term * w->cosmphi[mabs] * w->Plm[pidx];
            }
        }
    }

  return s;
} /* mfield_green() */

/*
mfield_extidx()
  Return external field coefficient index for a given time

Inputs: t - CDF_EPOCH (needed for doy)
        w - workspace

Return: external coefficient corresponding to doy
*/

size_t
mfield_extidx(const double t, const mfield_workspace *w)
{
  double fday = satdata_epoch2fday(t);
  double fday0 = satdata_epoch2fday(w->t0_data);
  int daynum = (int) (fday - fday0);

  return (w->ext_offset + daynum);
}

/*
mfield_coeff_nmidx()
  This function returns a unique index in [0,w->p-1] corresponding
to a given (l,m) pair. The array will look like:

[(1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) ...]

(the (0,0) coefficient is not solved for)

Inputs: n - SH degree (> 0)
        m - SH order (-l <= m <= l)

Return: index in [0,nnm-1]
*/

size_t
mfield_coeff_nmidx(const size_t n, const int m)
{
  size_t base = n * n; /* index of block for this n */
  int offset = m + n;  /* offset within block for this m */
  size_t nmidx;

  if (n == 0)
    {
      fprintf(stderr, "mfield_coeff_nmidx: error: n = 0\n");
      return 0;
    }

  nmidx = base + offset;

  /* subtract 1 to exclude (0,0) coefficient */
  return nmidx - 1;
} /* mfield_coeff_nmidx() */

double
mfield_get_mf(const gsl_vector *c, const size_t idx,
              const mfield_workspace *w)
{
  return gsl_vector_get(c, idx);
}

double
mfield_get_sv(const gsl_vector *c, const size_t idx,
              const mfield_workspace *w)
{
  if (idx < w->nnm_sv)
    return gsl_vector_get(c, idx + w->sv_offset);
  else
    return 0.0;
}

double
mfield_get_sa(const gsl_vector *c, const size_t idx,
              const mfield_workspace *w)
{
  if (idx < w->nnm_sa)
    return gsl_vector_get(c, idx + w->sa_offset);
  else
    return 0.0;
}

int
mfield_set_mf(gsl_vector *c, const size_t idx,
              const double x,
              const mfield_workspace *w)
{
  gsl_vector_set(c, idx, x);
  return GSL_SUCCESS;
}

int
mfield_set_sv(gsl_vector *c, const size_t idx,
              const double x,
              const mfield_workspace *w)
{
#if MFIELD_FIT_SECVAR
  /* check idx in case nmax_sv is less than nmax_mf */
  if (idx < w->nnm_sv)
    gsl_vector_set(c, idx + w->sv_offset, x);
#endif
  return GSL_SUCCESS;
}

int
mfield_set_sa(gsl_vector *c, const size_t idx,
              const double x,
              const mfield_workspace *w)
{
#if MFIELD_FIT_SECACC
  /* check idx in case nmax_sa is less than nmax_mf */
  if (idx < w->nnm_sa)
    gsl_vector_set(c, idx + w->sa_offset, x);
#endif
  return GSL_SUCCESS;
}

static int
mfield_print_residuals(const int header, FILE *fp,
                       const gsl_vector *r,
                       const gsl_vector *w_rob,
                       const gsl_vector *w_spatial)
{
  size_t i;

  if (!fp)
    return GSL_SUCCESS;

  if (header)
    {
      i = 1;
      fprintf(fp, "# Field %zu: residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: robust weight\n", i++);
      fprintf(fp, "# Field %zu: spatial weight\n", i++);
      fprintf(fp, "# Field %zu: final weight (robust x spatial)\n", i++);
      return GSL_SUCCESS;
    }

  for (i = 0; i < r->size; ++i)
    {
      double robi = gsl_vector_get(w_rob, i);
      double spati = gsl_vector_get(w_spatial, i);

      fprintf(fp, "%.5e %.5e %.5e %.5e\n",
              gsl_vector_get(r, i),
              robi,
              spati,
              robi * spati);
    }

  return GSL_SUCCESS;
} /* mfield_print_residuals() */

static int
mfield_update_histogram(const gsl_vector *r, gsl_histogram *h)
{
  const size_t n = r->size;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double ri = gsl_vector_get(r, i);
      gsl_histogram_increment(h, ri);
    }

  return GSL_SUCCESS;
} /* mfield_update_histogram() */

/*
mfield_print_histogram()
  Print histogram and normalize to unit area
*/

static int
mfield_print_histogram(FILE *fp, gsl_histogram *h)
{
  const size_t n = gsl_histogram_bins(h);
  const double sum = gsl_histogram_sum(h);
  size_t i;

  fprintf(fp, "# Histogram variable mean: %f\n", gsl_histogram_mean(h));
  fprintf(fp, "# Histogram variable sigma: %f\n", gsl_histogram_sigma(h));

  for (i = 0; i < n; ++i)
    {
      double hi = gsl_histogram_get(h, i);
      double lower, upper, width;

      gsl_histogram_get_range(h, i, &lower, &upper);
      width = upper - lower;

      fprintf(fp, "%g %.12e\n",
              0.5*(lower + upper),
              hi / (width * sum));
    }

  return GSL_SUCCESS;
} /* mfield_print_histogram() */
