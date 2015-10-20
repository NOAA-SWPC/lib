/*
 * mag_sqfilt.c
 *
 * Fit spherical harmonic model to satellite scalar residuals to remove
 * Sq and external fields
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_test.h>

#include "coord.h"
#include "common.h"
#include "interp.h"
#include "mag.h"

#include "Gdef.h"

static int sqfilt_linear_init(mag_workspace *mag_p, mag_sqfilt_workspace *w);
static int sqfilt_linear_matrix_row(const double r, const double thetaq,
                                    const double phi,
                                    const double b_int[3],
                                    gsl_vector *v, mag_workspace *w);
static int sqfilt_calc_F2(const gsl_vector *c, mag_workspace *w);
static size_t sqfilt_nmidx(const size_t type, const size_t n, const int m,
                           const mag_sqfilt_workspace *w);

mag_sqfilt_workspace *
mag_sqfilt_alloc(const size_t nmax_int, const size_t mmax_int,
                 const size_t nmax_ext, const size_t mmax_ext)
{
  const size_t ndata = 10000;
  size_t l;
  mag_sqfilt_workspace *w;

  w = calloc(1, sizeof(mag_sqfilt_workspace));
  if (!w)
    return 0;

  w->nmax_int = nmax_int;
  w->mmax_int = mmax_int;
  w->nmax_ext = nmax_ext;
  w->mmax_ext = mmax_ext;

  w->base_int = calloc(1, (w->nmax_int + 1) * sizeof(size_t));
  w->base_ext = calloc(1, (w->nmax_ext + 1) * sizeof(size_t));

  /*
   * precompute base offsets for internal (n,m) indexing, and count total
   * internal coefficients
   */
  w->p_int = 0;
  for (l = 1; l <= w->nmax_int; ++l)
    {
      int ns = (int) GSL_MIN(l, w->mmax_int);

      w->base_int[l] = w->p_int;
      w->p_int += 2*ns + 1;
    }

  /*
   * precompute base offsets for external (n,m) indexing, and count total
   * external coefficients
   */
  w->p_ext = 0;
  for (l = 1; l <= w->nmax_ext; ++l)
    {
      int ns = (int) GSL_MIN(l, w->mmax_ext);

      w->base_ext[l] = w->p_ext;
      w->p_ext += 2*ns + 1;
    }

  /* total model parameters */
  w->p = w->p_int + w->p_ext;

  w->int_offset = 0;
  w->ext_offset = w->int_offset + w->p_int;

  w->X = gsl_matrix_alloc(ndata, w->p);
  w->cov = gsl_matrix_alloc(w->p, w->p);
  w->c = gsl_vector_alloc(w->p);
  w->rhs = gsl_vector_alloc(ndata);
  w->work_p = gsl_vector_alloc(w->p);

  w->ntot = ndata;
  w->n = 0;

  w->nreg = 200;
  w->rho = gsl_vector_alloc(w->nreg);
  w->eta = gsl_vector_alloc(w->nreg);
  w->reg_param = gsl_vector_alloc(w->nreg);

  return w;
} /* mag_sqfilt_alloc() */

void
mag_sqfilt_free(mag_sqfilt_workspace *w)
{
  if (w->base_int)
    free(w->base_int);

  if (w->base_ext)
    free(w->base_ext);

  if (w->X)
    gsl_matrix_free(w->X);

  if (w->cov)
    gsl_matrix_free(w->cov);

  if (w->c)
    gsl_vector_free(w->c);

  if (w->rhs)
    gsl_vector_free(w->rhs);

  if (w->work_p)
    gsl_vector_free(w->work_p);

  if (w->multifit_workspace_p)
    gsl_multifit_linear_free(w->multifit_workspace_p);

  if (w->rho)
    gsl_vector_free(w->rho);

  if (w->eta)
    gsl_vector_free(w->eta);

  if (w->reg_param)
    gsl_vector_free(w->reg_param);

  free(w);
}

/*
mag_sqfilt()
  Fit a simple spherical harmonic model to a single satellite
track of scalar F^(1) residual data, excluding the EEJ at low-latitudes.

Notes:
1) On output, w->rnorm and w->snorm contain the
residual norm and solution norm

2) On output, mag_p->track.F2 contains the F^(2) residuals and
mag_p->track.Sq_model contains the b . (M + K) Sq model
*/

int
mag_sqfilt(mag_workspace *mag_p, mag_sqfilt_workspace *w)
{
  int s;
  gsl_vector_view bv;
  gsl_matrix_view Xv;
  double h; /* final damping parameter */

  s = sqfilt_linear_init(mag_p, w);
  if (s)
    return s;

  Xv = gsl_matrix_submatrix(w->X, 0, 0, w->n, w->p);
  bv = gsl_vector_subvector(w->rhs, 0, w->n);

  /* compute SVD of LS matrix */
  s = gsl_multifit_linear_svd(&Xv.matrix, w->multifit_workspace_p);
  if (s)
    return s;

  /* compute L-curve */
  s = gsl_multifit_linear_ridge_lcurve(&bv.vector, w->reg_param, w->rho, w->eta, w->multifit_workspace_p);
  if (s)
    return s;

  /* compute L-curve corner
   * 2015-10-01: found that the corner2 method has worse performance during storms
   */
  s = gsl_multifit_linear_ridge_lcorner(w->rho, w->eta, &(w->reg_idx));
  /*s = gsl_multifit_linear_ridge_lcorner2(w->reg_param, w->eta, &(w->reg_idx));*/
  if (s)
    return s;

  h = gsl_vector_get(w->reg_param, w->reg_idx);

  /* perform damped least squares fit */
  s = gsl_multifit_linear_ridge_solve(h, &bv.vector, w->c,
                                      w->cov, &(w->rnorm), &(w->snorm), w->multifit_workspace_p);
  if (s)
    return s;

  fprintf(stderr, "mag_sqfilt: final regularization factor = %g\n", h);
  fprintf(stderr, "mag_sqfilt: multifit residual norm = %.6e\n", w->rnorm);
  fprintf(stderr, "mag_sqfilt: multifit solution norm = %.6e\n", w->snorm);

  /* calculate F^(2) residuals (eq 9 of paper) */
  sqfilt_calc_F2(w->c, mag_p);

  return s;
} /* mag_sqfilt() */

/****************************************************
 *       INTERNAL ROUTINES                          *
 ****************************************************/

/*
sqfilt_linear_init()
  Fit a simple spherical harmonic model to a single satellite
track of scalar F^(1) residual data, excluding the EEJ at low-latitudes.
This function initializes the LS design matrix and rhs vector

Inputs: mag_p - mag workspace
        w     - sqfilt workspace

Notes:
1) mag_p->track must be initialized on input with a complete satellite track

2) w->n contains the number of data points to be fitted
*/

static int
sqfilt_linear_init(mag_workspace *mag_p, mag_sqfilt_workspace *w)
{
  int s = 0;

  /*
   * exclude points below this QD latitude from fit, except a few to
   * stabilize filter
   */
  const double exclude_lat = 12.0;

  const size_t ntot = mag_p->track.n;
  mag_track *track = &(mag_p->track);
  size_t i;
  size_t ndata = 0;
  double qd_min, qd_max, F1_min, F1_max;
  int found_lower = 0, found_upper = 0;

  for (i = 1; i < ntot; ++i)
    {
      double qdprev = track->qdlat[i - 1];
      double qdcur = track->qdlat[i];

      if (fabs(qdprev) > exclude_lat && fabs(qdcur) <= exclude_lat)
        {
          assert(found_lower == 0);

          found_lower = 1;
          qd_min = qdcur;
          F1_min = track->F1[i];
        }
      else if (fabs(qdprev) < exclude_lat && fabs(qdcur) >= exclude_lat)
        {
          assert(found_upper == 0);

          found_upper = 1;
          qd_max = qdcur;
          F1_max = track->F1[i];
        }
    }

  if (!found_lower || !found_upper)
    {
      fprintf(stderr, "sqfilt_linear_init: unable to find linear interpolation points\n");
      return -1;
    }

  /* build least squares matrix and RHS */
  for (i = 0; i < ntot; ++i)
    {
      double r = track->r[i];
      double thetaq = track->thetaq[i]; /* use QD colatitude */
      double phi = track->phi[i];
      double rhsval = track->F1[i];
      double b_int[3];
      gsl_vector_view v;

      /* exclude points in the equatorial electrojet region */
      if (fabs(track->qdlat[i]) < exclude_lat)
        {
          /* keep every 10th point in the equatorial region */
          if (i % 10 == 0)
            {
              /* interpolate linearly across EEJ region */
              rhsval = interp1d(qd_min, qd_max, F1_min, F1_max,
                                track->qdlat[i]);
            }
          else
            continue;
        }

      /* compute unit vector in internal field direction */
      b_int[0] = track->Bx_int[i] / track->F_int[i];
      b_int[1] = track->By_int[i] / track->F_int[i];
      b_int[2] = track->Bz_int[i] / track->F_int[i];

      /* add green basis functions to current row of matrix */
      v = gsl_matrix_row(w->X, ndata);
      sqfilt_linear_matrix_row(r, thetaq, phi, b_int, &v.vector, mag_p);

      /* set RHS value */
      gsl_vector_set(w->rhs, ndata, rhsval);

      ++ndata;
    }

  if (ndata != w->n)
    {
      if (w->multifit_workspace_p)
        gsl_multifit_linear_free(w->multifit_workspace_p);

      w->multifit_workspace_p = gsl_multifit_linear_alloc(ndata, w->p);
    }

  w->n = ndata;

  return s;
} /* sqfilt_linear_init() */

/*
sqfilt_linear_matrix_row()
  Add a row to the least squares matrix

Inputs: r      - radius (km)
        thetaq - QD colatitude (radians)
        phi    - geographic longitude (radians)
        b_int  - unit vector for internal field
        v      - row of LS matrix
        w      - workspace
*/

static int
sqfilt_linear_matrix_row(const double r, const double thetaq,
                         const double phi, const double b_int[3],
                         gsl_vector *v, mag_workspace *w)
{
  int s = 0;
  const double R = 6371.2;
  green_workspace *green_p = w->green_workspace_p;
  mag_sqfilt_workspace *sqfilt_p = w->sqfilt_workspace_p;
  double phi_sm, theta_sm, lat_sm;
  double lat = M_PI / 2.0 - thetaq;
  time_t unix_time = satdata_epoch2timet(w->track.t_eq);
  double fday = time2fday(unix_time);
  size_t i, n;

  /* compute basis functions for M(r,thetaq,phi) */
  green_calc(r, thetaq, phi, R, green_p);

  /* compute basis functions for K(r,thetaq_SM,phi_SM) */
  trans(GEO2SM, fday, phi, lat, &phi_sm, &lat_sm);
  theta_sm = M_PI / 2.0 - lat_sm;

  green_calc_ext(r, theta_sm, phi_sm, R, green_p);

  for (i = 0; i < green_p->nnm; ++i)
    {
      double phi_ss, lat_ss;
      double X, Y, Z;

      trans_vec(SM2GEO, fday, phi_sm, lat_sm, green_p->dX_ext[i],
                green_p->dY_ext[i], green_p->dZ_ext[i], &phi_ss,
                &lat_ss, &X, &Y, &Z);

      green_p->dX_ext[i] = X;
      green_p->dY_ext[i] = Y;
      green_p->dZ_ext[i] = Z;
    }

  /* internal coefficients */
  for (n = 1; n <= sqfilt_p->nmax_int; ++n)
    {
      int M = (int) GSL_MIN(n, sqfilt_p->mmax_int);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t idx = sqfilt_nmidx(0, n, m, sqfilt_p);
          size_t gidx = green_nmidx(n, m);
          double val = b_int[0] * green_p->dX[gidx] +
                       b_int[1] * green_p->dY[gidx] +
                       b_int[2] * green_p->dZ[gidx];

          gsl_vector_set(v, idx, val);
        }
    }

  /* external coefficients */
  for (n = 1; n <= sqfilt_p->nmax_ext; ++n)
    {
      int M = (int) GSL_MIN(n, sqfilt_p->mmax_ext);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t idx = sqfilt_nmidx(1, n, m, sqfilt_p);
          size_t gidx = green_nmidx(n, m);
          double val = b_int[0] * green_p->dX_ext[gidx] +
                       b_int[1] * green_p->dY_ext[gidx] +
                       b_int[2] * green_p->dZ_ext[gidx];

          gsl_vector_set(v, idx, val);
        }
    }

  return s;
} /* sqfilt_linear_matrix_row() */

/*
sqfilt_calc_F2()
  Calculate F^(2) residual (eq. 9 of paper) by subtracting
the model b . (M + K)

Inputs: c - model coefficients
        w - workspace

Notes:
1) On output, w->track.F2 is modified to contain F^(2) residuals
*/

static int
sqfilt_calc_F2(const gsl_vector *c, mag_workspace *w)
{
  mag_sqfilt_workspace *sqfilt_p = w->sqfilt_workspace_p;
  size_t i;
  mag_track *track = &(w->track);
  gsl_vector *v = sqfilt_p->work_p;

  gsl_vector_const_view c_int = gsl_vector_const_subvector(c, sqfilt_p->int_offset, sqfilt_p->p_int);
  gsl_vector_const_view c_ext = gsl_vector_const_subvector(c, sqfilt_p->ext_offset, sqfilt_p->p_ext);
  gsl_vector_const_view v_int = gsl_vector_const_subvector(v, sqfilt_p->int_offset, sqfilt_p->p_int);
  gsl_vector_const_view v_ext = gsl_vector_const_subvector(v, sqfilt_p->ext_offset, sqfilt_p->p_ext);

  for (i = 0; i < track->n; ++i)
    {
      double r = track->r[i];
      double thetaq = track->thetaq[i]; /* use QD colatitude */
      double phi = track->phi[i];
      double b_int[3];
      double val; /* Sq model value b . (M + K) */

      /* compute unit vector in internal field direction */
      b_int[0] = track->Bx_int[i] / track->F_int[i];
      b_int[1] = track->By_int[i] / track->F_int[i];
      b_int[2] = track->Bz_int[i] / track->F_int[i];

      /* compute basis functions */
      sqfilt_linear_matrix_row(r, thetaq, phi, b_int, v, w);

      /* compute internal model = b . M */
      gsl_blas_ddot(&v_int.vector, &c_int.vector, &val);
      track->Sq_int[i] = val;

      /* compute external model = b . K */
      gsl_blas_ddot(&v_ext.vector, &c_ext.vector, &val);
      track->Sq_ext[i] = val;

      track->F2[i] = track->F1[i] - (track->Sq_int[i] + track->Sq_ext[i]);
    }

  return GSL_SUCCESS;
} /* sqfilt_calc_F2() */

/*
sqfilt_nmidx()
  This function returns a unique index in [0,p-1] corresponding
to a given (n,m) pair. The array will look like:

[(1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) ...]

(the (0,0) coefficient is not solved for)

Things get a little more tricky when mmax != nmax, so the
base indices of each 'n' block are precomputed in _alloc and
stored for easy reference here. The offset of (n,m) in a given
'n' block is:

offset = m + min(n,mmax)

which defaults to the standard m + n for the case mmax = nmax

Inputs: type - 0 for internal index, 1 for external index
        n    - SH degree (> 0)
        m    - SH order (-n <= m <= n)
        w    - workspace

Return: index in [0,nnm-1]
*/

static size_t
sqfilt_nmidx(const size_t type, const size_t n, const int m, const mag_sqfilt_workspace *w)
{
  int ns;
  size_t mmax, *baseptr;
  size_t base; /* index of block for this n */
  size_t type_offset; /* offset for coefficient type (internal/external) */
  int offset;  /* offset within block for this m */
  size_t nmidx;

  if (type == 0)
    {
      mmax = w->mmax_int;
      baseptr = w->base_int;
      type_offset = w->int_offset;
    }
  else if (type == 1)
    {
      mmax = w->mmax_ext;
      baseptr = w->base_ext;
      type_offset = w->ext_offset;
    }
  else
    {
      fprintf(stderr, "sqfilt_nmidx: error: unknown coefficient type\n");
      return 0;
    }

  ns = (int) GSL_MIN(n, mmax);

  if (n == 0)
    {
      fprintf(stderr, "sqfilt_nmidx: error: n = 0\n");
      return 0;
    }
  else if (abs(m) > (int) mmax)
    {
      fprintf(stderr, "sqfilt_nmidx: error: m = %d\n", m);
      return 0;
    }

  base = baseptr[n]; /* precomputed */
  offset = m + ns;

  nmidx = base + offset + type_offset;

  return nmidx;
} /* sqfilt_nmidx() */
