/*
 * mag_sqfilt_vector.c
 *
 * Fit spherical harmonic model to satellite scalar or vector residuals
 * to remove Sq and external fields
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

static int sqfilt_linear_init(mag_workspace *mag_p, mag_sqfilt_vector_workspace *w);
static int sqfilt_linear_matrix_row(const double r, const double thetaq,
                                    const double phi, gsl_vector *vx,
                                    gsl_vector *vy, gsl_vector *vz, mag_workspace *w);
static int sqfilt_calc_B2(const gsl_vector *c, mag_workspace *w);
static size_t sqfilt_nmidx(const size_t type, const size_t n, const int m,
                           const mag_sqfilt_vector_workspace *w);

mag_sqfilt_vector_workspace *
mag_sqfilt_vector_alloc(const size_t nmax_int, const size_t mmax_int,
                        const size_t nmax_ext, const size_t mmax_ext)
{
  const size_t ndata = 30000; /* maximum data expected in 1 track */
  size_t l, n;
  mag_sqfilt_vector_workspace *w;

  w = calloc(1, sizeof(mag_sqfilt_vector_workspace));
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
  w->c = gsl_vector_alloc(w->p);
  w->rhs = gsl_vector_alloc(ndata);
  w->work1 = gsl_vector_alloc(w->p);
  w->work2 = gsl_vector_alloc(w->p);
  w->work3 = gsl_vector_alloc(w->p);
  w->L = gsl_vector_alloc(w->p);

  w->ntot = ndata;
  w->n = 0;

  w->nreg = 200;
  w->rho = gsl_vector_alloc(w->nreg);
  w->eta = gsl_vector_alloc(w->nreg);
  w->reg_param = gsl_vector_alloc(w->nreg);

  w->multifit_workspace_p = gsl_multifit_linear_alloc(ndata, w->p);

  {
    const double beta = -0.1;

    for (n = 1; n <= nmax_int; ++n)
      {
        int M = (int) GSL_MIN(n, mmax_int);
        int m;
        double val = pow((double) n, beta);

        for (m = -M; m <= M; ++m)
          {
            size_t idx = sqfilt_nmidx(0, n, m, w);
            gsl_vector_set(w->L, idx, val);
          }
      }

    for (n = 1; n <= nmax_ext; ++n)
      {
        int M = (int) GSL_MIN(n, mmax_ext);
        int m;
        double val = pow((double) n, beta);

        for (m = -M; m <= M; ++m)
          {
            size_t idx = sqfilt_nmidx(1, n, m, w);
            gsl_vector_set(w->L, idx, val);
          }
      }
  }

  return w;
} /* mag_sqfilt_vector_alloc() */

void
mag_sqfilt_vector_free(mag_sqfilt_vector_workspace *w)
{
  if (w->base_int)
    free(w->base_int);

  if (w->base_ext)
    free(w->base_ext);

  if (w->X)
    gsl_matrix_free(w->X);

  if (w->c)
    gsl_vector_free(w->c);

  if (w->rhs)
    gsl_vector_free(w->rhs);

  if (w->work1)
    gsl_vector_free(w->work1);

  if (w->work2)
    gsl_vector_free(w->work2);

  if (w->work3)
    gsl_vector_free(w->work3);

  if (w->multifit_workspace_p)
    gsl_multifit_linear_free(w->multifit_workspace_p);

  if (w->rho)
    gsl_vector_free(w->rho);

  if (w->eta)
    gsl_vector_free(w->eta);

  if (w->reg_param)
    gsl_vector_free(w->reg_param);

  if (w->L)
    gsl_vector_free(w->L);

  free(w);
}

/*
mag_sqfilt_vector()
  Fit a simple spherical harmonic model to a single satellite
track of scalar B^(1) residual data, excluding the EEJ at low-latitudes.

Notes:
1) On output, w->rnorm and w->snorm contain the
residual norm and solution norm

2) On output, mag_p->track.{X,Y,Z}2 contains the B^(2) residuals and
mag_p->track.Sq_model contains the (M + K) Sq model
*/

int
mag_sqfilt_vector(mag_workspace *mag_p, mag_sqfilt_vector_workspace *w)
{
  int s;
  gsl_vector_view bv;
  gsl_matrix_view Xv;
  double lambda; /* final damping parameter */
  double smax;

  s = sqfilt_linear_init(mag_p, w);
  if (s)
    return s;

  Xv = gsl_matrix_submatrix(w->X, 0, 0, w->n, w->p);
  bv = gsl_vector_subvector(w->rhs, 0, w->n);

  /* convert to standard form */
  s = gsl_multifit_linear_stdform1(w->L, &Xv.matrix, &bv.vector,
                                   &Xv.matrix, &bv.vector, w->multifit_workspace_p);
  if (s)
    return s;

  /* compute SVD of LS matrix */
  s = gsl_multifit_linear_svd(&Xv.matrix, w->multifit_workspace_p);
  if (s)
    return s;

  /* compute L-curve */
  s = gsl_multifit_linear_lcurve(&bv.vector, w->reg_param, w->rho, w->eta, w->multifit_workspace_p);
  if (s)
    return s;

  s = gsl_multifit_linear_lcorner(w->rho, w->eta, &(w->reg_idx));
  if (s)
    return s;

  smax = gsl_vector_get(w->multifit_workspace_p->S, 0);
  lambda = GSL_MAX(gsl_vector_get(w->reg_param, w->reg_idx), 0.05*smax);
  w->reg_idx = bsearch_desc_double(w->reg_param->data, lambda, 0, w->nreg - 1);
  lambda = gsl_vector_get(w->reg_param, w->reg_idx);

  /* solve system with optimal lambda */
  s = gsl_multifit_linear_solve(lambda, &Xv.matrix, &bv.vector, w->c,
                                &(w->rnorm), &(w->snorm), w->multifit_workspace_p);
  if (s)
    return s;

  /* convert back to general form */
  s = gsl_multifit_linear_genform1(w->L, w->c, w->c, w->multifit_workspace_p);

  fprintf(stderr, "mag_sqfilt_vector: final regularization factor = %g (smax = %g)\n", lambda, smax);
  fprintf(stderr, "mag_sqfilt_vector: multifit residual norm = %.6e\n", w->rnorm);
  fprintf(stderr, "mag_sqfilt_vector: multifit solution norm = %.6e\n", w->snorm);

#if 0
  /* calculate F^(2) residuals (eq 9 of paper) */
  sqfilt_calc_F2(w->c, mag_p);
#endif

  return s;
} /* mag_sqfilt_vector() */

/****************************************************
 *       INTERNAL ROUTINES                          *
 ****************************************************/

/*
sqfilt_linear_init()
  Fit a simple spherical harmonic model to a single satellite
track of scalar B^(1) residual data, excluding the EEJ at low-latitudes.
This function initializes the LS design matrix and rhs vector

Inputs: mag_p - mag workspace
        w     - sqfilt workspace

Notes:
1) mag_p->track must be initialized on input with a complete satellite track

2) On output, w->n contains the number of data points to be fitted
*/

static int
sqfilt_linear_init(mag_workspace *mag_p, mag_sqfilt_vector_workspace *w)
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
  double qd_min, qd_max;
  double X1_min, X1_max;
  double Y1_min, Y1_max;
  double Z1_min, Z1_max;
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
          X1_min = track->X1[i];
          Y1_min = track->Y1[i];
          Z1_min = track->Z1[i];
        }
      else if (fabs(qdprev) < exclude_lat && fabs(qdcur) >= exclude_lat)
        {
          assert(found_upper == 0);

          found_upper = 1;
          qd_max = qdcur;
          X1_max = track->X1[i];
          Y1_max = track->Y1[i];
          Z1_max = track->Z1[i];
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
      double rhs_X = track->X1[i];
      double rhs_Y = track->Y1[i];
      double rhs_Z = track->Z1[i];
      gsl_vector_view vx, vy, vz;

      /* exclude points in the equatorial electrojet region */
      if (fabs(track->qdlat[i]) < exclude_lat)
        {
          /* keep every 10th point in the equatorial region */
          if (i % 10 == 0)
            {
              /* interpolate linearly across EEJ region */
              rhs_X = interp1d(qd_min, qd_max, X1_min, X1_max, track->qdlat[i]);
              rhs_Y = interp1d(qd_min, qd_max, Y1_min, Y1_max, track->qdlat[i]);
              rhs_Z = interp1d(qd_min, qd_max, Z1_min, Z1_max, track->qdlat[i]);
            }
          else
            continue;
        }

      /* set RHS values */
      gsl_vector_set(w->rhs, ndata, rhs_X);
      gsl_vector_set(w->rhs, ndata + 1, rhs_Y);
      gsl_vector_set(w->rhs, ndata + 2, rhs_Z);

      /* add green basis functions to current row of matrix */
      vx = gsl_matrix_row(w->X, ndata);
      vy = gsl_matrix_row(w->X, ndata + 1);
      vz = gsl_matrix_row(w->X, ndata + 2);
      sqfilt_linear_matrix_row(r, thetaq, phi,
                               &vx.vector, &vy.vector, &vz.vector,
                               mag_p);

      ndata += 3;
    }

  w->n = ndata;

  return s;
} /* sqfilt_linear_init() */

/*
sqfilt_linear_matrix_row()
  Add three rows to the least squares matrix

Inputs: r      - radius (km)
        thetaq - QD colatitude (radians)
        phi    - geographic longitude (radians)
        vx     - row of LS matrix for X basis functions
        vy     - row of LS matrix for Y basis functions
        vz     - row of LS matrix for Z basis functions
        w      - workspace
*/

static int
sqfilt_linear_matrix_row(const double r, const double thetaq,
                         const double phi, gsl_vector *vx,
                         gsl_vector *vy, gsl_vector *vz, mag_workspace *w)
{
  int s = 0;
  const double R = 6371.2;
  green_workspace *green_p = w->green_workspace_p;
  mag_sqfilt_vector_workspace *sqfilt_p = w->sqfilt_vector_workspace_p;
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

          gsl_vector_set(vx, idx, green_p->dX[gidx]);
          gsl_vector_set(vy, idx, green_p->dY[gidx]);
          gsl_vector_set(vz, idx, green_p->dZ[gidx]);
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

          gsl_vector_set(vx, idx, green_p->dX_ext[gidx]);
          gsl_vector_set(vy, idx, green_p->dY_ext[gidx]);
          gsl_vector_set(vz, idx, green_p->dZ_ext[gidx]);
        }
    }

  return s;
} /* sqfilt_linear_matrix_row() */

/*
sqfilt_calc_B2()
  Calculate B^(2) residual (eq. 9 of paper) by subtracting
the model (M + K)

Inputs: c - model coefficients
        w - workspace

Notes:
1) On output, w->track.{X,Y,Z}2 is modified to contain B^(2) residuals
*/

static int
sqfilt_calc_B2(const gsl_vector *c, mag_workspace *w)
{
#if 0
  mag_sqfilt_vector_workspace *sqfilt_p = w->sqfilt_vector_workspace_p;
  size_t i;
  mag_track *track = &(w->track);
  gsl_vector *vx = sqfilt_p->work1;
  gsl_vector *vy = sqfilt_p->work2;
  gsl_vector *vz = sqfilt_p->work3;

  for (i = 0; i < track->n; ++i)
    {
      double r = track->r[i];
      double thetaq = track->thetaq[i]; /* use QD colatitude */
      double phi = track->phi[i];
      double b_int[3];
      double val; /* Sq model value (M + K) */

      /* compute unit vector in internal field direction */
      b_int[0] = track->Bx_int[i] / track->F_int[i];
      b_int[1] = track->By_int[i] / track->F_int[i];
      b_int[2] = track->Bz_int[i] / track->F_int[i];

      /* compute basis functions */
      sqfilt_linear_matrix_row(r, thetaq, phi, b_int, v, w);

      /* compute internal model = b . M */
      if (sqfilt_p->p_int > 0)
        {
          gsl_vector_const_view c_int = gsl_vector_const_subvector(c, sqfilt_p->int_offset, sqfilt_p->p_int);
          gsl_vector_const_view vx_int = gsl_vector_const_subvector(vx, sqfilt_p->int_offset, sqfilt_p->p_int);
          gsl_vector_const_view vy_int = gsl_vector_const_subvector(vy, sqfilt_p->int_offset, sqfilt_p->p_int);
          gsl_vector_const_view vz_int = gsl_vector_const_subvector(vz, sqfilt_p->int_offset, sqfilt_p->p_int);

          gsl_blas_ddot(&vx_int.vector, &c_int.vector, &val);
          track->Sq_int_X[i] = val;

          gsl_blas_ddot(&vy_int.vector, &c_int.vector, &val);
          track->Sq_int_Y[i] = val;

          gsl_blas_ddot(&vz_int.vector, &c_int.vector, &val);
          track->Sq_int_Z[i] = val;
        }
      else
          track->Sq_int[i] = 0.0;

      /* compute external model = b . K */
      if (sqfilt_p->p_ext > 0)
        {
          gsl_vector_const_view c_ext = gsl_vector_const_subvector(c, sqfilt_p->ext_offset, sqfilt_p->p_ext);
          gsl_vector_const_view v_ext = gsl_vector_const_subvector(v, sqfilt_p->ext_offset, sqfilt_p->p_ext);

          gsl_blas_ddot(&v_ext.vector, &c_ext.vector, &val);
          track->Sq_ext[i] = val;
        }
      else
          track->Sq_ext[i] = 0.0;

      track->F2[i] = track->F1[i] - (track->Sq_int[i] + track->Sq_ext[i]);
    }
#endif

  return GSL_SUCCESS;
} /* sqfilt_calc_B2() */

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
sqfilt_nmidx(const size_t type, const size_t n, const int m, const mag_sqfilt_vector_workspace *w)
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
