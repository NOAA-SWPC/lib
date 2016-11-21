/*
 * secs1d.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <satdata/satdata.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

#include "common.h"
#include "oct.h"
#include "secs1d.h"
#include "track.h"

static int build_matrix_row_df(const double r, const double theta,
                               gsl_vector *X, gsl_vector *Z, secs1d_workspace *w);
static int build_matrix_row_cf(const double r, const double theta,
                               gsl_vector *Y, secs1d_workspace *w);
static int build_matrix_row_df_J(const double theta, gsl_vector *Y, secs1d_workspace *w);

/*
secs1d_alloc()
  Allocate secs 1d workspace

Inputs: flags        - SECS1D_FLG_xxx
        lmax         - maximum degree for Legendre functions in expansion
        R_iono       - radius of ionosphere (km)
        pole_spacing - along-orbit latitude spacing of SECS poles (degrees)

Return: pointer to workspace
*/

secs1d_workspace *
secs1d_alloc(const size_t flags, const size_t lmax, const double R_iono, const double pole_spacing)
{
  secs1d_workspace *w;
  const double max_lat = 60.0;
  const size_t npoles = (size_t) ((2.0 * max_lat) / pole_spacing + 1.0);
  const double dtheta = pole_spacing * M_PI / 180.0;
  const double theta_min = M_PI / 2.0 - max_lat * M_PI / 180.0;
  const double theta_max = M_PI / 2.0 + max_lat * M_PI / 180.0;
  size_t i;

  w = calloc(1, sizeof(secs1d_workspace));
  if (!w)
    return 0;

  w->nmax = 30000;
  w->n = 0;
  w->p = 0;
  w->npoles = npoles;
  w->R_iono = R_iono;
  w->lmax = lmax;
  w->flags = flags;
  w->df_offset = 0;
  w->cf_offset = 0;
  w->dtheta = dtheta;

  if (flags & SECS1D_FLG_FIT_DF)
    {
      w->df_offset = w->p;
      w->p += npoles;
    }

  if (flags & SECS1D_FLG_FIT_CF)
    {
      w->cf_offset = w->p;
      w->p += npoles;
    }

  w->X = gsl_matrix_alloc(w->nmax, w->p);
  w->c = gsl_vector_alloc(w->p);
  w->rhs = gsl_vector_alloc(w->nmax);
  w->wts = gsl_vector_alloc(w->nmax);
  w->cov = gsl_matrix_alloc(w->p, w->p);
  w->multifit_p = gsl_multifit_linear_alloc(w->nmax, w->p);

  /* regularization matrix L */
  {
    const size_t k = 1;
    w->L = gsl_matrix_alloc(w->p - k, w->p);
    w->Ltau = gsl_vector_alloc(w->p - k);
    w->M = gsl_matrix_alloc(w->nmax, w->p);
    w->Xs = gsl_matrix_alloc(w->nmax, w->p);
    w->bs = gsl_vector_alloc(w->nmax);
    w->cs = gsl_vector_alloc(w->p);

    gsl_multifit_linear_Lk(w->p, k, w->L);
    print_octave(w->L, "L");
    gsl_multifit_linear_L_decomp(w->L, w->Ltau);
  }

  w->theta0 = malloc(npoles * sizeof(double));
  w->Pltheta0 = malloc((w->lmax + 1) * sizeof(double));
  w->Pltheta = malloc((w->lmax + 1) * sizeof(double));

  /* fill in theta0 array */
  for (i = 0; i < npoles; ++i)
    w->theta0[i] = theta_min + i * dtheta;

  return w;
}

void
secs1d_free(secs1d_workspace *w)
{
  if (w->X)
    gsl_matrix_free(w->X);

  if (w->c)
    gsl_vector_free(w->c);

  if (w->rhs)
    gsl_vector_free(w->rhs);

  if (w->wts)
    gsl_vector_free(w->wts);

  if (w->cov)
    gsl_matrix_free(w->cov);

  if (w->L)
    gsl_matrix_free(w->L);

  if (w->Ltau)
    gsl_vector_free(w->Ltau);

  if (w->M)
    gsl_matrix_free(w->M);

  if (w->Xs)
    gsl_matrix_free(w->Xs);

  if (w->bs)
    gsl_vector_free(w->bs);

  if (w->cs)
    gsl_vector_free(w->cs);

  if (w->theta0)
    free(w->theta0);

  if (w->Pltheta0)
    free(w->Pltheta0);

  if (w->Pltheta)
    free(w->Pltheta);

  if (w->multifit_p)
    gsl_multifit_linear_free(w->multifit_p);

  free(w);
}

int
secs1d_add_track(const track_data *tptr, const satdata_mag *data,
                 secs1d_workspace *w)
{
  size_t rowidx = w->n;
  size_t i;

  if (w->flags & SECS1D_FLG_FIT_DF)
    {
      for (i = 0; i < tptr->n; ++i)
        {
          size_t didx = i + tptr->start_idx;
          double r = data->altitude[didx] + data->R;
          double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
          double wi = 1.0;
          gsl_vector_view vx = gsl_matrix_row(w->X, rowidx);
          gsl_vector_view vz = gsl_matrix_row(w->X, rowidx + 1);

          if (!SATDATA_AvailableData(data->flags[didx]))
            continue;

          /* fit only low-latitude data */
          if (fabs(data->qdlat[didx]) > SECS1D_QDMAX)
            continue;

          /* upweight equatorial data */
          if (fabs(data->qdlat[didx]) < 10.0)
            wi *= 1.0;

          /* set rhs vector */
          gsl_vector_set(w->rhs, rowidx, tptr->Bx[i]);
          gsl_vector_set(w->rhs, rowidx + 1, tptr->Bz[i]);

          /* set weight vector */
          gsl_vector_set(w->wts, rowidx, wi);
          gsl_vector_set(w->wts, rowidx + 1, wi);

          /* build 2 rows of the LS matrix for DF SECS */
          build_matrix_row_df(r, theta, &vx.vector, &vz.vector, w);

          rowidx += 2;
        }
    }

  if (w->flags & SECS1D_FLG_FIT_CF)
    {
      for (i = 0; i < tptr->n; ++i)
        {
          size_t didx = i + tptr->start_idx;
          double r = data->altitude[didx] + data->R;
          double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
          double wi = 1.0;
          gsl_vector_view vy = gsl_matrix_row(w->X, rowidx);

          if (!SATDATA_AvailableData(data->flags[didx]))
            continue;

          /* fit only low-latitude data */
          if (fabs(data->qdlat[didx]) > SECS1D_QDMAX)
            continue;

          /* set rhs vector */
          gsl_vector_set(w->rhs, rowidx, tptr->By[i]);

          /* set weight vector */
          gsl_vector_set(w->wts, rowidx, wi);

          /* build 1 row of the LS matrix for CF SECS */
          build_matrix_row_cf(r, theta, &vy.vector, w);

          rowidx += 1;
        }
    }

  w->n = rowidx;

  return 0;
}

/*
secs1d_fit()
  Fit 1D SECS to previously added tracks

Inputs: w - workspace

Return: success/error

Notes:
1) Data must be added to workspace via
secs1d_add_track()
*/

int
secs1d_fit(secs1d_workspace *w)
{
  const size_t npts = 200;
  const double tol = 1.0e-4;
  gsl_vector *reg_param = gsl_vector_alloc(npts);
  gsl_vector *rho = gsl_vector_alloc(npts);
  gsl_vector *eta = gsl_vector_alloc(npts);
  gsl_vector *G = gsl_vector_alloc(npts);
  gsl_matrix_view A = gsl_matrix_submatrix(w->X, 0, 0, w->n, w->p);
  gsl_vector_view b = gsl_vector_subvector(w->rhs, 0, w->n);
  gsl_vector_view wts = gsl_vector_subvector(w->wts, 0, w->n);
  double lambda_gcv, lambda_l, G_gcv;
  double rnorm, snorm;
  size_t i;
  const size_t m = w->L->size1;
  gsl_matrix_view M = gsl_matrix_submatrix(w->M, 0, 0, w->n, w->p);
  gsl_matrix_view As = gsl_matrix_submatrix(w->Xs, 0, 0, w->n - w->p + m, m);
  gsl_vector_view bs = gsl_vector_subvector(w->bs, 0, w->n - w->p + m);
  gsl_vector_view cs = gsl_vector_subvector(w->cs, 0, m);

  print_octave(&A.matrix, "A");
  printv_octave(&b.vector, "b");

#if 0 /* TSVD */

  {
    double chisq;
    size_t rank;

    gsl_multifit_wlinear_tsvd(&A.matrix, &wts.vector, &b.vector, tol, w->c, w->cov,
                              &chisq, &rank, w->multifit_p);

    rnorm = sqrt(chisq);
    snorm = gsl_blas_dnrm2(w->c);

    fprintf(stderr, "secs1d_fit: rank = %zu/%zu\n", rank, w->p);
  }

#else /* Tikhonov / L-curve */

#if 0
  /* convert to standard form */
  gsl_multifit_linear_applyW(&A.matrix, &wts.vector, &b.vector, &A.matrix, &b.vector);
#else

  gsl_multifit_linear_wstdform2(w->L, w->Ltau, &A.matrix, &wts.vector, &b.vector,
                                &As.matrix, &bs.vector, &M.matrix, w->multifit_p);

#endif

  /* compute SVD of A */
  gsl_multifit_linear_svd(&As.matrix, w->multifit_p);

  /* compute L-curve */
  gsl_multifit_linear_lcurve(&bs.vector, reg_param, rho, eta, w->multifit_p);
  gsl_multifit_linear_lcorner(rho, eta, &i);
  lambda_l = gsl_vector_get(reg_param, i);

  /* compute GCV curve */
  gsl_multifit_linear_gcv(&bs.vector, reg_param, G, &lambda_gcv, &G_gcv, w->multifit_p);

  /* solve regularized system with lambda_l */
  lambda_l = tol * gsl_vector_get(w->multifit_p->S, 0);
  gsl_multifit_linear_solve(lambda_l, &As.matrix, &bs.vector, &cs.vector, &rnorm, &snorm, w->multifit_p);

  /* convert back to general form */
  gsl_multifit_linear_wgenform2(w->L, w->Ltau, &A.matrix, &wts.vector, &b.vector, &cs.vector, &M.matrix,
                                w->c, w->multifit_p);

  fprintf(stderr, "lambda_l = %.12e\n", lambda_l);
  fprintf(stderr, "lambda_gcv = %.12e\n", lambda_gcv);

  for (i = 0; i < npts; ++i)
    {
      fprintf(stdout, "%e %e %e %e\n",
              gsl_vector_get(reg_param, i),
              gsl_vector_get(rho, i),
              gsl_vector_get(eta, i),
              gsl_vector_get(G, i));
    }

#endif

  fprintf(stderr, "rnorm = %.12e\n", rnorm);
  fprintf(stderr, "snorm = %.12e\n", snorm);

  printv_octave(w->c, "c");

  gsl_vector_free(reg_param);
  gsl_vector_free(rho);
  gsl_vector_free(eta);
  gsl_vector_free(G);

  return 0;
}

/*
secs1d_eval_B()
  Evaluate magnetic field at a given (r,theta) using
previously computed 1D SECS coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        B     - (output) magnetic field vector (nT)
        w     - workspace

Notes:
1) w->c must contain 1D SECS coefficients
*/

int
secs1d_eval_B(const double r, const double theta,
              double B[3], secs1d_workspace *w)
{
  gsl_vector_view vx = gsl_matrix_row(w->X, 0);
  gsl_vector_view vy = gsl_matrix_row(w->X, 1);
  gsl_vector_view vz = gsl_matrix_row(w->X, 2);

  B[0] = 0.0;
  B[1] = 0.0;
  B[2] = 0.0;

  if (w->flags & SECS1D_FLG_FIT_DF)
    {
      build_matrix_row_df(r, theta, &vx.vector, &vz.vector, w);

      gsl_blas_ddot(&vx.vector, w->c, &B[0]);
      gsl_blas_ddot(&vz.vector, w->c, &B[2]);
    }

  if (w->flags & SECS1D_FLG_FIT_CF)
    {
      build_matrix_row_cf(r, theta, &vy.vector, w);

      gsl_blas_ddot(&vy.vector, w->c, &B[1]);
    }

  return 0;
}

/*
secs1d_eval_J()
  Evaluate current density at a given (r,theta) using
previously computed 1D SECS coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        J     - (output) current density vector [A/km]
        w     - workspace

Notes:
1) w->c must contain 1D SECS coefficients
*/

int
secs1d_eval_J(const double r, const double theta,
              double J[3], secs1d_workspace *w)
{
  gsl_vector_view vx = gsl_matrix_row(w->X, 0);
  gsl_vector_view vy = gsl_matrix_row(w->X, 1);
  gsl_vector_view vz = gsl_matrix_row(w->X, 2);
  size_t i;

  J[0] = 0.0;
  J[1] = 0.0;
  J[2] = 0.0;

  if (w->flags & SECS1D_FLG_FIT_DF)
    {
      build_matrix_row_df_J(theta, &vy.vector, w);

      gsl_blas_ddot(&vy.vector, w->c, &J[1]);
    }

#if 0
  if (w->flags & SECS1D_FLG_FIT_CF)
    {
      build_matrix_row_cf(r, theta, &vy.vector, w);

      gsl_blas_ddot(&vy.vector, w->c, &B[1]);
    }
#endif

  /* convert to units of A/km */
  for (i = 0; i < 3; ++i)
    J[i] *= 1.0e3;

  return 0;
}

/*
secs1d_print_track()
  Print a track with SECS model to a file
*/

int
secs1d_print_track(const int header, FILE *fp, const track_data *tptr,
                   const satdata_mag *data, secs1d_workspace *w)
{
  size_t i;

  if (header)
    {
      i = 1;
      fprintf(fp, "# Field %zu: timestamp\n", i++);
      fprintf(fp, "# Field %zu: radius (km)\n", i++);
      fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: X measurement (nT)\n", i++);
      fprintf(fp, "# Field %zu: Y measurement (nT)\n", i++);
      fprintf(fp, "# Field %zu: Z measurement (nT)\n", i++);
      fprintf(fp, "# Field %zu: SECS1D B_X (nT)\n", i++);
      fprintf(fp, "# Field %zu: SECS1D B_Y (nT)\n", i++);
      fprintf(fp, "# Field %zu: SECS1D B_Z (nT)\n", i++);
      fprintf(fp, "# Field %zu: SECS1D J_Y\n", i++);
      return 0;
    }

  for (i = 0; i < tptr->n; ++i)
    {
      size_t didx = i + tptr->start_idx;
      double r = data->altitude[didx] + data->R;
      double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
      double B_secs[3], J_secs[3];

      if (!SATDATA_AvailableData(data->flags[didx]))
        continue;

      /* fit only low-latitude data */
      if (fabs(data->qdlat[didx]) > SECS1D_QDMAX)
        continue;

      /* compute SECS model prediction */
      secs1d_eval_B(r, theta, B_secs, w);
      secs1d_eval_J(r, theta, J_secs, w);

      fprintf(fp, "%ld %f %f %f %f %f %f %f %f %f %f %f\n",
              satdata_epoch2timet(data->t[didx]),
              r,
              data->longitude[didx],
              data->latitude[didx],
              data->qdlat[didx],
              tptr->Bx[i],
              tptr->By[i],
              tptr->Bz[i],
              B_secs[0],
              B_secs[1],
              B_secs[2],
              J_secs[1]);
    }

  fprintf(fp, "\n\n");

  return 0;
}

int
secs1d_green_df_init(const double theta, secs1d_workspace *w)
{
  /* compute P_l(cos(theta)) */
  gsl_sf_legendre_Pl_array(w->lmax, cos(theta), w->Pltheta);

  return 0;
}

/*
secs1d_green_df()
  Compute magnetic field Green's function for a single divergence-free
1D SECS

Inputs: r      - radius (km)
        theta  - colatitude (radians)
        theta0 - pole position (radians)
        B      - (output) magnetic field Green's function (X,Y,Z)
        w      - workspace

Return: success/error

Notes:
1) The array w->Pltheta must be filled with P_l(cos(theta))
   prior to calling this function (see secs1d_green_df_init)

2) DF 1D SECS has no Y component, so B[1] = 0 always
*/

int
secs1d_green_df(const double r, const double theta, const double theta0,
                double B[3], secs1d_workspace *w)
{
  const double costheta = cos(theta);
  B[0] = 0.0;
  B[1] = 0.0;
  B[2] = 0.0;

  /* compute P_l(cos(theta0)) */
  gsl_sf_legendre_Pl_array(w->lmax, cos(theta0), w->Pltheta0);

  if (r > w->R_iono)
    {
      double ratio = w->R_iono / r;
      double rterm = ratio; /* (R / r)^{l+1} */
      size_t l;

      for (l = 1; l <= w->lmax; ++l)
        {
          double Pl1 = gsl_sf_legendre_Plm(l, 1, costheta);

          /* (R / r)^{l+1} */
          rterm *= ratio;

          B[0] += rterm / (l + 1.0) * w->Pltheta0[l] * Pl1;
          B[2] -= rterm * w->Pltheta0[l] * w->Pltheta[l];
        }

      B[0] *= SECS1D_MU_0 / (2.0 * r);
      B[2] *= SECS1D_MU_0 / (2.0 * r);
    }
  else
    {
    }

  return GSL_SUCCESS;
}

/*
secs1d_green_df_J()
  Compute current density Green's function for a single divergence-free
1D SECS

Inputs: theta  - colatitude (radians)
        theta0 - pole position (radians)
        K      - (output) current density Green's function (X,Y,Z)
        w      - workspace

Return: success/error
*/

int
secs1d_green_df_J(const double theta, const double theta0, double K[3], secs1d_workspace *w)
{
  const double dtheta_2 = 0.5 * w->dtheta;
  const double d = theta - theta0;

  K[0] = 0.0;
  K[2] = 0.0;

  if (fabs(d) > dtheta_2)
    {
      double cs = cos(theta);
      double sn = sin(theta);

      K[1] = (cs + GSL_SIGN(d)) / sn;
    }
  else
    {
      /* use linear interpolation near 1D SECS pole */
      double north = -tan(0.5 * (theta0 - dtheta_2));
      double south = 1.0 / tan(0.5 * (theta0 + dtheta_2));
      double f = (south - north) / w->dtheta;

      K[1] = 0.5 * (north + south) + d * f;
    }

  K[1] /= 2.0 * w->R_iono;

  return GSL_SUCCESS;
}

/*
secs1d_green_cf()
  Compute magnetic field Green's function for a single divergence-free
1D SECS

Inputs: r      - radius (km)
        theta  - colatitude (radians)
        theta0 - pole position (radians)
        B      - (output) magnetic field Green's function (X,Y,Z)
        w      - workspace

Return: success/error
*/

int
secs1d_green_cf(const double r, const double theta, const double theta0,
                double B[3], secs1d_workspace *w)
{
  B[0] = 0.0;
  B[1] = 0.0;
  B[2] = 0.0;

  if (r > w->R_iono)
    {
      if (theta > theta0)
        B[1] = -tan(0.5 * (M_PI - theta));
      else
        B[1] = tan(0.5 * theta);

      B[1] *= SECS1D_MU_0 / (2.0 * r);
    }

  return GSL_SUCCESS;
}

/*
secs1d_green_cf_J()
  Compute current density Green's function for a single curl-free
1D SECS

Inputs: r      - radius (km)
        theta  - colatitude (radians)
        theta0 - pole position (radians)
        K      - (output) current density Green's function (X,Y,Z)
        w      - workspace

Return: success/error
*/

int
secs1d_green_cf_J(const double r, const double theta, const double theta0, double K[3], secs1d_workspace *w)
{
  K[1] = 0.0;
  K[2] = 0.0;

  if (theta < theta0)
    K[0] = tan(0.5 * theta);
  else
    K[0] = -tan(0.5 * (M_PI - theta));

  K[0] /= (2.0 * w->R_iono);

  if (r >= w->R_iono)
    K[2] = -1.0 / (2.0 * r * r);

  return GSL_SUCCESS;
}

/*
build_matrix_row_df()
  Build matrix rows corresponding to DF SECS

Inputs: r     - radius (km)
        theta - colatitude (radians)
        X     - (output) X Green's functions for B_df
        Z     - (output) Z Green's functions for B_df
        w     - workspace

Notes:
1) B field for DF 1D SECS does not have Y component
*/

static int
build_matrix_row_df(const double r, const double theta,
                    gsl_vector *X, gsl_vector *Z, secs1d_workspace *w)
{
  size_t i;

  gsl_vector_set_zero(X);
  gsl_vector_set_zero(Z);

  /* initialize needed arrays */
  secs1d_green_df_init(theta, w);

  for (i = 0; i < w->npoles; ++i)
    {
      double theta0 = w->theta0[i];
      double B[3];

      /* compute divergence-free 1D SECS Green's functions */
      secs1d_green_df(r, theta, theta0, B, w);

      gsl_vector_set(X, w->df_offset + i, B[0]);
      gsl_vector_set(Z, w->df_offset + i, B[2]);
    }

  return 0;
}

/*
build_matrix_row_cf()
  Build matrix rows corresponding to CF SECS

Inputs: r - radius (km)
        theta - colatitude (radians)
        Y     - (output) Y Green's functions for B_cf
        w     - workspace

Notes:
1) B field for CF 1D SECS only has Y component
*/

static int
build_matrix_row_cf(const double r, const double theta,
                    gsl_vector *Y, secs1d_workspace *w)
{
  size_t i;

  gsl_vector_set_zero(Y);

  for (i = 0; i < w->npoles; ++i)
    {
      double theta0 = w->theta0[i];
      double B[3];

      /* compute curl-free 1D SECS Green's functions */
      secs1d_green_cf(r, theta, theta0, B, w);

      gsl_vector_set(Y, w->cf_offset + i, B[1]);
    }

  return 0;
}

/*
build_matrix_row_df_J()
  Build matrix rows corresponding to DF SECS (current)

Inputs: theta - colatitude (radians)
        Y     - (output) Y Green's functions for J_df
        w     - workspace

Notes:
1) J vector for DF 1D SECS does not have X, Z components
*/

static int
build_matrix_row_df_J(const double theta, gsl_vector *Y, secs1d_workspace *w)
{
  size_t i;

  gsl_vector_set_zero(Y);

  for (i = 0; i < w->npoles; ++i)
    {
      double theta0 = w->theta0[i];
      double J[3];

      /* compute divergence-free 1D SECS Green's functions */
      secs1d_green_df_J(theta, theta0, J, w);

      gsl_vector_set(Y, w->df_offset + i, J[1]);
    }

  return 0;
}
