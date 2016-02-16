/*
 * stage1.c
 *
 * 1. Read tiegcm data file
 * 2. For each time step t_k, invert B(t_k) for SH coefficients
 *
 * ./stage1 <-i tiegcm_nc_file>
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <assert.h>
#include <errno.h>

#include <lapacke/lapacke.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "apex.h"
#include "common.h"
#include "magdata.h"
#include "oct.h"
#include "poltor.h"
#include "tiegcm.h"

#define USE_SYNTH_DATA              0

int
print_coefficients(poltor_workspace *w)
{
#if USE_SYNTH_DATA

  /* print and check synthetic coefficients */
  poltor_synth_print(w);

#else

  const double b = w->b;
  const double R = w->R;
  const size_t nmax_int = GSL_MIN(3, w->nmax_int);
  const size_t mmax_int = GSL_MIN(3, w->mmax_int);
  const size_t nmax_ext = GSL_MIN(2, w->nmax_ext);
  const size_t mmax_ext = GSL_MIN(2, w->mmax_ext);
  const size_t nmax_sh = GSL_MIN(2, w->nmax_sh);
  const size_t mmax_sh = GSL_MIN(2, w->mmax_sh);
  const size_t nmax_tor = GSL_MIN(3, w->nmax_tor);
  const size_t mmax_tor = GSL_MIN(3, w->mmax_tor);
  size_t n;

  /* print internal poloidal coefficients */
  fprintf(stderr, "Internal poloidal coefficients:\n");
  for (n = 1; n <= nmax_int; ++n)
    {
      int ni = (int) GSL_MIN(n, mmax_int);
      int m;

      /*
       * only need to print positive m coefficients since
       * c_{n,-m} = c_{nm}
       */
      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          gsl_complex coef = poltor_get(cidx, w);
          double qnm = GSL_REAL(coef);
          /* convert from poloidal to potential representation */
          double gnm = -(double)n / (2.0*n + 1.0) * qnm *
                       pow(b / R, n + 3.0);

          fprintf(stderr, "g(%2zu,%2d) = %12g q(%2zu,%2d) = %12g [nT]\n",
                  n, m, gnm,
                  n, m, qnm);
        }
    }

  /* print external poloidal coefficients */
  fprintf(stderr, "External poloidal coefficients:\n");
  for (n = 1; n <= nmax_ext; ++n)
    {
      int ni = (int) GSL_MIN(n, mmax_ext);
      int m;

      /*
       * only need to print positive m coefficients since
       * c_{n,-m} = c_{nm}
       */
      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
          gsl_complex coef = poltor_get(cidx, w);
          double knm = GSL_REAL(coef);

          fprintf(stderr, "k(%2zu,%2d) = %12g [nT]\n",
                  n,
                  m,
                  knm);
        }
    }

#endif /* USE_SYNTH_DATA */

  return 0;
} /* print_coefficients() */

int
print_residuals(const char *filename, poltor_workspace *w)
{
  size_t i, j;
  FILE *fp;
  magdata *data = w->data;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_residuals: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: time (decimal years)\n", i++);
  fprintf(fp, "# Field %zu: local time (hours)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: longitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: geocentric latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: X observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: X residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal X model (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal Y model (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal Z model (nT)\n", i++);
  fprintf(fp, "# Field %zu: external poloidal X model (nT)\n", i++);
  fprintf(fp, "# Field %zu: external poloidal Y model (nT)\n", i++);
  fprintf(fp, "# Field %zu: external poloidal Z model (nT)\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      double B_obs[3], B_model[3], B_int[4], B_sh[4], B_ext[4], B_tor[4], B_res[3];
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double lt = get_localtime(unix_time, data->phi[i]);

      /* store observation vector */
      magdata_residual(i, B_obs, data);

      /* compute individual magnetic field models at this point */
      poltor_eval_B_all(data->r[i], data->theta[i], data->phi[i], B_int, B_ext, B_sh, B_tor, w);

      /* compute total model and residual vector */
      for (j = 0; j < 3; ++j)
        {
          B_model[j] = B_int[j] + B_ext[j] + B_sh[j] + B_tor[j];
          B_res[j] = B_obs[j] - B_model[j];
        }

      fprintf(fp, "%f %.2f %.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
              satdata_epoch2year(data->t[i]),
              lt,
              data->r[i] - w->R,
              data->phi[i] * 180.0 / M_PI,
              90.0 - data->theta[i] * 180.0 / M_PI,
              data->qdlat[i],
              B_obs[0],
              B_obs[1],
              B_obs[2],
              B_res[0],
              B_res[1],
              B_res[2],
              B_int[0],
              B_int[1],
              B_int[2],
              B_ext[0],
              B_ext[1],
              B_ext[2]);
    }

  fclose(fp);

  return 0;
} /* print_residuals() */

/* plot current stream function grid */
int
print_chi(const char *filename, poltor_workspace *w)
{
  const double b = w->R + 110.0;
  const double lat_min = -80.0;
  const double lat_max = 80.0;
  const double lon_min = -180.0;
  const double lon_max = 180.0;
  const double phi_min = lon_min * M_PI / 180.0;
  const double phi_max = lon_max * M_PI / 180.0;
  const double theta_min = M_PI / 2.0 - lat_max * M_PI / 180.0;
  const double theta_max = M_PI / 2.0 - lat_min * M_PI / 180.0;
  double theta, phi;
  FILE *fp;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_chi: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: chi (external) (kA)\n", i++);

  fprintf(stderr, "print_chi: writing current stream function to %s...", filename);

  for (phi = phi_min; phi < phi_max; phi += 2.0 * M_PI / 180.0)
    {
      for (theta = theta_min; theta < theta_max; theta += 2.0 * M_PI / 180.0)
        {
          double chi; /* current stream function */

          poltor_eval_chi_ext(b, theta, phi, &chi, w);

          fprintf(fp, "%f %f %f\n",
                  phi * 180.0 / M_PI,
                  90.0 - theta * 180.0 / M_PI,
                  chi);
        }
      fprintf(fp, "\n");
    }

  fprintf(stderr, "done\n");

  fclose(fp);

  return 0;
}

magdata *
tiegcm_magdata(const size_t tidx, tiegcm_data *data)
{
  const size_t grid_size = data->nlon * data->nlat;
  const double eps = 1.0e-6;
  magdata *mdata;
  magdata_datum datum;
  size_t ilat, ilon;
  apex_workspace *apex_p;

  mdata = magdata_alloc(grid_size, R_EARTH_KM);
  if (!mdata)
    return 0;

  apex_p = apex_alloc(2016);

  magdata_datum_init(&datum);

  datum.t = satdata_timet2epoch(data->t[tidx]);
  datum.r = R_EARTH_KM;
  datum.flags = MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z;

  fprintf(stderr, "tiegcm_magdata: building magdata structure for time index %zu...", tidx);

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      double phi = data->glon[ilon] * M_PI / 180.0;

      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          double theta = M_PI / 2.0 - data->glat[ilat] * M_PI / 180.0;
          size_t idx = TIEGCM_BIDX(tidx, ilat, ilon, data);
          double qdlat, alon, alat;

          if (theta < eps)
            theta = eps;
          if (theta > M_PI - eps)
            theta = M_PI - eps;

          apex_transform_geodetic(theta, phi, 0.0, &alon, &alat, &qdlat,
                                  NULL, NULL, NULL, apex_p);

          datum.theta = theta;
          datum.phi = phi;
          datum.qdlat = qdlat;
          datum.B_nec[0] = data->Bx[idx] * 1.0e9;
          datum.B_nec[1] = data->By[idx] * 1.0e9;
          datum.B_nec[2] = data->Bz[idx] * 1.0e9;

          magdata_add(&datum, mdata);
        }
    }

  fprintf(stderr, "done\n");

  apex_free(apex_p);

  return mdata;
}

poltor_workspace *
main_invert(magdata *mdata)
{
  size_t nmax_int = 0;
  size_t mmax_int = 0;
  size_t nmax_ext = 60;
  size_t mmax_ext = 60;
  size_t nmax_sh = 0;
  size_t mmax_sh = 0;
  size_t nmax_tor = 0;
  size_t mmax_tor = 0;
  double alpha_int = 0.0;
  double alpha_sh = 0.0;
  double alpha_tor = 0.0;
  size_t robust_maxit = 1;
  const double R = R_EARTH_KM;
  const double b = R + 110.0;   /* radius of internal current shell (Sq+EEJ) */
  const double d = R + 350.0;   /* radius of current shell for gravity/diamag */
  double universal_time = 11.0; /* UT in hours for data selection */
  char *datamap_file = "datamap.dat";
  char *data_file = "data.dat";
  char *spectrum_file = "tiegcm.s";
  char *corr_file = "corr.dat";
  char *residual_file = "res.dat";
  char *output_file = NULL;
  char *chisq_file = NULL;
  char *lls_file = NULL;
  char *Lcurve_file = NULL;
  poltor_workspace *poltor_p;
  poltor_parameters params;
  struct timeval tv0, tv1;
  int print_data = 0;

  mmax_int = GSL_MIN(mmax_int, nmax_int);
  mmax_ext = GSL_MIN(mmax_ext, nmax_ext);
  mmax_sh = GSL_MIN(mmax_sh, nmax_sh);
  mmax_tor = GSL_MIN(mmax_tor, nmax_tor);

  fprintf(stderr, "main_invert: initializing spatial weighting histogram...");
  gettimeofday(&tv0, NULL);
  magdata_init(mdata);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* re-compute weights, nvec, nres based on flags update */
  fprintf(stderr, "main_invert: computing spatial weighting of data...");
  gettimeofday(&tv0, NULL);
  magdata_calc(mdata);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main_invert: setting unit spatial weights...");
  magdata_unit_weights(mdata);
  fprintf(stderr, "done\n");

  params.R = R;
  params.b = b;
  params.d = d;
  params.rmin = GSL_MAX(mdata->rmin, mdata->R + 250.0);
  params.rmax = GSL_MIN(mdata->rmax, mdata->R + 450.0);
  params.nmax_int = nmax_int;
  params.mmax_int = mmax_int;
  params.nmax_ext = nmax_ext;
  params.mmax_ext = mmax_ext;
  params.nmax_sh = nmax_sh;
  params.mmax_sh = mmax_sh;
  params.nmax_tor = nmax_tor;
  params.mmax_tor = mmax_tor;
  params.shell_J = 0;
  params.data = mdata;
  params.alpha_int = alpha_int;
  params.alpha_sh = alpha_sh;
  params.alpha_tor = alpha_tor;
  params.flags = POLTOR_FLG_QD_HARMONICS*0;

  poltor_p = poltor_alloc(&params);

#if USE_SYNTH_DATA
  fprintf(stderr, "main_invert: replacing with synthetic data...");
  gettimeofday(&tv0, NULL);
  poltor_synth(poltor_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
#endif

  {
    size_t maxiter = robust_maxit;
    size_t iter = 0;
    char buf[2048];

    while (iter++ < maxiter)
      {
        fprintf(stderr, "main_invert: ROBUST ITERATION %zu/%zu\n", iter, maxiter);

        /* build LS system */
        poltor_calc(poltor_p);

        /* solve LS system */
        poltor_solve(poltor_p);

        sprintf(buf, "%s.iter%zu", spectrum_file, iter);
        fprintf(stderr, "main_invert: printing spectrum to %s...", buf);
        poltor_print_spectrum(buf, poltor_p);
        fprintf(stderr, "done\n");
      }
  }

  print_coefficients(poltor_p);

#if 0
  print_chi("chi.dat", poltor_p);

  fprintf(stderr, "main_invert: printing correlation data to %s...", corr_file);
  print_correlation(corr_file, poltor_p);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main_invert: printing spectrum to %s...", spectrum_file);
  poltor_print_spectrum(spectrum_file, poltor_p);
  fprintf(stderr, "done\n");

  if (Lcurve_file)
    {
      fprintf(stderr, "main_invert: writing L-curve data to %s...", Lcurve_file);
      print_Lcurve(Lcurve_file, poltor_p);
      fprintf(stderr, "done\n");
    }

  if (output_file)
    {
      fprintf(stderr, "main_invert: writing output coefficients to %s...", output_file);
      poltor_write(output_file, poltor_p);
      fprintf(stderr, "done\n");
    }

  if (chisq_file)
    {
      fprintf(stderr, "main_invert: printing chisq/dof to %s...", chisq_file);
      print_chisq(chisq_file, poltor_p);
      fprintf(stderr, "done\n");
    }

  if (residual_file)
    {
      fprintf(stderr, "main_invert: printing residuals to %s...", residual_file);
      print_residuals(residual_file, poltor_p);
      fprintf(stderr, "done\n");
    }
#endif

  return poltor_p;
}

int
main_build_matrix(const tiegcm_data *data, msynth_workspace *msynth_p,
                  gsl_matrix *A)
{
  const size_t n = A->size1;
  const size_t p = A->size2;
  const double r = R_EARTH_KM;
  const double eps = 1.0e-6;
  size_t rowidx = 0;
  size_t ilon, ilat;

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      double phi = data->glon[ilon] * M_PI / 180.0;

      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          double theta = M_PI / 2.0 - data->glat[ilat] * M_PI / 180.0;
          gsl_vector_view v, row;

          if (theta < eps)
            theta = eps;
          if (theta > M_PI - eps)
            theta = M_PI - eps;

          /* compute external Green's functions */
          msynth_green_ext(r, theta, phi, msynth_p);

          /* copy Green's functions into matrix */

          v = gsl_vector_view_array(msynth_p->dX_ext, p);
          row = gsl_matrix_row(A, rowidx++);
          gsl_vector_memcpy(&row.vector, &v.vector);

          v = gsl_vector_view_array(msynth_p->dY_ext, p);
          row = gsl_matrix_row(A, rowidx++);
          gsl_vector_memcpy(&row.vector, &v.vector);

          v = gsl_vector_view_array(msynth_p->dZ_ext, p);
          row = gsl_matrix_row(A, rowidx++);
          gsl_vector_memcpy(&row.vector, &v.vector);
        }
    }

  assert(rowidx == n);

  return 0;
}

/*
main_build_rhs()
  Construct RHS vector for a given time index
*/

int
main_build_rhs(const size_t tidx, const tiegcm_data *data,
               gsl_vector *b)
{
  const size_t n = b->size;
  size_t ilon, ilat;
  size_t rowidx = 0;

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          size_t idx = TIEGCM_BIDX(tidx, ilat, ilon, data);

          gsl_vector_set(b, rowidx++, data->Bx[idx] * 1.0e9);
          gsl_vector_set(b, rowidx++, data->By[idx] * 1.0e9);
          gsl_vector_set(b, rowidx++, data->Bz[idx] * 1.0e9);
        }
    }

  assert(rowidx == n);

  return 0;
}

int
lapack_lls(const gsl_matrix * A, const gsl_matrix * B, gsl_matrix * X)
{
  int s;
  lapack_int m = A->size1;
  lapack_int n = A->size2;
  lapack_int nrhs = B->size2;
  lapack_int lda = A->size1;
  lapack_int ldb = B->size1;
  lapack_int rank;
  lapack_int *jpvt = malloc(n * sizeof(lapack_int));
  gsl_matrix *work_A = gsl_matrix_alloc(A->size2, A->size1);
  gsl_matrix *work_B = gsl_matrix_alloc(B->size2, B->size1);
  double rcond = 1.0e-6;

  gsl_matrix_transpose_memcpy(work_A, A);
  gsl_matrix_transpose_memcpy(work_B, B);

  s = LAPACKE_dgelsy(LAPACK_COL_MAJOR,
                     m,
                     n,
                     nrhs,
                     work_A->data,
                     lda,
                     work_B->data,
                     ldb,
                     jpvt,
                     rcond,
                     &rank);

  /* store solution in X */
  {
    gsl_matrix_view m = gsl_matrix_submatrix(work_B, 0, 0, X->size2, X->size1);
    gsl_matrix_transpose_memcpy(X, &m.matrix);
  }

  gsl_matrix_free(work_A);
  gsl_matrix_free(work_B);
  free(jpvt);

  return s;
}

int
main_proc(const char *filename, tiegcm_data *data)
{
  int status;
  const size_t nmax = 60;
  msynth_workspace *msynth_p = msynth_alloc2(nmax, nmax, 1, NULL);
  const size_t n = 3 * data->nlon * data->nlat; /* number of residuals */
  const size_t p = msynth_p->nnm;               /* number of external coefficients */
  const size_t nrhs = data->nt;                 /* number of right hand sides */
  gsl_matrix *A = gsl_matrix_alloc(n, p);       /* least squares matrix */
  gsl_matrix *B = gsl_matrix_alloc(n, nrhs);    /* right hand sides */
  gsl_matrix *X = gsl_matrix_alloc(p, nrhs);    /* solution vectors */
  gsl_vector *r = gsl_vector_alloc(n);          /* residual vector */
  size_t k;
  FILE *fp;
  struct timeval tv0, tv1;

  fprintf(stderr, "main_proc: %zu observations\n", n);
  fprintf(stderr, "main_proc: %zu model coefficients\n", p);

  /* construct least squares matrix (common for all timestamps) */
  fprintf(stderr, "main_proc: building least squares matrix A...");
  gettimeofday(&tv0, NULL);
  status = main_build_matrix(data, msynth_p, A);
  if (status)
    return status;
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fp = fopen(filename, "w");

  k = 1;
  fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", k++);
  fprintf(fp, "# Field %zu: k(1,0) (nT)\n", k++);
  fprintf(fp, "# Field %zu: k(1,1) (nT)\n", k++);
  fprintf(fp, "# Field %zu: k(2,0) (nT)\n", k++);
  fprintf(fp, "# Field %zu: k(2,1) (nT)\n", k++);
  fprintf(fp, "# Field %zu: k(2,2) (nT)\n", k++);

#if 1

  fprintf(stderr, "main_proc: building rhs vectors...");
  gettimeofday(&tv0, NULL);

  /* construct right hand side vectors */
  for (k = 0; k < data->nt; ++k)
    {
      gsl_vector_view b = gsl_matrix_column(B, k);

      /* construct rhs vector for time t_k */
      main_build_rhs(k, data, &b.vector);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* solve least squares system for all rhs vectors */
  fprintf(stderr, "main_proc: solving LS system with QR decomposition of A...");
  gettimeofday(&tv0, NULL);
  status = lapack_lls(A, B, X);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, s = %d)\n", time_diff(tv0, tv1), status);

  for (k = 0; k < data->nt; ++k)
    {
      size_t N;
#if 0
      gsl_vector_view b = gsl_matrix_column(B, k);
      gsl_vector_view x = gsl_matrix_column(X, k);

      /* compute r = b - A x */
      gsl_vector_memcpy(r, &b.vector);
      gsl_blas_dgemv(CblasNoTrans, -1.0, A, &x.vector, 1.0, r);

      fprintf(stderr, "main_proc: residual for timestamp (%zu/%zu): %.12e [nT]\n",
              k + 1, data->nt, gsl_blas_dnrm2(r));
#endif

      fprintf(fp, "%ld ", data->t[k]);

      for (N = 1; N <= 2; ++N)
        {
          int M = (int) N;
          int m;

          for (m = 0; m <= M; ++m)
            {
              size_t cidx = msynth_nmidx(N, m, msynth_p);
              double knm = gsl_matrix_get(X, cidx, k);

              fprintf(fp, "%f ", knm);
            }
        }

      putc('\n', fp);
      fflush(fp);
    }

#else

  for (k = 0; k < data->nt; ++k)
    {
      magdata *mdata;

      fprintf(stderr, "main_proc: processing grid for timestamp (%zu/%zu): %s",
              k + 1, data->nt, ctime(&data->t[k]));
      
      mdata = tiegcm_magdata(k, data);

      if (mdata)
        {
          poltor_workspace *poltor_p = main_invert(mdata);
          size_t n;

          fprintf(fp, "%ld ", data->t[k]);

          for (n = 1; n <= 2; ++n)
            {
              int M = (int) n;
              int m;

              for (m = 0; m <= M; ++m)
                {
                  size_t cidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, poltor_p);
                  gsl_complex coef = poltor_get(cidx, poltor_p);
                  double knm = GSL_REAL(coef);

                  fprintf(fp, "%f ", knm);
                }
            }

          putc('\n', fp);
          fflush(fp);

          poltor_free(poltor_p);
        }

      magdata_free(mdata);
    }

#endif

  msynth_free(msynth_p);
  gsl_matrix_free(A);
  gsl_matrix_free(B);
  gsl_matrix_free(X);

  fclose(fp);

  return 0;
}

int
main(int argc, char *argv[])
{
  tiegcm_data *data;
  struct timeval tv0, tv1;
  char *infile = NULL;
  char *outfile = "knm.dat";

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i tiegcm_nc_file>\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = tiegcm_read(infile, NULL);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->nt,
          time_diff(tv0, tv1));

  main_proc(outfile, data);

  tiegcm_free(data);

  return 0;
}
