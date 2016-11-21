/*
 * test.c
 *
 * 1. Synthesize magnetic field data for a single unit 1D SECS
 *    (Eqs 6-7 of Liisa's paper)
 * 2. Invert magnetic field data for internal Gauss coefficients
 * 3. Construct sheet current at 110km from Gauss coefficients and
 *    compare with analytic current (Eq. 5)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>

#include "common.h"
#include "green.h"
#include "lapack_wrapper.h"
#include "secs1d.h"

/* test df 1D SECS */
int
test_df(const double theta0, secs1d_workspace *w)
{
  int s = 0;
  const double eps = 1.0e-4;
  const double r = R_EARTH_KM + 450.0;
  const double dtheta = 5.0 * M_PI / 180.0;
  const double dphi = 5.0 * M_PI / 180.0;
  double theta, phi;
  const double nmax = 80;
  green_workspace *green_p = green_alloc(nmax, nmax, R_EARTH_KM);
  const size_t n = 100000;
  const size_t p = green_nnm(green_p);
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_vector *rhs = gsl_vector_alloc(n);
  gsl_vector *c = gsl_vector_alloc(p);
  gsl_matrix *cov = gsl_matrix_alloc(p, p);
  size_t rowidx = 0;

  for (theta = eps; theta <= M_PI - eps; theta += dtheta)
    {
      secs1d_green_df_init(theta, w);

      for (phi = 0.0; phi <= 2.0 * M_PI; phi += dphi)
        {
          double B[3];
          gsl_vector_view vx = gsl_matrix_row(X, rowidx);
          gsl_vector_view vy = gsl_matrix_row(X, rowidx + 1);
          gsl_vector_view vz = gsl_matrix_row(X, rowidx + 2);

          /* compute SECS Green's functions */
          secs1d_green_df(r, theta, theta0, B, w);

          /* compute internal potential field Green's functions */
          green_calc_int(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data, green_p);

          gsl_vector_set(rhs, rowidx, B[0]);
          gsl_vector_set(rhs, rowidx + 1, B[1]);
          gsl_vector_set(rhs, rowidx + 2, B[2]);

          rowidx += 3;
        }
    }

  {
    gsl_matrix_view m = gsl_matrix_submatrix(X, 0, 0, rowidx, p);
    gsl_vector_view b = gsl_vector_subvector(rhs, 0, rowidx);
    gsl_vector *res = gsl_vector_alloc(rowidx);
    int rank = 0;
    struct timeval tv0, tv1;
    double rms = 0.0;
    size_t nrms = 0;

    fprintf(stderr, "test_df: size of LS system: (%zu,%zu)\n", rowidx, p);

    fprintf(stderr, "test_df: solving LS system...");

    gettimeofday(&tv0, NULL);
    lapack_lls2(&m.matrix, &b.vector, c, &rank);
    gettimeofday(&tv1, NULL);
    gsl_multifit_linear_residuals(&m.matrix, &b.vector, c, res);

    fprintf(stderr, "done (%g seconds, rank = %d, residual = %.12e)\n",
            time_diff(tv0, tv1), rank, gsl_blas_dnrm2(res));

    phi = 0.0;
    for (theta = eps; theta <= M_PI - eps; theta += dtheta)
      {
        double B[3], B_model[3], K[3], K_model[3];
        gsl_vector_view vx = gsl_matrix_row(X, 0);
        gsl_vector_view vy = gsl_matrix_row(X, 1);
        gsl_vector_view vz = gsl_matrix_row(X, 2);

        /* compute SECS Green's functions */
        secs1d_green_df(r, theta, theta0, B, w);

        /* compute internal potential field Green's functions */
        green_calc_int(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data, green_p);

        gsl_blas_ddot(&vx.vector, c, &B_model[0]);
        gsl_blas_ddot(&vy.vector, c, &B_model[1]);
        gsl_blas_ddot(&vz.vector, c, &B_model[2]);

        /* compute current densities */
        secs1d_green_df_J(theta, theta0, K, w);
        green_eval_sheet_int(R_EARTH_KM + 110.0, theta, phi, c->data, K_model, green_p);

        rms += pow(K[1] - K_model[1], 2.0);
        ++nrms;

#if 0
        fprintf(stdout, "%f %f %f %f %f %f %f\n",
                90.0 - theta * 180.0 / M_PI,
                B[0],
                B[1],
                B[2],
                B_model[0],
                B_model[1],
                B_model[2]);
#elif 1
        fprintf(stdout, "%f %e %e\n",
                90.0 - theta * 180.0 / M_PI,
                K[1],
                K_model[1]);
#endif
      }

    if (nrms > 0)
      rms = sqrt(rms / nrms);

    gsl_vector_free(res);
  }

  green_free(green_p);
  gsl_matrix_free(X);
  gsl_vector_free(rhs);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  return s;
}

/* test cf 1D SECS */
int
test_cf(const double theta0, secs1d_workspace *w)
{
  int s = 0;
  const double eps = 1.0e-4;
  const double r = R_EARTH_KM + 450.0;
  const double dtheta = 5.0 * M_PI / 180.0;
  const double dphi = 5.0 * M_PI / 180.0;
  double theta, phi;
  const double nmax = 80;
  green_workspace *green_p = green_alloc(nmax, nmax, R_EARTH_KM);
  const size_t n = 100000;
  const size_t p = green_nnm(green_p);
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_vector *rhs = gsl_vector_alloc(n);
  gsl_vector *c = gsl_vector_alloc(p);
  gsl_matrix *cov = gsl_matrix_alloc(p, p);
  size_t rowidx = 0;

  for (theta = eps; theta <= M_PI - eps; theta += dtheta)
    {
      for (phi = 0.0; phi <= 2.0 * M_PI; phi += dphi)
        {
          double B[3];
          gsl_vector_view vx = gsl_matrix_row(X, rowidx);
          gsl_vector_view vy = gsl_matrix_row(X, rowidx + 1);
          gsl_vector_view vz = gsl_matrix_row(X, rowidx + 2);

          /* compute SECS Green's functions */
          secs1d_green_cf(r, theta, theta0, B, w);

          /* compute internal potential field Green's functions */
          green_calc_int(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data, green_p);

          gsl_vector_set(rhs, rowidx, B[0]);
          gsl_vector_set(rhs, rowidx + 1, B[1]);
          gsl_vector_set(rhs, rowidx + 2, B[2]);

          rowidx += 3;
        }
    }

  {
    gsl_matrix_view m = gsl_matrix_submatrix(X, 0, 0, rowidx, p);
    gsl_vector_view b = gsl_vector_subvector(rhs, 0, rowidx);
    gsl_vector *res = gsl_vector_alloc(rowidx);
    int rank = 0;
    struct timeval tv0, tv1;
    double rms = 0.0;
    size_t nrms = 0;

    fprintf(stderr, "test_cf: size of LS system: (%zu,%zu)\n", rowidx, p);

    fprintf(stderr, "test_cf: solving LS system...");

    gettimeofday(&tv0, NULL);
    lapack_lls2(&m.matrix, &b.vector, c, &rank);
    gettimeofday(&tv1, NULL);
    gsl_multifit_linear_residuals(&m.matrix, &b.vector, c, res);

    fprintf(stderr, "done (%g seconds, rank = %d, residual = %.12e)\n",
            time_diff(tv0, tv1), rank, gsl_blas_dnrm2(res));

    phi = 1.0;
    for (theta = eps; theta <= M_PI - eps; theta += dtheta)
      {
        double B[3], B_model[3], K[3], K_model[3];
        gsl_vector_view vx = gsl_matrix_row(X, 0);
        gsl_vector_view vy = gsl_matrix_row(X, 1);
        gsl_vector_view vz = gsl_matrix_row(X, 2);

        /* compute SECS Green's functions */
        secs1d_green_cf(r, theta, theta0, B, w);

        /* compute internal potential field Green's functions */
        green_calc_int(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data, green_p);

        gsl_blas_ddot(&vx.vector, c, &B_model[0]);
        gsl_blas_ddot(&vy.vector, c, &B_model[1]);
        gsl_blas_ddot(&vz.vector, c, &B_model[2]);

        /* compute current densities */
        secs1d_green_cf_J(R_EARTH_KM + 110.0, theta, theta0, K, w);
        green_eval_sheet_int(R_EARTH_KM + 110.0, theta, phi, c->data, K_model, green_p);

        rms += pow(K[1] - K_model[1], 2.0);
        ++nrms;

#if 0
        fprintf(stdout, "%f %f %f %f %f %f %f\n",
                90.0 - theta * 180.0 / M_PI,
                B[0],
                B[1],
                B[2],
                B_model[0],
                B_model[1],
                B_model[2]);
#elif 1
        fprintf(stdout, "%f %e %e %e %e %e\n",
                90.0 - theta * 180.0 / M_PI,
                K[0],
                K[2],
                K_model[0],
                K_model[1],
                K_model[2]);
#endif
      }

    if (nrms > 0)
      rms = sqrt(rms / nrms);

    gsl_vector_free(res);
  }

  green_free(green_p);
  gsl_matrix_free(X);
  gsl_vector_free(rhs);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  return s;
}

int
main()
{
  const double pole_spacing = 1.0;
  secs1d_workspace *w;

  w = secs1d_alloc(SECS1D_FLG_FIT_DF, SECS1D_LMAX, R_EARTH_KM + 110.0, pole_spacing);

  test_df(M_PI / 3.0, w);
  /*test_cf(M_PI / 3.0, w);*/

  secs1d_free(w);

  return 0;
}
