/*
 * test.c
 *
 * This program tests for orthogonality of the Green's functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "common.h"
#include "green.h"

int
main()
{
  const size_t nmax = 60;
  const size_t mmax = 60;
  const size_t npoints = 10000;
  const double rmin = R_EARTH_KM;
  const double rmax = rmin + 500.0;
  green_workspace *green_p = green_alloc(nmax, mmax, R_EARTH_KM);
  size_t nnm = green_nnm(green_p);
  double *r = malloc(npoints * sizeof(double));
  double *theta = malloc(npoints * sizeof(double));
  double *phi = malloc(npoints * sizeof(double));
  gsl_matrix *dX = gsl_matrix_alloc(npoints, nnm);
  gsl_matrix *dY = gsl_matrix_alloc(npoints, nnm);
  gsl_matrix *dZ = gsl_matrix_alloc(npoints, nnm);
  gsl_matrix *dX2 = gsl_matrix_alloc(npoints, nnm);
  gsl_matrix *dY2 = gsl_matrix_alloc(npoints, nnm);
  gsl_matrix *dZ2 = gsl_matrix_alloc(npoints, nnm);
  gsl_rng *rng_p = gsl_rng_alloc(gsl_rng_default);
  size_t i;
  struct timeval tv0, tv1;

  fprintf(stderr, "main: nmax = %zu\n", nmax);
  fprintf(stderr, "main: mmax = %zu\n", mmax);
  fprintf(stderr, "main: npoints = %zu\n", npoints);

  fprintf(stderr, "main: selecting observation points...");

  for (i = 0; i < npoints; ++i)
    {
      double r_i = gsl_rng_uniform(rng_p) * (rmax - rmin) + rmin;
      double theta_i = gsl_rng_uniform(rng_p) * M_PI;
      double phi_i = gsl_rng_uniform(rng_p) * 2.0 * M_PI;

      r[i] = r_i;
      theta[i] = theta_i;
      phi[i] = phi_i;
    }

  fprintf(stderr, "done\n");

  fprintf(stderr, "main: computing Green's functions...");
  gettimeofday(&tv0, NULL);

  for (i = 0; i < npoints; ++i)
    {
      gsl_vector_view X = gsl_matrix_row(dX, i);
      gsl_vector_view Y = gsl_matrix_row(dY, i);
      gsl_vector_view Z = gsl_matrix_row(dZ, i);

      green_calc_int(r[i], theta[i], phi[i],
                     X.vector.data, Y.vector.data, Z.vector.data,
                     green_p);

      {
        double r_i = gsl_rng_uniform(rng_p) * (rmax - rmin) + rmin;
        double theta_i = gsl_rng_uniform(rng_p) * M_PI;
        double phi_i = gsl_rng_uniform(rng_p) * 2.0 * M_PI;
        X = gsl_matrix_row(dX2, i);
        Y = gsl_matrix_row(dY2, i);
        Z = gsl_matrix_row(dZ2, i);

        green_calc_int(r_i, theta_i, phi[i],
                       X.vector.data, Y.vector.data, Z.vector.data,
                       green_p);
      }
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: testing orthogonality of Green's functions...");
  gettimeofday(&tv0, NULL);

  for (i = 0; i < npoints; ++i)
    {
      const double eps = 1.0e-10;
      const double eps2 = 1.0e-8;
      gsl_vector_view X = gsl_matrix_row(dX, i);
      gsl_vector_view Y = gsl_matrix_row(dY, i);
      gsl_vector_view Z = gsl_matrix_row(dZ, i);
      gsl_vector_view X2 = gsl_matrix_row(dX2, i);
      gsl_vector_view Y2 = gsl_matrix_row(dY2, i);
      gsl_vector_view Z2 = gsl_matrix_row(dZ2, i);
      double XY, XZ, YZ;

      if (theta[i] < 1.0e-2 || theta[i] > M_PI - 1.0e-2)
        continue; /* precision is lost near the poles */

      gsl_blas_ddot(&X.vector, &Y2.vector, &XY);
      gsl_blas_ddot(&X.vector, &Z2.vector, &XZ);
      gsl_blas_ddot(&Y.vector, &Z2.vector, &YZ);

      if (fabs(XY) > eps || fabs(XZ) > eps2 || fabs(YZ) > eps)
        {
          printf("%f %f %f %e %e %e\n",
                 r[i],
                 90.0 - theta[i] * 180.0 / M_PI,
                 phi[i] * 180.0 / M_PI,
                 XY,
                 XZ,
                 YZ);
        }
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  green_free(green_p);
  free(r);
  free(theta);
  free(phi);
  gsl_matrix_free(dX);
  gsl_matrix_free(dY);
  gsl_matrix_free(dZ);
  gsl_rng_free(rng_p);

  return 0;
}
