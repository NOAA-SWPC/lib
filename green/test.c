#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

#include "common.h"
#include "green.h"

int
main()
{
  const char *outfile = "mat/green.dat";
  const size_t nmax = 60;
  const size_t mmax = 60;
  const size_t npoints = 500000;
  green_workspace *green_p = green_alloc(nmax, mmax, R_EARTH_KM);
  size_t nnm = green_nnm(green_p);
  double *r = malloc(npoints * sizeof(double));
  double *theta = malloc(npoints * sizeof(double));
  double *phi = malloc(npoints * sizeof(double));
  gsl_matrix *dX = gsl_matrix_alloc(npoints, nnm);
  gsl_matrix *dY = gsl_matrix_alloc(npoints, nnm);
  gsl_matrix *dZ = gsl_matrix_alloc(npoints, nnm);
  gsl_rng *rng_p = gsl_rng_alloc(gsl_rng_default);
  size_t i;
  struct timeval tv0, tv1;
  FILE *fp;

  fprintf(stderr, "main: nmax = %zu\n", nmax);
  fprintf(stderr, "main: mmax = %zu\n", mmax);
  fprintf(stderr, "main: npoints = %zu\n", npoints);

  fprintf(stderr, "main: selecting observation points...");

  for (i = 0; i < npoints; ++i)
    {
      double r_i = gsl_rng_uniform(rng_p) * 6500.0;
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
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fp = fopen(outfile, "w");
  fprintf(stderr, "writing Green's functions to %s...", outfile);
  gettimeofday(&tv0, NULL);

  gsl_matrix_fwrite(fp, dX);
  gsl_matrix_fwrite(fp, dY);
  gsl_matrix_fwrite(fp, dZ);

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
  fclose(fp);

  fp = fopen(outfile, "r");
  fprintf(stderr, "reading Green's functions from %s...", outfile);
  gettimeofday(&tv0, NULL);

  gsl_matrix_fread(fp, dX);
  gsl_matrix_fread(fp, dY);
  gsl_matrix_fread(fp, dZ);

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
  fclose(fp);

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
