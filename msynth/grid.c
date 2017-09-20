/*
 * grid.c
 *
 * Compute grid of field model
 *
 * Usage: ./grid [args]
 * [-m mf7_file]
 * [-b NGDC720_file]
 * [-f EMM_crust_file]
 * [-h swarm_file]
 * [-l bggm_file]
 * [-n nmax]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>
#include <omp.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_test.h>

#include "common.h"
#include "msynth.h"
#include "msynth_grid.h"

/*
build_grid()

Inputs: naive   - 1 = use naive approach, 0 = use fft approach
        nr      - number of radial grid points
        ntheta  - number of theta grid points
        nphi    - number of phi grid points
        t       - (output) time elapsed
        grid_X  - (output) if not NULL, store computed X grid here
        grid_Y  - (output) if not NULL, store computed Y grid here
        grid_Z  - (output) if not NULL, store computed Z grid here
        w       - msynth workspace
        grid_p  - msynth grid workspace

Return: success/error
*/

int
build_grid(const int naive, const size_t nr, const size_t ntheta, const size_t nphi,
           double *t, double *grid_X, double *grid_Y, double *grid_Z,
           msynth_workspace *w)
{
  const int max_threads = omp_get_max_threads();
  const double rmin = R_EARTH_KM - 10.0;
  const double rmax = R_EARTH_KM + 10.0;
  const double rstep = (rmax - rmin) / nr;
  const double theta_min = 0.01;
  const double theta_max = M_PI - 0.01;
  const double theta_step = (theta_max - theta_min) / (ntheta - 1.0);
  gsl_matrix **B = malloc(max_threads * sizeof(gsl_matrix *));
  msynth_workspace **msynth_p = malloc(max_threads * sizeof(msynth_workspace *));
  msynth_grid_workspace **grid_p = malloc(max_threads * sizeof(msynth_grid_workspace *));
  size_t j;
  struct timeval tv0, tv1;
  double *theta_array = malloc(ntheta * sizeof(double));
  size_t *omp_ntheta = calloc(max_threads, sizeof(size_t));

  for (j = 0; j < (size_t) max_threads; ++j)
    {
      B[j] = gsl_matrix_alloc(nphi, 3);
      msynth_p[j] = msynth_copy(w);
      grid_p[j] = msynth_grid_alloc(nphi, w->eval_nmin, w->eval_nmax, w->c);
    }

  for (j = 0; j < ntheta; ++j)
    theta_array[j] = theta_min + j * theta_step;

  fprintf(stderr, "\n");

  if (naive)
    {
      double *g = w->c;
      double *dg = g + w->sv_offset;
      double *ddg = g + w->sa_offset;

      gettimeofday(&tv0, NULL);

#pragma omp parallel for private(j)
      for (j = 0; j < ntheta; ++j)
        {
          int thread_id = omp_get_thread_num();
          double progress = (double)j / (double)ntheta;
          double theta = theta_array[j];
          size_t i, k;

          /* compute Legendre functions for this theta */
          msynth_green_init_ALF(theta, msynth_p[thread_id]);

          for (k = 0; k < nphi; ++k)
            {
              double phi = 2.0 * M_PI * k / (double) nphi;

              /* compute phi arrays */
              msynth_green_init_phi(phi, msynth_p[thread_id]);

              for (i = 0; i < nr; ++i)
                {
                  double r = rmin + i * rstep;
                  double B_naive[4];

                  /* compute internal Green's functions */
                  msynth_green_calc_int(r, theta, msynth_p[thread_id]);

                  msynth_eval_sum(2017.0, 2017.0, g, dg, ddg, B_naive, msynth_p[thread_id]);

                  if (grid_X)
                    {
                      grid_X[CIDX3(i,nr,j,ntheta,k,nphi)] = B_naive[0];
                      grid_Y[CIDX3(i,nr,j,ntheta,k,nphi)] = B_naive[1];
                      grid_Z[CIDX3(i,nr,j,ntheta,k,nphi)] = B_naive[2];
                    }
                }
            }

          if (thread_id == 0)
            {
              double progress = 0.0;
              size_t l;

              for (l = 0; l < (size_t) max_threads; ++l)
                progress += omp_ntheta[l];

              progress /= (double) ntheta;

              /* print progress bar */
              progress_bar(stderr, progress, 80);
            }

          omp_ntheta[thread_id]++;
        }

      gettimeofday(&tv1, NULL);
      *t = time_diff(tv0, tv1);

      progress_bar(stderr, 1.0, 80);
      fprintf(stderr, "\n");
    }
  else /* fast method */
    {
      gettimeofday(&tv0, NULL);

#pragma omp parallel for private(j)
      for (j = 0; j < ntheta / 2; ++j)
        {
          int thread_id = omp_get_thread_num();
          double progress = 2.0 * (double)j / (double)ntheta;
          double theta = theta_array[j];
          size_t i, k;

          /* compute associated Legendre functions */
          msynth_grid_init_theta(theta, grid_p[thread_id]);

          for (i = 0; i < nr; ++i)
            {
              double r = rmin + i * rstep;

              /* precompute radial terms */
              msynth_grid_init_r(r, grid_p[thread_id]);

              /* calculate grid */
              msynth_grid_calc(0, r, theta, B[thread_id], grid_p[thread_id]);

              if (grid_X)
                {
                  for (k = 0; k < nphi; ++k)
                    {
                      grid_X[CIDX3(i,nr,j,ntheta,k,nphi)] = gsl_matrix_get(B[thread_id], k, 0);
                      grid_Y[CIDX3(i,nr,j,ntheta,k,nphi)] = gsl_matrix_get(B[thread_id], k, 1);
                      grid_Z[CIDX3(i,nr,j,ntheta,k,nphi)] = gsl_matrix_get(B[thread_id], k, 2);
                    }
                }

              /*
               * now, using the same Plm/dPlm arrays, calculate the mirror point about
               * the equator (r, pi - theta, phi_k)
               */
              msynth_grid_calc(1, r, theta, B[thread_id], grid_p[thread_id]);

              if (grid_X)
                {
                  for (k = 0; k < nphi; ++k)
                    {
                      grid_X[CIDX3(i,nr,ntheta - j - 1,ntheta,k,nphi)] = gsl_matrix_get(B[thread_id], k, 0);
                      grid_Y[CIDX3(i,nr,ntheta - j - 1,ntheta,k,nphi)] = gsl_matrix_get(B[thread_id], k, 1);
                      grid_Z[CIDX3(i,nr,ntheta - j - 1,ntheta,k,nphi)] = gsl_matrix_get(B[thread_id], k, 2);
                    }
                }
            }

          if (thread_id == 0)
            {
              double progress = 0.0;
              size_t l;

              for (l = 0; l < (size_t) max_threads; ++l)
                progress += omp_ntheta[l];

              progress *= 2.0 / (double) ntheta;

              /* print progress bar */
              progress_bar(stderr, progress, 80);
            }

          omp_ntheta[thread_id]++;
        }

      gettimeofday(&tv1, NULL);
      *t = time_diff(tv0, tv1);

      progress_bar(stderr, 1.0, 80);
      fprintf(stderr, "\n");
    }

  for (j = 0; j < (size_t) max_threads; ++j)
    {
      gsl_matrix_free(B[j]);
      msynth_free(msynth_p[j]);
      msynth_grid_free(grid_p[j]);
    }

  free(B);
  free(msynth_p);
  free(grid_p);
  free(theta_array);

  return 0;
}

int
main(int argc, char *argv[])
{
  msynth_workspace *w = NULL;
  char *outfile = NULL;
  FILE *fp_out = NULL;
  int c;
  size_t nmax = 0;
  size_t nphi = 2048;
  size_t nr = 7;
  size_t ntheta = 2;
  size_t i;

  while ((c = getopt(argc, argv, "b:f:l:m:o:n:")) != (-1))
    {
      switch (c)
        {
          case 'b':
            w = msynth_ngdc720_read(optarg);
            break;

          case 'f':
            w = msynth_emm_read(optarg);
            break;

          case 'l':
            w = msynth_bggm_read(optarg);
            break;

          case 'm':
            w = msynth_mf7_read(optarg);
            break;

          case 'n':
            nmax = (size_t) atoi(optarg);
            break;

          case 'o':
            outfile = optarg;
            break;
        }
    }

  if (!w)
    {
      fprintf(stderr, "Usage: %s [-m mf7_file] [-l bggm_file] [-b NGDC720_file] [-f EMM_crust_file] [-n nmax] [-o output_file]\n", argv[0]);
      exit(1);
    }

  if (nmax != 0)
    msynth_set(1, nmax, w);

  fprintf(stderr, "main: spectrum nmin = %zu\n", w->eval_nmin);
  fprintf(stderr, "main: spectrum nmax = %zu\n", w->eval_nmax);

  if (outfile)
    {
      fp_out = fopen(outfile, "w");
      fprintf(stderr, "main: output file   = %s\n", outfile);

      i = 1;
      fprintf(fp_out, "# Field %zu: total grid points\n", i++);
      fprintf(fp_out, "# Field %zu: radial grid points\n", i++);
      fprintf(fp_out, "# Field %zu: theta grid points\n", i++);
      fprintf(fp_out, "# Field %zu: phi grid points\n", i++);
      fprintf(fp_out, "# Field %zu: time elapsed for FFT approach (sec)\n", i++);
      fprintf(fp_out, "# Field %zu: time elapsed for naive approach (sec)\n", i++);
      fprintf(fp_out, "# Field %zu: X residual between FFT and naive grids (nT)\n", i++);
      fprintf(fp_out, "# Field %zu: Y residual between FFT and naive grids (nT)\n", i++);
      fprintf(fp_out, "# Field %zu: Z residual between FFT and naive grids (nT)\n", i++);
    }

#if 1

  {
    double t_fft;
    int naive = 1;

    nr = 100;
    ntheta = 8192;
    nphi = 16384;

    fprintf(stderr, "main: building %zu-by-%zu-by-%zu grid [%s method]...",
            nr, ntheta, nphi,
            (naive == 1) ? "naive" : "FFT");
    build_grid(naive, nr, ntheta, nphi, &t_fft, NULL, NULL, NULL, w);
    fprintf(stderr, "done (%g seconds)\n", t_fft);
  }

#else

  nr = 1;

  /* loop over phi in the outer loop */
  for (i = 11; i <= 16; ++i)
    {
      double t_fft, t_naive;

      nphi = 1 << i;

      for (ntheta = 10; ntheta <= 100; ntheta += 10)
        {
          size_t ntot = nr * ntheta * nphi;
          double *grid_fft_X = malloc(ntot * sizeof(double));
          double *grid_fft_Y = malloc(ntot * sizeof(double));
          double *grid_fft_Z = malloc(ntot * sizeof(double));
          double *grid_naive_X = malloc(ntot * sizeof(double));
          double *grid_naive_Y = malloc(ntot * sizeof(double));
          double *grid_naive_Z = malloc(ntot * sizeof(double));

          fprintf(stderr, "main: building %zu-by-%zu-by-%zu grid [FFT method]...", nr, ntheta, nphi);
          build_grid(0, nr, ntheta, nphi, &t_fft, grid_fft_X, grid_fft_Y, grid_fft_Z, w );
          fprintf(stderr, "done (%g seconds)\n", t_fft);

          fprintf(stderr, "main: building %zu-by-%zu-by-%zu grid [naive method]...", nr, ntheta, nphi);
          build_grid(1, nr, ntheta, nphi, &t_naive, grid_naive_X, grid_naive_Y, grid_naive_Z, w);
          fprintf(stderr, "done (%g seconds)\n", t_naive);

          if (fp_out)
            {
              gsl_vector_view v1, v2;
              double res[3];

              /* compute residuals between grids */

              v1 = gsl_vector_view_array(grid_fft_X, ntot);
              v2 = gsl_vector_view_array(grid_naive_X, ntot);
              gsl_vector_sub(&v1.vector, &v2.vector);
              res[0] = gsl_blas_dnrm2(&v1.vector);

              v1 = gsl_vector_view_array(grid_fft_Y, ntot);
              v2 = gsl_vector_view_array(grid_naive_Y, ntot);
              gsl_vector_sub(&v1.vector, &v2.vector);
              res[1] = gsl_blas_dnrm2(&v1.vector);

              v1 = gsl_vector_view_array(grid_fft_Z, ntot);
              v2 = gsl_vector_view_array(grid_naive_Z, ntot);
              gsl_vector_sub(&v1.vector, &v2.vector);
              res[2] = gsl_blas_dnrm2(&v1.vector);

              fprintf(fp_out, "%zu %zu %zu %zu %g %g %.12e %.12e %.12e\n",
                      ntot,
                      nr,
                      ntheta,
                      nphi,
                      t_fft,
                      t_naive,
                      res[0],
                      res[1],
                      res[2]);
              fflush(fp_out);
            }

          free(grid_fft_X);
          free(grid_fft_Y);
          free(grid_fft_Z);
          free(grid_naive_X);
          free(grid_naive_Y);
          free(grid_naive_Z);
        }
    }

#endif

  msynth_free(w);

  if (fp_out)
    fclose(fp_out);

  return 0;
}
