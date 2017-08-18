/*
 * grid.c
 *
 * Compute grid of field model
 *
 * Usage: ./grid [args]
 * [-c coef_file]            - ASCII coefficient file
 * [-m mf7_file]
 * [-w wmm_file]
 * [-a arnaud_file]
 * [-b NGDC720_file]
 * [-f EMM_crust_file]
 * [-h swarm_file]
 * [-p pomme_file]
 * [-i igrf12_mf_candidate]
 * [-g ipgp_file]
 * [-l bggm_file]
 * [-n nmax]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>

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

Inputs: naive    - 1 = use naive approach, 0 = use fft approach
        nr       - number of radial grid points
        ntheta   - number of theta grid points
        nphi     - number of phi grid points
        t        - (output) time elapsed
        grid_X   - (output) if not NULL, store computed X grid here
        grid_Y   - (output) if not NULL, store computed Y grid here
        grid_Z   - (output) if not NULL, store computed Z grid here
        msynth_p - msynth workspace
        grid_p   - msynth grid workspace

Return: success/error
*/

int
build_grid(const int naive, const size_t nr, const size_t ntheta, const size_t nphi,
           double *t, double *grid_X, double *grid_Y, double *grid_Z,
           msynth_workspace *msynth_p, msynth_grid_workspace *grid_p)
{
  gsl_matrix *B = gsl_matrix_alloc(nphi, 3);
  const double rmin = R_EARTH_KM - 10.0;
  const double rmax = R_EARTH_KM + 10.0;
  const double rstep = (rmax - rmin) / nr;
  const double theta_min = 0.01;
  const double theta_max = M_PI - 0.01;
  const double theta_step = (theta_max - theta_min) / (ntheta - 1.0);
  size_t i, j, k;
  struct timeval tv0, tv1;
  double *theta_array = malloc(ntheta * sizeof(double));

  for (j = 0; j < ntheta; ++j)
    theta_array[j] = theta_min + j * theta_step;

  fprintf(stderr, "\n");

  if (naive)
    {
      double *g = msynth_p->c;
      double *dg = g + msynth_p->sv_offset;
      double *ddg = g + msynth_p->sa_offset;

      gettimeofday(&tv0, NULL);

      for (j = 0; j < ntheta; ++j)
        {
          double progress = (double)j / (double)ntheta;
          double theta = theta_array[j];

          /* print progress bar */
          progress_bar(stderr, progress, 80);

          /* compute Legendre functions for this theta */
          msynth_green_init_ALF(theta, msynth_p);

          for (k = 0; k < nphi; ++k)
            {
              double phi = 2.0 * M_PI * k / (double) nphi;

              /* compute phi arrays */
              msynth_green_init_phi(phi, msynth_p);

              for (i = 0; i < nr; ++i)
                {
                  double r = rmin + i * rstep;
                  double B_naive[4];

                  /* compute internal Green's functions */
                  msynth_green_calc_int(r, theta, msynth_p);

                  msynth_eval_sum(2017.0, 2017.0, g, dg, ddg, B_naive, msynth_p);

                  if (grid_X)
                    {
                      grid_X[CIDX3(i,nr,j,ntheta,k,nphi)] = B_naive[0];
                      grid_Y[CIDX3(i,nr,j,ntheta,k,nphi)] = B_naive[1];
                      grid_Z[CIDX3(i,nr,j,ntheta,k,nphi)] = B_naive[2];
                    }
                }
            }
        }

      gettimeofday(&tv1, NULL);
      *t = time_diff(tv0, tv1);
    }
  else /* fast method */
    {
      gettimeofday(&tv0, NULL);

      for (j = 0; j < ntheta / 2; ++j)
        {
          double progress = 2.0 * (double)j / (double)ntheta;
          double theta = theta_array[j];
          double cost = cos(theta);

          /* print progress bar */
          progress_bar(stderr, progress, 80);

          /* compute associated Legendre functions */
          gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, grid_p->nmax, cost,
                                          msynth_p->Plm, msynth_p->dPlm);

          for (i = 0; i < nr; ++i)
            {
              double r = rmin + i * rstep;

              msynth_grid_calc(r, theta, B, grid_p);

              if (grid_X)
                {
                  for (k = 0; k < nphi; ++k)
                    {
                      grid_X[CIDX3(i,nr,j,ntheta,k,nphi)] = gsl_matrix_get(B, k, 0);
                      grid_Y[CIDX3(i,nr,j,ntheta,k,nphi)] = gsl_matrix_get(B, k, 1);
                      grid_Z[CIDX3(i,nr,j,ntheta,k,nphi)] = gsl_matrix_get(B, k, 2);
                    }
                }

              /*
               * now, using the same Plm/dPlm arrays, calculate the mirror point about
               * the equator (r, pi - theta, phi_k)
               */
              msynth_grid_calc2(r, theta, B, grid_p);

              if (grid_X)
                {
                  for (k = 0; k < nphi; ++k)
                    {
                      grid_X[CIDX3(i,nr,ntheta - j - 1,ntheta,k,nphi)] = gsl_matrix_get(B, k, 0);
                      grid_Y[CIDX3(i,nr,ntheta - j - 1,ntheta,k,nphi)] = gsl_matrix_get(B, k, 1);
                      grid_Z[CIDX3(i,nr,ntheta - j - 1,ntheta,k,nphi)] = gsl_matrix_get(B, k, 2);
                    }
                }
            }
        }

      gettimeofday(&tv1, NULL);
      *t = time_diff(tv0, tv1);
    }

  gsl_matrix_free(B);
  free(theta_array);

  return 0;
}

int
main(int argc, char *argv[])
{
  msynth_workspace *w = NULL;
  msynth_grid_workspace *grid_p;
  char *outfile = NULL;
  FILE *fp_out = NULL;
  int c;
  size_t nmax = 0;
  size_t nphi = 2048;
  size_t nr = 7;
  size_t ntheta = 2;
  double epoch = -1.0;
  size_t i;

  while ((c = getopt(argc, argv, "a:b:c:e:f:g:l:m:w:o:h:p:n:")) != (-1))
    {
      switch (c)
        {
          case 'c':
            w = msynth_read(optarg);
            break;

          case 'a':
            w = msynth_arnaud_read(optarg);
            break;

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

          case 'w':
            w = msynth_wmm_read(optarg);
            break;

          case 'h':
            w = msynth_swarm_read(optarg);
            break;

          case 'p':
            w = msynth_pomme_read(optarg);
            break;

          case 'g':
            w = msynth_ipgp_read(optarg);
            break;

          case 'n':
            nmax = (size_t) atoi(optarg);
            break;

          case 'e':
            epoch = atof(optarg);
            break;

          case 'o':
            outfile = optarg;
            break;
        }
    }

  if (!w)
    {
      fprintf(stderr, "Usage: %s [-c coef_file] [-m mf7_file] [-l bggm_file] [-w wmm_file] [-a arnaud_file] [-b NGDC720_file] [-f EMM_crust_file] [-h swarm_shc_file] [-p pomme_file] [-g ipgp_file] [-n nmax] [-e epoch] [-o output_file]\n", argv[0]);
      exit(1);
    }

  if (epoch < 0.0)
    epoch = w->epochs[0];

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

#if 0

  {
    double t_fft;
    int naive = 0;

    nr = 100;
    ntheta = 8192;
    nphi = 16384;

    grid_p = msynth_grid_alloc(nphi, w);

    fprintf(stderr, "main: building %zu-by-%zu-by-%zu grid [%s method]...",
            nr, ntheta, nphi,
            (naive == 1) ? "naive" : "FFT");
    build_grid(naive, nr, ntheta, nphi, &t_fft, NULL, NULL, NULL, w, grid_p);
    fprintf(stderr, "done (%g seconds)\n", t_fft);
  }

#else

  nr = 1;

  /* loop over phi in the outer loop */
  for (i = 11; i <= 16; ++i)
    {
      double t_fft, t_naive;

      nphi = 1 << i;
      grid_p = msynth_grid_alloc(nphi, w);

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
          build_grid(0, nr, ntheta, nphi, &t_fft, grid_fft_X, grid_fft_Y, grid_fft_Z, w, grid_p);
          fprintf(stderr, "done (%g seconds)\n", t_fft);

          fprintf(stderr, "main: building %zu-by-%zu-by-%zu grid [naive method]...", nr, ntheta, nphi);
          build_grid(1, nr, ntheta, nphi, &t_naive, grid_naive_X, grid_naive_Y, grid_naive_Z, w, grid_p);
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

      msynth_grid_free(grid_p);
    }

#endif

  msynth_free(w);

  if (fp_out)
    fclose(fp_out);

  return 0;
}
