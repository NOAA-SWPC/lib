/*
 * print_map.c
 *
 * Print model lat/lon map
 *
 * Usage: ./print_spectrum [args]
 * [-c coef_file]            - ASCII coefficient file
 * [-m mf7_file]
 * [-w wmm_file]
 * [-a arnaud_file]
 * [-h swarm_shc_file]
 * [-p pomme_file]
 * [-i igrf12_mf_candidate]
 * [-g ipgp_file]
 * [-z]
 * [-n nmax]
 * [-t epoch]
 * [-o output_file]
 * [-d]                      - Compute model difference with CHAOS
 *
 * Note: if two coefficient files are specified, the model difference is computed
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include "common.h"
#include "msynth.h"
#include "msynth_grid.h"

#define GRID        1

#if GRID

int
print_map(const char *filename, const double epoch, msynth_workspace *w)
{
  const double r = 6371.2;
  const double c = 3485.0; /* CMB radius */
  const size_t nphi = 512;
  const size_t ntheta = 256;
  const size_t nmin = w->eval_nmin;
  const size_t nmax = w->eval_nmax;
  const size_t epoch_idx = msynth_epoch_idx(epoch, w);
  const double *g = w->c + epoch_idx * w->p;
  const double *dg = g + w->sv_offset;
  const double *ddg = g + w->sa_offset;
  const double theta_min = 0.01;
  const double theta_max = M_PI - 0.01;
  const double theta_step = (theta_max - theta_min) / (ntheta - 1.0);
  msynth_grid_workspace *grid_main_p;
  msynth_grid_workspace *grid_sa_p;
  gsl_matrix *B_main = gsl_matrix_alloc(nphi, 3);
  gsl_matrix *B_sa = gsl_matrix_alloc(nphi, 3);
  double *coef = malloc(w->nnm * sizeof(double));
  double lon, lat;
  FILE *fp;
  size_t i, j;

  fp = fopen(filename, "w");
  if (!fp)
    return -1;

  msynth_gnm(epoch, coef, w);

  grid_main_p = msynth_grid_alloc(nphi, nmin, nmax, g);
  grid_sa_p = msynth_grid_alloc(nphi, nmin, nmax, ddg);

  msynth_grid_init_r(r, grid_main_p);
  msynth_grid_init_r(c, grid_sa_p);

  i = 1;
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: B_x (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_y (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_z (nT)\n", i++);
  fprintf(fp, "# Field %zu: SA B_r at CMB (uT)\n", i++);

  for (j = 0; j < ntheta; ++j)
    {
      double theta = theta_min + j * theta_step;

      msynth_grid_init_theta(theta, grid_main_p);
      msynth_grid_init_theta(theta, grid_sa_p);

      msynth_grid_calc(0, r, theta, B_main, grid_main_p);
      msynth_grid_calc(0, c, theta, B_sa, grid_sa_p);

      /* print values for phi \in [-180:0] */
      for (i = nphi / 2 + 1; i < nphi; ++i)
        {
          double Bx_main = gsl_matrix_get(B_main, i, 0);
          double By_main = gsl_matrix_get(B_main, i, 1);
          double Bz_main = gsl_matrix_get(B_main, i, 2);
          double Bz_sa = gsl_matrix_get(B_sa, i, 2);
          double phi = 2.0 * M_PI * i / (double) nphi;

          fprintf(fp, "%f %f %f %f %f %f\n",
                  wrap180(phi * 180.0 / M_PI),
                  90.0 - theta * 180.0 / M_PI,
                  Bx_main,
                  By_main,
                  Bz_main,
                  -Bz_sa * 1.0e-3);
        }

      /* print values for phi \in [0:180] */
      for (i = 0; i < nphi / 2; ++i)
        {
          double Bx_main = gsl_matrix_get(B_main, i, 0);
          double By_main = gsl_matrix_get(B_main, i, 1);
          double Bz_main = gsl_matrix_get(B_main, i, 2);
          double Bz_sa = gsl_matrix_get(B_sa, i, 2);
          double phi = 2.0 * M_PI * i / (double) nphi;

          fprintf(fp, "%f %f %f %f %f %f\n",
                  wrap180(phi * 180.0 / M_PI),
                  90.0 - theta * 180.0 / M_PI,
                  Bx_main,
                  By_main,
                  Bz_main,
                  -Bz_sa * 1.0e-3);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  free(coef);
  gsl_matrix_free(B_main);
  gsl_matrix_free(B_sa);
  msynth_grid_free(grid_main_p);
  msynth_grid_free(grid_sa_p);

  return 0;
}

#elif 1

int
print_map(const char *filename, const double epoch, msynth_workspace *w)
{
  const double r = 6371.2;
  const double c = 3485.0; /* CMB radius */
  const size_t nphi = 512;
  const size_t ntheta = 256;
  const size_t nmin = w->eval_nmin;
  const size_t nmax = w->eval_nmax;
  const size_t epoch_idx = msynth_epoch_idx(epoch, w);
  const double *g = w->c + epoch_idx * w->p;
  const double *dg = g + w->sv_offset;
  const double *ddg = g + w->sa_offset;
  const double theta_min = 0.01;
  const double theta_max = M_PI - 0.01;
  const double theta_step = (theta_max - theta_min) / (ntheta - 1.0);
  gsl_matrix *B_main = gsl_matrix_alloc(nphi, 3);
  gsl_matrix *B_sa = gsl_matrix_alloc(nphi, 3);
  double *coef = malloc(w->nnm * sizeof(double));
  double lon, lat;
  FILE *fp;
  size_t i, j;

  fp = fopen(filename, "w");
  if (!fp)
    return -1;

  msynth_gnm(epoch, coef, w);

  i = 1;
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: B_x (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_y (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_z (nT)\n", i++);
  fprintf(fp, "# Field %zu: SA B_r at CMB (uT)\n", i++);

  for (j = 0; j < ntheta; ++j)
    {
      double theta = theta_min + j * theta_step;

      /* print values for phi \in [-180:0] */
      for (i = nphi / 2 + 1; i < nphi; ++i)
        {
          double phi = 2.0 * M_PI * i / (double) nphi;
          double B[4], B_mf[3], B_sv[3], B_sa[3];

          msynth_eval(epoch, r, theta, phi, B, w);
          msynth_eval2(epoch, c, theta, phi, B_mf, B_sv, B_sa, w);

          fprintf(fp, "%f %f %f %f %f %f\n",
                  wrap180(phi * 180.0 / M_PI),
                  90.0 - theta * 180.0 / M_PI,
                  B[0],
                  B[1],
                  B[2],
                  -B_sa[2] * 1.0e-3);
        }

      for (i = 0; i < nphi / 2; ++i)
        {
          double phi = 2.0 * M_PI * i / (double) nphi;
          double B[4], B_mf[3], B_sv[3], B_sa[3];

          msynth_eval(epoch, r, theta, phi, B, w);
          msynth_eval2(epoch, c, theta, phi, B_mf, B_sv, B_sa, w);

          fprintf(fp, "%f %f %f %f %f %f\n",
                  wrap180(phi * 180.0 / M_PI),
                  90.0 - theta * 180.0 / M_PI,
                  B[0],
                  B[1],
                  B[2],
                  -B_sa[2] * 1.0e-3);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  free(coef);
  gsl_matrix_free(B_main);
  gsl_matrix_free(B_sa);

  return 0;
}

#endif

int
main(int argc, char *argv[])
{
  msynth_workspace *w = NULL;
  msynth_workspace *msynth1 = NULL;
  msynth_workspace *msynth2 = NULL;
  char *outfile = "map.dat";
  int c;
  size_t nmax = 0;
  double epoch = -1.0;
  int diff = 0;

  while ((c = getopt(argc, argv, "c:m:w:o:a:h:p:n:e:i:g:d")) != (-1))
    {
      switch (c)
        {
          case 'c':
            if (!msynth1)
              msynth1 = msynth_read(optarg);
            else
              msynth2 = msynth_read(optarg);
            break;

          case 'a':
            if (!msynth1)
              msynth1 = msynth_arnaud_read(optarg);
            else
              msynth2 = msynth_arnaud_read(optarg);
            break;

          case 'm':
            if (!msynth1)
              msynth1 = msynth_mf7_read(optarg);
            else
              msynth2 = msynth_mf7_read(optarg);
            break;

          case 'w':
            if (!msynth1)
              msynth1 = msynth_wmm_read(optarg);
            else
              msynth2 = msynth_wmm_read(optarg);
            break;

          case 'h':
            if (!msynth1)
              msynth1 = msynth_swarm_read(optarg);
            else
              msynth2 = msynth_swarm_read(optarg);
            break;

          case 'p':
            if (!msynth1)
              msynth1 = msynth_pomme_read(optarg);
            else
              msynth2 = msynth_pomme_read(optarg);
            break;

          case 'g':
            if (!msynth1)
              msynth1 = msynth_ipgp_read(optarg);
            else
              msynth2 = msynth_ipgp_read(optarg);
            break;

          case 'i':
            if (!msynth1)
              {
                msynth1 = msynth_igrf_read_mf(optarg);
                msynth_igrf_read_sv("cof/IGRF12/SV-2015-2020-E.cof", msynth1);
              }
            else
              {
                msynth2 = msynth_igrf_read_mf(optarg);
                msynth_igrf_read_sv("cof/IGRF12/SV-2015-2020-E.cof", msynth2);
              }
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

          case 'd':
            diff = 1;
            msynth2 = msynth_swarm_read(MSYNTH_CHAOS_FILE);
            break;
        }
    }

  if (!msynth1)
    {
      fprintf(stderr, "Usage: %s [-c coef_file] [-m mf7_file] [-w wmm_file] [-a arnaud_file] [-h swarm_shc_file] [-p pomme_file] [-i igrf12_mf_candidate] [-g ipgp_file] [-z] [-n nmax] [-t epoch] [-o output_file] [-d]\n", argv[0]);
      exit(1);
    }

  if (epoch < 0.0)
    epoch = msynth1->epochs[0];

  if (msynth2)
    {
      /* compute model difference */
      fprintf(stderr, "main: computing coefficient difference...");
      w = msynth_diff(epoch, msynth1, msynth2);
      fprintf(stderr, "done\n");
    }
  else
    w = msynth_copy(msynth1);

  if (nmax != 0)
    msynth_set(1, nmax, w);

  fprintf(stderr, "main: epoch = %g\n", epoch);
  fprintf(stderr, "main: spectrum nmax = %zu\n", w->eval_nmax);

  fprintf(stderr, "main: printing map to %s...", outfile);
  print_map(outfile, epoch, w);
  fprintf(stderr, "done\n");

  /*XXX*/
  {
    const char *sfile = "F15_new.txt";
    msynth_swarm_write(sfile, w);
  }

  msynth_free(w);
  msynth_free(msynth1);

  if (msynth2)
    msynth_free(msynth2);

  return 0;
} /* main() */
