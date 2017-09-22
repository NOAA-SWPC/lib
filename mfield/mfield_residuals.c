/*
 * mfield_residuals
 *
 * Compute X,Y,Z,F residuals of a given dataset against
 * a specified core field model
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rstat.h>

#include <msynth/msynth.h>

#include "euler.h"
#include "magdata.h"

#define QDLAT_CUTOFF       (55.0)

#define PRINT_RES          0

void
print_residual_stats(const magdata *data, msynth_workspace *core_p)
{
  size_t i, j;
  msynth_workspace *crust_p = msynth_mf7_read(MSYNTH_MF7_FILE);
  gsl_rstat_workspace *rstat_x = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_y = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_z = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_f = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_lowx = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_highx = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_lowy = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_highy = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_lowz = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_highz = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_lowf = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_highf = gsl_rstat_alloc();

#if PRINT_RES
  i = 1;
  printf("# Field %zu: time (years)\n", i++);
  printf("# Field %zu: QD latitude (degrees)\n", i++);
  printf("# Field %zu: residual B_x (nT)\n", i++);
  printf("# Field %zu: residual B_y (nT)\n", i++);
  printf("# Field %zu: residual B_z (nT)\n", i++);
  printf("# Field %zu: residual F (nT)\n", i++);
  printf("# Field %zu: fit X data\n", i++);
  printf("# Field %zu: fit Y data\n", i++);
  printf("# Field %zu: fit Z data\n", i++);
  printf("# Field %zu: fit F data\n", i++);
#endif

  for (i = 0; i < data->n; ++i)
    {
      double t = satdata_epoch2year(data->t[i]);
      double r = data->r[i];
      double theta = data->theta[i];
      double phi = data->phi[i];
      int fit_scal = MAGDATA_ExistScalar(data->flags[i]) && MAGDATA_FitMF(data->flags[i]);
      int fit_X = MAGDATA_ExistX(data->flags[i]) && MAGDATA_FitMF(data->flags[i]);
      int fit_Y = MAGDATA_ExistY(data->flags[i]) && MAGDATA_FitMF(data->flags[i]);
      int fit_Z = MAGDATA_ExistZ(data->flags[i]) && MAGDATA_FitMF(data->flags[i]);
      double B_core[4], B_crust[4], B_ext[4], B_total[4], B_res[4];
      double B_nec[3];

      if (MAGDATA_Discarded(data->flags[i]))
        continue;

      /* XXX */
      fit_X = fit_Y = fit_Z;

      B_nec[0] = data->Bx_nec[i];
      B_nec[1] = data->By_nec[i];
      B_nec[2] = data->Bz_nec[i];

#if 0
      B_ext[0] = data->Bx_model[i];
      B_ext[1] = data->By_model[i];
      B_ext[2] = data->Bz_model[i];
#else
      B_ext[0] = B_ext[1] = B_ext[2] = 0.0;
#endif

      msynth_eval(t, r, theta, phi, B_core, core_p);
#if 0
      msynth_eval(t, r, theta, phi, B_crust, crust_p);
#else
      B_crust[0] = B_crust[1] = B_crust[2] = 0.0;
#endif

      for (j = 0; j < 3; ++j)
        {
          B_total[j] = B_core[j] + B_crust[j] + B_ext[j];
          B_res[j] = B_nec[j] - B_total[j];
        }

      B_total[3] = gsl_hypot3(B_total[0], B_total[1], B_total[2]);
      B_res[3] = data->F[i] - B_total[3];

      if (fit_X)
        {
          gsl_rstat_add(B_res[0], rstat_x);

          if (fabs(data->qdlat[i]) <= QDLAT_CUTOFF)
            gsl_rstat_add(B_res[0], rstat_lowx);
          else
            gsl_rstat_add(B_res[0], rstat_highx);
        }

      if (fit_Y)
        {
          gsl_rstat_add(B_res[1], rstat_y);

          if (fabs(data->qdlat[i]) <= QDLAT_CUTOFF)
            gsl_rstat_add(B_res[1], rstat_lowy);
          else
            gsl_rstat_add(B_res[1], rstat_highy);
        }

      if (fit_Z)
        {
          gsl_rstat_add(B_res[2], rstat_z);

          if (fabs(data->qdlat[i]) <= QDLAT_CUTOFF)
            gsl_rstat_add(B_res[2], rstat_lowz);
          else
            gsl_rstat_add(B_res[2], rstat_highz);
        }

      if (fit_scal)
        {
          gsl_rstat_add(B_res[3], rstat_f);

          if (fabs(data->qdlat[i]) <= QDLAT_CUTOFF)
            gsl_rstat_add(B_res[3], rstat_lowf);
          else
            gsl_rstat_add(B_res[3], rstat_highf);
        }

#if PRINT_RES
      printf("%f %f %.3f %.3f %.3f %.3f %d %d %d %d\n",
             t,
             data->qdlat[i],
             B_res[0],
             B_res[1],
             B_res[2],
             B_res[3],
             fit_X,
             fit_Y,
             fit_Z,
             fit_scal);
#endif
    }

  fprintf(stderr, "%8s %10s %12s %12s %12s\n",
          "", "N", "mean (nT)", "sigma (nT)", "rms (nT)");

  fprintf(stderr, "%8s %10zu %12.2f %12.2f %12.2f\n",
          "X",
          gsl_rstat_n(rstat_x),
          gsl_rstat_mean(rstat_x),
          gsl_rstat_sd(rstat_x),
          gsl_rstat_rms(rstat_x));

  fprintf(stderr, "%8s %10zu %12.2f %12.2f %12.2f\n",
          "Y",
          gsl_rstat_n(rstat_y),
          gsl_rstat_mean(rstat_y),
          gsl_rstat_sd(rstat_y),
          gsl_rstat_rms(rstat_y));

  fprintf(stderr, "%8s %10zu %12.2f %12.2f %12.2f\n",
          "Z",
          gsl_rstat_n(rstat_z),
          gsl_rstat_mean(rstat_z),
          gsl_rstat_sd(rstat_z),
          gsl_rstat_rms(rstat_z));

  fprintf(stderr, "%8s %10zu %12.2f %12.2f %12.2f\n",
          "F",
          gsl_rstat_n(rstat_f),
          gsl_rstat_mean(rstat_f),
          gsl_rstat_sd(rstat_f),
          gsl_rstat_rms(rstat_f));

  fprintf(stderr, "%8s %10zu %12.2f %12.2f %12.2f\n",
          "low X",
          gsl_rstat_n(rstat_lowx),
          gsl_rstat_mean(rstat_lowx),
          gsl_rstat_sd(rstat_lowx),
          gsl_rstat_rms(rstat_lowx));

  fprintf(stderr, "%8s %10zu %12.2f %12.2f %12.2f\n",
          "high X",
          gsl_rstat_n(rstat_highx),
          gsl_rstat_mean(rstat_highx),
          gsl_rstat_sd(rstat_highx),
          gsl_rstat_rms(rstat_highx));

  fprintf(stderr, "%8s %10zu %12.2f %12.2f %12.2f\n",
          "low Y",
          gsl_rstat_n(rstat_lowy),
          gsl_rstat_mean(rstat_lowy),
          gsl_rstat_sd(rstat_lowy),
          gsl_rstat_rms(rstat_lowy));

  fprintf(stderr, "%8s %10zu %12.2f %12.2f %12.2f\n",
          "high Y",
          gsl_rstat_n(rstat_highy),
          gsl_rstat_mean(rstat_highy),
          gsl_rstat_sd(rstat_highy),
          gsl_rstat_rms(rstat_highy));

  fprintf(stderr, "%8s %10zu %12.2f %12.2f %12.2f\n",
          "low Z",
          gsl_rstat_n(rstat_lowz),
          gsl_rstat_mean(rstat_lowz),
          gsl_rstat_sd(rstat_lowz),
          gsl_rstat_rms(rstat_lowz));

  fprintf(stderr, "%8s %10zu %12.2f %12.2f %12.2f\n",
          "high Z",
          gsl_rstat_n(rstat_highz),
          gsl_rstat_mean(rstat_highz),
          gsl_rstat_sd(rstat_highz),
          gsl_rstat_rms(rstat_highz));

  fprintf(stderr, "%8s %10zu %12.2f %12.2f %12.2f\n",
          "low F",
          gsl_rstat_n(rstat_lowf),
          gsl_rstat_mean(rstat_lowf),
          gsl_rstat_sd(rstat_lowf),
          gsl_rstat_rms(rstat_lowf));

  fprintf(stderr, "%8s %10zu %12.2f %12.2f %12.2f\n",
          "high F",
          gsl_rstat_n(rstat_highf),
          gsl_rstat_mean(rstat_highf),
          gsl_rstat_sd(rstat_highf),
          gsl_rstat_rms(rstat_highf));

  msynth_free(crust_p);
  gsl_rstat_free(rstat_x);
  gsl_rstat_free(rstat_y);
  gsl_rstat_free(rstat_z);
  gsl_rstat_free(rstat_f);
}

int
main(int argc, char *argv[])
{
  euler_workspace *euler_p = NULL;
  magdata *mdata = NULL;
  msynth_workspace *msynth_p = NULL;
  int c;
  double tmin = -1.0;
  double tmax = -1.0;

  while ((c = getopt(argc, argv, "a:b:i:c:e:h")) != (-1))
    {
      switch (c)
        {
          case 'a':
            tmin = atof(optarg);
            break;

          case 'b':
            tmax = atof(optarg);
            break;

          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            mdata = magdata_read(optarg, NULL);
            fprintf(stderr, "done (%zu data read)\n", mdata->n);
            break;

          case 'c':
            fprintf(stderr, "main: reading coefficients from %s...", optarg);
            msynth_p = msynth_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'e':
            fprintf(stderr, "main: reading Euler angles from %s...", optarg);
            euler_p = euler_read(optarg);
            fprintf(stderr, "done (%zu angles read)\n", euler_p->n);
            break;

          case 'h':
            fprintf(stderr, "main: reading coefficients from %s...", MSYNTH_CHAOS_FILE);
            msynth_p = msynth_swarm_read(MSYNTH_CHAOS_FILE);
            fprintf(stderr, "done\n");
            break;

          default:
            break;
        }
    }

  if (!mdata || !msynth_p)
    {
      fprintf(stderr, "%s <-i magdata_file> [-c mysnth_coeff_file] [-h] [-e euler_file] [-a tmin] [-b tmax]\n", argv[0]);
      exit(1);
    }

  msynth_set(1, 15, msynth_p);

  if (euler_p)
    {
      fprintf(stderr, "main: applying Euler angle rotations...");
      euler_magdata_apply(mdata, euler_p);
      fprintf(stderr, "done\n");
    }

  {
    size_t i;

    for (i = 0; i < mdata->n; ++i)
      {
        double t = satdata_epoch2year(mdata->t[i]);

        if (tmin > 0.0 && t < tmin)
          mdata->flags[i] |= MAGDATA_FLG_DISCARD;
        else if (tmax > 0.0 && t > tmax)
          mdata->flags[i] |= MAGDATA_FLG_DISCARD;
      }
  }

  print_residual_stats(mdata, msynth_p);

  msynth_free(msynth_p);

  return 0;
} /* main() */
