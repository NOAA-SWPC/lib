/*
 * storm_rms.c
 *
 * Compute rms between Weimer and ground observatory data
 * for various storms
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include <common/bsearch.h>

#include "grobs.h"

/* IAGA data has bad baselines for DED in 2015, so compute mean
 * of each component and subtract from data */
#define SUBTRACT_MEAN           1

/*
compute_rms

Inputs: t0          - starting time
        t1          - ending time
        iaga_data   - IAGA measurement
        weimer_data - Weimer model
        rms         - (output) rms between IAGA data and Weimer prediction
                      rms[0] = rms X [nT]
                      rms[1] = rms Y [nT]
                      rms[2] = rms Z [nT]
                      rms[3] = rms D [deg]
                      rms[4] = rms I [deg]
                      rms[5] = rms F [nT]
        rms_data    - (output) rms of IAGA data minus baseline
*/

int
compute_rms(const char *filename, const time_t t0, const time_t t1,
            const grobs_data *iaga_data,
            const grobs_data *weimer_data, double rms[6], double rms_data[6], size_t * nrms)
{
  size_t i, j;
  size_t n = 0;
  size_t ncomp = 6;
  double mean_iaga[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double mean_weimer[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  FILE *fp = fopen(filename, "w");

  for (i = 0; i < ncomp; ++i)
    {
      rms[i] = 0.0;
      rms_data[i] = 0.0;
    }

#if SUBTRACT_MEAN
  /* XXX IAGA data has bad baselines for DED in 2015, so compute mean
   * of each component and subtract
   */
  {
    for (i = 0; i < iaga_data->n; ++i)
      {
        if (iaga_data->t[i] < t0 || iaga_data->t[i] > t1)
          continue;

        j = bsearch_timet(weimer_data->t, iaga_data->t[i], 0, weimer_data->n);

        if (weimer_data->t[j] != iaga_data->t[i])
          continue;

        mean_iaga[0] += iaga_data->X[i];
        mean_iaga[1] += iaga_data->Y[i];
        mean_iaga[2] += iaga_data->Z[i];
        mean_iaga[3] += iaga_data->D[i];
        mean_iaga[4] += iaga_data->I[i];
        mean_iaga[5] += gsl_hypot3(iaga_data->X[i], iaga_data->Y[i], iaga_data->Z[i]);

        mean_weimer[0] += weimer_data->X[j];
        mean_weimer[1] += weimer_data->Y[j];
        mean_weimer[2] += weimer_data->Z[j];
        mean_weimer[3] += weimer_data->D[j];
        mean_weimer[4] += weimer_data->I[j];
        mean_weimer[5] += gsl_hypot3(weimer_data->X[j], weimer_data->Y[j], weimer_data->Z[j]);

        ++n;

        if (n >= 720)
          break;
      }

    for (i = 0; i < ncomp; ++i)
      {
        mean_iaga[i] /= (double) n;
        mean_weimer[i] /= (double) n;
      }

    n = 0;
  }
#endif

  i = 1;
  fprintf(fp, "# Field %zu: timestamp\n", i++);
  fprintf(fp, "# Field %zu: IAGA X measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: IAGA Y measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: IAGA Z measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: IAGA D measurement (deg)\n", i++);
  fprintf(fp, "# Field %zu: IAGA I measurement (deg)\n", i++);
  fprintf(fp, "# Field %zu: IAGA F measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: Weimer X model (nT)\n", i++);
  fprintf(fp, "# Field %zu: Weimer Y model (nT)\n", i++);
  fprintf(fp, "# Field %zu: Weimer Z model (nT)\n", i++);
  fprintf(fp, "# Field %zu: Weimer D model (deg)\n", i++);
  fprintf(fp, "# Field %zu: Weimer I model (deg)\n", i++);
  fprintf(fp, "# Field %zu: Weimer F model (nT)\n", i++);

  for (i = 0; i < iaga_data->n; ++i)
    {
      double F_iaga, F_weimer;

      if (iaga_data->t[i] < t0 || iaga_data->t[i] > t1)
        continue;

      j = bsearch_timet(weimer_data->t, iaga_data->t[i], 0, weimer_data->n);

      if (weimer_data->t[j] != iaga_data->t[i])
        continue;

      F_iaga = gsl_hypot3(iaga_data->X[i], iaga_data->Y[i], iaga_data->Z[i]);
      F_weimer = gsl_hypot3(weimer_data->X[j], weimer_data->Y[j], weimer_data->Z[j]);

      /* XXX */
      iaga_data->X[i] -= mean_iaga[0];
      iaga_data->Y[i] -= mean_iaga[1];
      iaga_data->Z[i] -= mean_iaga[2];
      iaga_data->D[i] -= mean_iaga[3];
      iaga_data->I[i] -= mean_iaga[4];
      F_iaga -= mean_iaga[5];

      weimer_data->X[j] -= mean_weimer[0];
      weimer_data->Y[j] -= mean_weimer[1];
      weimer_data->Z[j] -= mean_weimer[2];
      weimer_data->D[j] -= mean_weimer[3];
      weimer_data->I[j] -= mean_weimer[4];
      F_weimer -= mean_weimer[5];

      rms[0] += pow(iaga_data->X[i] - weimer_data->X[j], 2.0);
      rms[1] += pow(iaga_data->Y[i] - weimer_data->Y[j], 2.0);
      rms[2] += pow(iaga_data->Z[i] - weimer_data->Z[j], 2.0);
      rms[3] += pow(iaga_data->D[i] - weimer_data->D[j], 2.0);
      rms[4] += pow(iaga_data->I[i] - weimer_data->I[j], 2.0);
      rms[5] += pow(F_iaga - F_weimer, 2.0);

      rms_data[0] += pow(iaga_data->X[i], 2.0);
      rms_data[1] += pow(iaga_data->Y[i], 2.0);
      rms_data[2] += pow(iaga_data->Z[i], 2.0);
      rms_data[3] += pow(iaga_data->D[i], 2.0);
      rms_data[4] += pow(iaga_data->I[i], 2.0);
      rms_data[5] += pow(F_iaga, 2.0);

      fprintf(fp, "%ld %f %f %f %f %f %f %f %f %f %f %f %f\n",
              iaga_data->t[i],
              iaga_data->X[i],
              iaga_data->Y[i],
              iaga_data->Z[i],
              iaga_data->D[i],
              iaga_data->I[i],
              F_iaga,
              weimer_data->X[j],
              weimer_data->Y[j],
              weimer_data->Z[j],
              weimer_data->D[j],
              weimer_data->I[j],
              F_weimer);

      ++n;
    }

  if (n > 0)
    {
      size_t i;

      for (i = 0; i < ncomp; ++i)
        {
          rms[i] = sqrt(rms[i] / (double)n);
          rms_data[i] = sqrt(rms_data[i] / (double)n);
        }
    }

  *nrms = n;

  fclose(fp);

  return 0;
}

int
main(int argc, char *argv[])
{
  int c;
  char *outfile = "output.txt";
  grobs_data *iaga_data = NULL;
  grobs_data *weimer_data = NULL;
  time_t t_storm = 0;
  time_t t0, t1;
  double rms[6], rms_data[6];
  size_t nrms;

  while ((c = getopt(argc, argv, "i:t:w:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            iaga_data = grobs_iaga_read(optarg, NULL);
            break;

          case 't':
            t_storm = atoi(optarg);
            break;

          case 'w':
            weimer_data = grobs_iaga_read(optarg, NULL);
            break;

          case '?':
          default:
            printf("usage: %s [options]\n", argv[0]);
            exit(1);
            break;
        }
    }

  if (iaga_data == NULL || weimer_data == NULL)
    {
      fprintf(stderr, "usage: %s <-i iaga_data_file> <-w weimer_data_file> <-t t_storm>\n",
              argv[0]);
      exit(1);
    }

  /* use +/- 2 day window for rms */
  t0 = t_storm - 2*86400;
  t1 = t_storm + 2*86400;

  compute_rms(outfile, t0, t1, iaga_data, weimer_data, rms, rms_data, &nrms);
  fprintf(stderr, "main: data written to %s\n", outfile);

  fprintf(stderr, "%s", ctime(&t_storm));
  fprintf(stderr, "nrms  = %zu\n", nrms);
  fprintf(stderr, "rms X = %.2f (%.2f, delta = %.2f) [nT]\n", rms[0], rms_data[0], rms_data[0] - rms[0]);
  fprintf(stderr, "rms Y = %.2f (%.2f, delta = %.2f) [nT]\n", rms[1], rms_data[1], rms_data[1] - rms[1]);
  fprintf(stderr, "rms Z = %.2f (%.2f, delta = %.2f) [nT]\n", rms[2], rms_data[2], rms_data[2] - rms[2]);
  fprintf(stderr, "rms D = %.4f (%.4f, delta = %.4f) [deg]\n", rms[3], rms_data[3], rms_data[3] - rms[3]);
  fprintf(stderr, "rms I = %.4f (%.4f, delta = %.4f) [deg]\n", rms[4], rms_data[4], rms_data[4] - rms[4]);
  fprintf(stderr, "rms F = %.2f (%.2f, delta = %.2f) [nT]\n", rms[5], rms_data[5], rms_data[5] - rms[5]);

  return 0;
}
