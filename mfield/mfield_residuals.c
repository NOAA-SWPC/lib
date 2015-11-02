/*
 * mfield_residuals
 *
 * Compute X,Y,Z,F residuals of a given dataset against
 * DMSP-MAG-1 model
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <gsl/gsl_math.h>

#include <satdata/satdata.h>

#include "mfield_eval.h"

int
main(int argc, char *argv[])
{
  char *infile = NULL;
  mfield_eval_workspace *w;
  satdata_mag *data;
  size_t i;
  int c;

  while ((c = getopt(argc, argv, "i:")) != (-1))
    {
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
      fprintf(stderr, "%s <-i dmsp_idx_file>\n", argv[0]);
      exit(1);
    }

  w = mfield_eval_read("dmsp_mag_1.txt");
  if (!w)
    exit(1);

  fprintf(stderr, "main: reading %s...", infile);
  data = satdata_dmsp_read_idx(infile, 1);
  if (!data)
    exit(1);
  fprintf(stderr, "done (%zu data read)\n", data->n);

  {
    const double thresh_low = 20.0;
    const double thresh_high = 60.0;
    double rms_low[4], rms_mid[4], rms_high[4];
    size_t nlow = 0, nmid = 0, nhigh = 0;

    for (i = 0; i < 4; ++i)
      {
        rms_low[i] = 0.0;
        rms_mid[i] = 0.0;
        rms_high[i] = 0.0;
      }

    for (i = 0; i < data->n; ++i)
      {
        double t = satdata_epoch2year(data->t[i]);
        double r = data->R + data->altitude[i];
        double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
        double phi = data->longitude[i] * M_PI / 180.0;
        double B[4], B_sat[4];

        /* get satellite measurement */
        B_sat[0] = SATDATA_VEC_X(data->B, i);
        B_sat[1] = SATDATA_VEC_Y(data->B, i);
        B_sat[2] = SATDATA_VEC_Z(data->B, i);
        B_sat[3] = data->F[i];

        mfield_eval(t, r, theta, phi, B, w);

        /* F17 has 3 or so screwy B_z measurements which distort the rms */
        if (fabs(B_sat[2] - B[2]) > 10000.0 ||
            fabs(B_sat[3] - B[3]) > 500.0)
          continue;

        if (fabs(data->latitude[i]) < thresh_low)
          {
            rms_low[0] += pow(B_sat[0] - B[0], 2.0);
            rms_low[1] += pow(B_sat[1] - B[1], 2.0);
            rms_low[2] += pow(B_sat[2] - B[2], 2.0);
            rms_low[3] += pow(B_sat[3] - B[3], 2.0);
            ++nlow;
          }
        else if (fabs(data->latitude[i]) > thresh_high)
          {
            rms_high[0] += pow(B_sat[0] - B[0], 2.0);
            rms_high[1] += pow(B_sat[1] - B[1], 2.0);
            rms_high[2] += pow(B_sat[2] - B[2], 2.0);
            rms_high[3] += pow(B_sat[3] - B[3], 2.0);
            ++nhigh;
          }
        else
          {
            rms_mid[0] += pow(B_sat[0] - B[0], 2.0);
            rms_mid[1] += pow(B_sat[1] - B[1], 2.0);
            rms_mid[2] += pow(B_sat[2] - B[2], 2.0);
            rms_mid[3] += pow(B_sat[3] - B[3], 2.0);
            ++nmid;
          }
      }
    
    for (i = 0; i < 4; ++i)
      {
        rms_low[i] = sqrt(rms_low[i] / nlow);
        rms_mid[i] = sqrt(rms_mid[i] / nmid);
        rms_high[i] = sqrt(rms_high[i] / nhigh);
      }

    fprintf(stderr, "=== Low Latitude stats (lat below %.2f [deg], %zu points) ===\n",
            thresh_low, nlow);
    fprintf(stderr, "rms X = %.2f [nT]\n", rms_low[0]);
    fprintf(stderr, "rms Y = %.2f [nT]\n", rms_low[1]);
    fprintf(stderr, "rms Z = %.2f [nT]\n", rms_low[2]);
    fprintf(stderr, "rms F = %.2f [nT]\n", rms_low[3]);

    fprintf(stderr, "=== Mid Latitude stats (lat between %.2f and %.2f [deg], %zu points) ===\n",
            thresh_low, thresh_high, nmid);
    fprintf(stderr, "rms X = %.2f [nT]\n", rms_mid[0]);
    fprintf(stderr, "rms Y = %.2f [nT]\n", rms_mid[1]);
    fprintf(stderr, "rms Z = %.2f [nT]\n", rms_mid[2]);
    fprintf(stderr, "rms F = %.2f [nT]\n", rms_mid[3]);

    fprintf(stderr, "=== High Latitude stats (lat above %.2f [deg], %zu points) ===\n",
            thresh_high, nhigh);
    fprintf(stderr, "rms X = %.2f [nT]\n", rms_high[0]);
    fprintf(stderr, "rms Y = %.2f [nT]\n", rms_high[1]);
    fprintf(stderr, "rms Z = %.2f [nT]\n", rms_high[2]);
    fprintf(stderr, "rms F = %.2f [nT]\n", rms_high[3]);
  }

  mfield_eval_free(w);
  satdata_mag_free(data);

  return 0;
} /* main() */
