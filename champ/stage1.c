/*
 * stage1.c
 *
 * 1. Read CHAMP CDF data file
 * 2. Compute main, crustal and external fields
 * 3. Create new CDF with main/crustal/external field outputs
 *
 * Usage:
 *
 * ./stage1 <-i champ_cdf_file>
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>

#include <satdata/satdata.h>

#include "chaos.h"
#include "common.h"
#include "msynth.h"
#include "track.h"

satdata_mag *
champ_read_diamag(const char *filename)
{
  FILE *fp = fopen(filename, "r");
  satdata_mag *data = satdata_mag_alloc(86400);
  char buf[4096];
  size_t n = 0;

  data->R = 6371.2;

  while (fgets(buf, 4096, fp) != NULL)
    {
      int c;
      int year, mm, dd, hh, min, sec, msec;
      double ss, glat, glon, alt, alat, qdlat, qdlon, mlt;
      double X, Y, Z, F, F_dia;

      c = sscanf(buf, "%d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                 &year,
                 &mm,
                 &dd,
                 &hh,
                 &min,
                 &ss,
                 &glat,
                 &glon,
                 &alt,
                 &alat,
                 &qdlat,
                 &qdlon,
                 &mlt,
                 &X,
                 &Y,
                 &Z,
                 &F,
                 &F_dia);
      if (c < 18)
        continue;

      sec = (int) ss;
      msec = (int) ((ss - sec) * 1000.0);

      data->t[n] = computeEPOCH(year, mm, dd, hh, min, sec, msec);
      data->latitude[n] = glat;
      data->longitude[n] = glon;
      data->altitude[n] = alt;
      data->r[n] = data->R + alt;
      SATDATA_VEC_X(data->B, n) = X;
      SATDATA_VEC_Y(data->B, n) = Y;
      SATDATA_VEC_Z(data->B, n) = Z;
      data->F[n] = F_dia;

      if (++n > data->ntot)
        {
          fprintf(stderr, "champ_read_diamag: ntot too small\n");
          break;
        }
    }

  data->n = n;

  fclose(fp);

  return data;
}

int
main(int argc, char *argv[])
{
  satdata_mag *data_in = NULL, *data_out;
  struct timeval tv0, tv1;
  int c;
  char *outfile = NULL;
  int down_sample = 1;
  msynth_workspace *msynth_core;
  int use_chaos = 0;

  while ((c = getopt(argc, argv, "i:o:d:y:c")) != (-1))
    {
      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);

            data_in = satdata_champ_mag_read_L3(optarg, NULL);
            if (!data_in)
              {
                fprintf(stderr, "main: error reading %s\n", optarg);
                exit(1);
              }

            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data_in->n,
                    time_diff(tv0, tv1));

            break;

          case 'o':
            outfile = optarg;
            break;

          case 'y':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);

            data_in = champ_read_diamag(optarg);
            if (!data_in)
              {
                fprintf(stderr, "main: error reading %s\n", optarg);
                exit(1);
              }

            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data_in->n,
                    time_diff(tv0, tv1));
            break;

          case 'd':
            down_sample = atoi(optarg);
            break;

          case 'c':
            use_chaos = 1;
            break;

          default:
            break;
        }
    }

  if (!data_in)
    {
      fprintf(stderr, "Usage: %s <-i champ_cdf_file> [-y yunliang_diamag_corr_file] [-o output_cdf_file] [-d down_sample] [-c]\n",
              argv[0]);
      exit(1);
    }

  if (use_chaos)
    chaos_init();

  if (outfile)
    fprintf(stderr, "output file = %s\n", outfile);

  if (data_in->n == 0)
    exit(1);

  data_out = satdata_mag_alloc(data_in->n);
  msynth_core = msynth_swarm_read(MSYNTH_CHAOS_FILE);

  fprintf(stderr, "main: computing magnetic field model (downsampling by %d)...", down_sample);
  gettimeofday(&tv0, NULL);

  if (use_chaos)
    track_synth_chaos(down_sample, data_in, data_out, msynth_core);
  else
    track_synth(down_sample, data_in, data_out, msynth_core);

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, %zu records kept)\n", time_diff(tv0, tv1), data_out->n);

  if (outfile)
    {
      fprintf(stderr, "Writing %s...", outfile);
      gettimeofday(&tv0, NULL);
      satdata_champ_write(outfile, data_out);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%zu records written, %g seconds)\n", data_out->n,
              time_diff(tv0, tv1));
    }

  if (use_chaos)
    chaos_term();

  satdata_mag_free(data_in);
  satdata_mag_free(data_out);

  return 0;
} /* main() */
