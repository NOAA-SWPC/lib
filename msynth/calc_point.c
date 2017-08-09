/*
 * calc_point.c
 *
 * Synthesize magnetic field for a specified point
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "geo.h"
#include "msynth.h"

int
main(int argc, char *argv[])
{
  msynth_workspace *w = NULL;
  time_t unix_time;
  double t, r, theta, phi;
  double latd, lond, altd, latc;
  double B[4];

  /* default values */
  latd = 0.0;
  lond = 0.0;
  altd = 0.0;
  unix_time = 1420070400;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "latitude", required_argument, NULL, 'A' },
          { "longitude", required_argument, NULL, 'B' },
          { "altitude", required_argument, NULL, 'C' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "A:B:C:c:w:a:h:i:l:t:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'A':
            latd = atof(optarg);
            break;

          case 'B':
            lond = atof(optarg);
            break;

          case 'C':
            altd = atof(optarg);
            break;

          case 'c':
            fprintf(stderr, "main: reading coefficient file %s...", optarg);
            w = msynth_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'a':
            fprintf(stderr, "main: reading Arnaud coefficient file %s...", optarg);
            w = msynth_arnaud_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'h':
            fprintf(stderr, "main: reading CHAOS coefficient file %s...", optarg);
            w = msynth_swarm_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'i':
            fprintf(stderr, "main: reading IGRF coefficient file %s...", optarg);
            w = msynth_igrf_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'l':
            fprintf(stderr, "main: reading BGGM coefficient file %s...", optarg);
            w = msynth_bggm_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'w':
            fprintf(stderr, "main: reading WMM coefficient file %s...", optarg);
            w = msynth_wmm_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 't':
            unix_time = atoi(optarg);
            break;
        }
    }

  if (!w)
    {
      fprintf(stderr, "Usage: %s [-c coef_file] [-w wmm_file] [-a arnaud_file] [-h swarm_shc_file] [-i igrf_file] [-l bggm_file] [--latitude geodetic_lat] [--longitude longitude] [--height height]\n", argv[0]);
      exit(1);
    }

  geodetic2geo(latd * M_PI / 180.0, altd, &latc, &r);
  theta = M_PI / 2.0 - latc;
  phi = lond * M_PI / 180.0;

  putenv("TZ=GMT");
  t = get_year(unix_time);

  fprintf(stderr, "date = %f / %s", t, ctime(&unix_time));

  fprintf(stderr, "geodetic lat = %f [deg]\n", latd);
  fprintf(stderr, "geodetic lon = %f [deg]\n", lond);
  fprintf(stderr, "geodetic alt = %f [km]\n", altd);

  fprintf(stderr, "geocentric lat = %f [deg]\n", latc * 180.0 / M_PI);
  fprintf(stderr, "geocentric lon = %f [deg]\n", lond);
  fprintf(stderr, "geocentric r   = %f [km]\n", r);

  msynth_eval(t, r, theta, phi, B, w);

  fprintf(stderr, "X = %f [nT]\n", B[0]);
  fprintf(stderr, "Y = %f [nT]\n", B[1]);
  fprintf(stderr, "Z = %f [nT]\n", B[2]);
  fprintf(stderr, "F = %f [nT]\n", B[3]);

  msynth_free(w);

  return 0;
}
