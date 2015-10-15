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

#include "geo.h"
#include "msynth.h"

int
main(int argc, char *argv[])
{
  msynth_workspace *w = NULL;
  int c;
  double t, r, theta, phi;
  double B[4];

  while ((c = getopt(argc, argv, "c:w:a:h:i:")) != (-1))
    {
      switch (c)
        {
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
            w = msynth_chaos_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'i':
            fprintf(stderr, "main: reading IGRF coefficient file %s...", optarg);
            w = msynth_igrf_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'w':
            fprintf(stderr, "main: reading WMM coefficient file %s...", optarg);
            w = msynth_wmm_read(optarg);
            fprintf(stderr, "done\n");
            break;
        }
    }

  if (!w)
    {
      fprintf(stderr, "Usage: %s [-c coef_file] [-w wmm_file] [-a arnaud_file] [-h chaos_file] [-i igrf_file]\n", argv[0]);
      exit(1);
    }

  {
    double latd = 0.0 * M_PI / 180.0;
    double h = 400.0;
    double latc; /* geocentric latitude */

    geodetic2geo(latd, h, &latc, &r);
    theta = M_PI / 2.0 - latc;
  }

  t = 2015.2;
  phi = -50.0 * M_PI / 180.0;

  fprintf(stderr, "geocentric lat = %f [deg]\n", 90.0 - theta * 180.0 / M_PI);
  fprintf(stderr, "geocentric lon = %f [deg]\n", phi * 180.0 / M_PI);
  fprintf(stderr, "geocentric r   = %f [km]\n", r);

  msynth_eval(t, r, theta, phi, B, w);

  fprintf(stderr, "X = %f [nT]\n", B[0]);
  fprintf(stderr, "Y = %f [nT]\n", B[1]);
  fprintf(stderr, "Z = %f [nT]\n", B[2]);
  fprintf(stderr, "F = %f [nT]\n", B[3]);

  msynth_free(w);

  return 0;
} /* main() */
