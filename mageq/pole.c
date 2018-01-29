/*
 * pole.c
 *
 * Calculate position of magnetic poles
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

#include <gsl/gsl_math.h>

#include "magpole.h"

int
main(int argc, char *argv[])
{
  magpole_workspace *w;
  size_t N = 100000;
  double r = 6371.2;
  double t = 2015.0;
  int c;

  while ((c = getopt(argc, argv, "t:n:")) != (-1))
    {
      switch (c)
        {
          case 't':
            t = atof(optarg);
            break;

          case 'n':
            N = (size_t) atol(optarg);
            break;

          case '?':
          default:
            printf("usage: %s [-t epoch (decimal year)]\n", argv[0]);
        }
    }

  w = magpole_alloc(N);

  /*XXX*/
  msynth_set(1, 1, w->igrf_workspace_p);

  magpole_calc(t, r, w);

  fprintf(stderr, "Epoch         : %g\n", t);
  fprintf(stderr, "Samples       : %zu\n", N);
  fprintf(stderr, "Pole latitude : %.4f [deg]\n", 90.0 - w->theta_pole * 180.0 / M_PI);
  fprintf(stderr, "Pole longitude: %.4f [deg]\n", w->phi_pole * 180.0 / M_PI);

  magpole_free(w);

  return 0;
}
