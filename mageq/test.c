/*
 * test.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "mageq.h"

int
main()
{
  mageq_workspace *w;
  double phi;
  double r = 6371.2 + 110.0;
  double lat, angle;
  double t = 2014.0;
  size_t i;

  w = mageq_alloc();

  i = 1;
  printf("# Field %zu: longitude (deg)\n", i++);
  printf("# Field %zu: geocentric latitude (deg)\n", i++);
  printf("# Field %zu: mageq angle (deg)\n", i++);

  for (phi = -M_PI; phi < M_PI; phi += 0.1)
    {
      lat = mageq_calc(phi, r, t, w);
      angle = mageq_angle(phi, r, t, w);
      printf("%f %f %f\n",
             phi * 180.0 / M_PI,
             lat * 180.0 / M_PI,
             angle * 180.0 / M_PI);
    }

  mageq_free(w);

  return 0;
}
