/*
 * test.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "pomme.h"

static char tz_str[] = "TZ=GMT";

int
main(int argc, char *argv[])

{
  pomme_workspace *pomme_p;
  double height;
  time_t t;
  double B[4], B_int[4], B_ext[4];
  time_t t0 = 1066629600;
  time_t t1 = 1068447600;
  double lat = 21.3166;
  double lon = -157.9996;
  double theta = M_PI / 2.0 - lat * M_PI / 180.0;
  double phi = lon * M_PI / 180.0;

  putenv(tz_str);

  pomme_p = pomme_alloc_default();

  for (t = t0; t <= t1; t += 3600)
    {
      double B_ext[4];

      pomme_calc_ext(theta, phi, t, 0.0, B_ext, pomme_p);

      printf("%ld %f %f %f\n",
             t,
             B_ext[0] * 1.0e9,
             B_ext[1] * 1.0e9,
             B_ext[2] * 1.0e9);
    }

  pomme_free(pomme_p);

  return 0;
} /* main() */
