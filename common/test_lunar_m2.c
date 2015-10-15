#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "lunar.h"

int
main()
{
  time_t t;
  double lon = 0.0 * M_PI / 180.0;

  /*for (t = 946684800; t < 1046684800; t += 86400)*/
  for (t = -315594000; t < 0; t += 86400)
    {
      double nu = lunarphase(t);
      double m2 = lunarphase_m2(t);
      double tau = lunartime(t, lon);
      double tau2 = lunartime_m2_r(t, lon) * 24.83333 / 2.0 / M_PI;

      printf("%d %f %f\n", t, nu, m2);
    }

  return 0;
}
