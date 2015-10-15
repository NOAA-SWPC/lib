#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "lunar.h"

int
main()
{
  time_t t;
  double longitude;

  longitude = -76.87 * M_PI / 180.0;

  for (t = 946684800; t < 1046684800; t += 60)
    {
      double nu = lunarphase(t);

      if (nu < 0.0005)
        {
          double tau = lunartime(t, longitude);
          double lt = get_localtime(t, longitude);

          printf("%f %f\n", lt, tau);
        }
    }

  return 0;
}
