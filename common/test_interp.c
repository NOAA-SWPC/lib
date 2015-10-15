#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "interp.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

int
main()
{
  double gridpts[4][2] = { {0.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}, {1.0, 1.0}};
  double z[4] = { 0.0, 1.0, 1.0, 0.5 };
  double x, y;

  for (x = 0.0; x < 1.0; x += 0.01)
    {
      for (y = 0.0; y < 1.0; y += 0.01)
        {
          double f = interp2d(x, y, gridpts, z);
          double expected;

          expected = z[0]*(1.-x)*(1.-y) +
                     z[1]*x*(1.-y) +
                     z[2]*(1.-x)*y +
                     z[3]*x*y;

          gsl_test_rel(f, expected, 1.0e-10, "x = %f, y = %f\n", x, y);
        }
    }

  return 0;
}
