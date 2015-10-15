/*
 * test_ut1.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "ut1.h"

int
main(int argc, char *argv[])
{
  ut1_workspace *w;
  time_t t0 = 946684800;
  time_t t1 = 1356998400;
  time_t t;

  w = ut1_alloc(UT1_DATA_FILE);

  for (t = t0; t < t1; t += 86400)
    {
      double dut1;
      int s = ut1_get(t, &dut1, w);

      if (s == 0)
        printf("%ld %f\n", t, dut1);
    }

  ut1_free(w);

  return 0;
} /* main() */
