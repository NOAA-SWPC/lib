/*
 * test.c
 *
 * Test the Est/Ist module by reading all data files, testing a few
 * specific known values, and checking for missing data
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include "common.h"
#include "estist2.h"

int
main(int argc, char *argv[])
{
  int s = 0;
  estist2_workspace *estist2_workspace_p;
  time_t t0_start = 1041379200; /* Jan 1 2003 00:00:00 */
  time_t t0_2011 = 1293840000;
  time_t t0_2013 = 1356998400;
  time_t t0_2014 = 1388534400;
  time_t t;
  double est, ist;

  estist2_workspace_p = estist2_alloc(ESTIST2_IDX_FILE);
  if (!estist2_workspace_p)
    exit(1);

  /* advance by 1 hour */
  for (t = t0_start; t < t0_2014; t += 3600)
    {
      s += estist2_get(t, &est, &ist, estist2_workspace_p);
      printf("%ld %f %f %f\n", t, time2fday(t), est, ist);
    }

  estist2_free(estist2_workspace_p);

  fprintf(stderr, "%s: %d errors\n", argv[0], s);

  return s;
} /* main() */
