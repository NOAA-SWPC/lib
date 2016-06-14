/*
 * test.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include <indices/indices.h>

#include "common.h"
#include "estist_calc.h"

int
main()
{
  dst_workspace *dst_p;
  estist_calc_workspace *estist_calc_p;
  double dst, est, ist;
  time_t t0, t1, t;
  int s = 0;

  dst_p = dst_alloc(DST_IDX_FILE);
  estist_calc_p = estist_calc_alloc(DST_IDX_FILE);

  /* start/end dates */
  t0 = 1112076000; /* Mar 29 06:00:00 2005 */
  t1 = 1112659200; /* Apr 5 00:00:00 2005 */

  for (t = t0; t <= t1; t += 3600)
    {
      s += estist_calc_get(t, &est, &ist, estist_calc_p);
      s += dst_get(t, &dst, dst_p);

      gsl_test_abs(dst - est - ist, 0.0, 1.0e-12, "t=%ld", t);

      printf("%f %f %f %f\n",
             time2fday(t) + 0.5 * 0.041670, /* add 1/2 hr to be consistent with GFZ timestamps */
             dst,
             est,
             ist);
      fflush(stdout);
    }

  dst_free(dst_p);
  estist_calc_free(estist_calc_p);

  return 0;
} /* main() */
