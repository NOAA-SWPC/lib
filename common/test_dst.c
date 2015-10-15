#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include "common.h"
#include "dst.h"

int
main(int argc, char *argv[])
{
  int s;
  dst_workspace *dst_workspace_p;
  double dst, est, ist;
  time_t t;
  time_t t0 = 820479600;
  time_t t_end = 1293865200;
  double fday = 0.0;

  if (argc > 1)
    fday = atof(argv[1]);

  fprintf(stderr, "reading data file...");
  dst_workspace_p = dst_alloc(DST_DATA_FILE);
  fprintf(stderr, "done\n");

  t = fday2timet(fday);
  dst_get(t, &dst, dst_workspace_p);
  est_get(t, &est, dst_workspace_p);
  ist_get(t, &ist, dst_workspace_p);

  fprintf(stderr, "dst = %f\n", dst);
  fprintf(stderr, "est = %f\n", est);
  fprintf(stderr, "ist = %f\n", ist);

  for (t = t0; t < t_end; t += 1)
    {
      s += dst_get(t, &dst, dst_workspace_p);
      s += est_get(t, &est, dst_workspace_p);
      s += ist_get(t, &ist, dst_workspace_p);
      if (s)
        printf("here\n");

      if (fabs(dst) < 1.0e-10)
        gsl_test_abs(est + ist, dst, 1.0e-10, "t = %ld\n", t);
      else
        gsl_test_rel(est + ist, dst, 1.0e-10, "t = %ld\n", t);
    }

  dst_free(dst_workspace_p);

  fprintf(stderr, "%d errors\n", s);

  exit (gsl_test_summary());
}
