/*
 * main.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include <indices/indices.h>

#include <common/common.h>

#include "estist_calc.h"

int
main(int argc, char *argv[])
{
  dst_workspace *dst_p;
  char *outfile = NULL;
  FILE *fp;
  estist_calc_workspace *estist_calc_p;
  double dst, est, ist;
  time_t t0 = 946684800; /* Jan 1 2000 00:00:00 UTC */
  time_t t1 = 978307199; /* Dec 31 2000 23:59:59 UTC */
  time_t t;
  int s = 0;
  int c;
  struct tm tmp;
  char buf[1024];

  while ((c = getopt(argc, argv, "a:b:o:")) != (-1))
    {
      switch (c)
        {
          case 'a':
            t0 = atoi(optarg);
            break;

          case 'b':
            t1 = atoi(optarg);
            break;

          case 'o':
            outfile = optarg;
            break;
        }
    }

  if (!outfile)
    {
      fprintf(stderr, "Usage: %s [-a start_timestamp] [-b end_timestamp] <-o output_file>\n",
              argv[0]);
      exit(1);
    }

  fp = fopen(outfile, "a");
  if (!fp)
    {
      fprintf(stderr, "main: error opening %s: %s\n",
              outfile, strerror(errno));
      exit(1);
    }

  dst_p = dst_alloc(DST_IDX_FILE);
  estist_calc_p = estist_calc_alloc(DST_IDX_FILE);

  gmtime_r(&t0, &tmp);
  asctime_r(&tmp, buf);
  buf[strlen(buf) - 1] = '\0';
  fprintf(stderr, "main: start time: %ld (%s)\n", t0, buf);

  gmtime_r(&t1, &tmp);
  asctime_r(&tmp, buf);
  buf[strlen(buf) - 1] = '\0';
  fprintf(stderr, "main: end time:   %ld (%s)\n", t1, buf);

  for (t = t0; t <= t1; t += 3600)
    {
      /* output status every 5 days */
      if (t %  432000 == 0)
        {
          gmtime_r(&t, &tmp);
          asctime_r(&tmp, buf);
          buf[strlen(buf) - 1] = '\0';
          fprintf(stderr, "main: processing day: %s\n", buf);
        }

      s += estist_calc_get(t, &est, &ist, estist_calc_p);
      s += dst_get(t, &dst, dst_p);

      gsl_test_abs(dst - est - ist, 0.0, 1.0e-12, "t=%ld", t);

      if (s == 0)
        {
          fprintf(fp, "%12f %12f %12f %12f\n",
                  time2fday(t) + 0.5 * 0.041670, /* add 1/2 hr to be consistent with GFZ timestamps */
                  dst,
                  est,
                  ist);
          fflush(fp);
        }
    }

  dst_free(dst_p);
  estist_calc_free(estist_calc_p);

  fclose(fp);

  return 0;
} /* main() */
