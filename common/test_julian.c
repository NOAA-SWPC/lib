#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "julian.h"

int putenv(char *string);

int
main(int argc, char *argv[])
{
  time_t t, t0;
  struct tm tm_p;
  double jd0, GAST0;
  double hh, s;
  int h, m;

  tm_p.tm_sec = 0;
  tm_p.tm_min = 0;
  tm_p.tm_hour = 0;
  tm_p.tm_mday = 1;
  tm_p.tm_mon = 0;
  tm_p.tm_year = 2000 - 1900;
  tm_p.tm_isdst = 0;

  putenv("TZ=GMT");
  t0 = mktime(&tm_p);

  /* compute julian day */
  jd0 = (t0 / 86400.0) + 2440587.5;

  if (argc > 1)
    jd0 = atof(argv[1]);

  GAST0 = julian2GAST(jd0);

  /* convert to hours, min, sec */
  hh = GAST0 * 12.0 / M_PI;
  h = (int) hh;
  m = (int) ((hh - h) * 60.0);
  s = (hh - h) * 3600.0 - m * 60.0;

  fprintf(stderr, "t = %ld\n", t0);
  fprintf(stderr, "jd = %7.12f\n", jd0);
  fprintf(stderr, "GAST = %1.12f radians\n", GAST0);
  fprintf(stderr, "GAST = %d hours, %d minutes, %f seconds\n", h, m, s);

  for (t = t0; t < t0 + 86400; ++t)
    {
      double jd = (t / 86400.0) + 2440587.5;
      double GAST = julian2GAST(jd);
      double ut = get_ut(t);

      printf("%f %1.12f %1.12f\n",
             ut,
             wrap2pi(GAST),
             wrap2pi(GAST - GAST0));
    }

  return 0;
} /* main() */
