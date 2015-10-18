/*
 * test.c
 * Patrick Alken
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "pomme.h"

static char tz_str[] = "TZ=GMT";

int
main(int argc, char *argv[])

{
  pomme_workspace *pomme_p;
  double lat,
         lon,
         theta;
  int year,
      month,
      day;
  double hour; /* local time hour */
  double height;
  time_t t;
  struct tm tm_p;
  double B[4], B_int[4], B_ext[4];

  lat = 12.53 * M_PI / 180.0;
  lon = 285.983 * M_PI / 180.0;
  year = 2013;
  month = 6;
  day = 3;
  hour = 12.16;

  theta = M_PI / 2.0 - lat;

  /* convert to UT */
  hour = hour - lon * 180.0 / M_PI / 15.0;

  putenv(tz_str);

  tm_p.tm_sec = 0;
  tm_p.tm_min = (int) ((hour - (int)hour) * 60.0);
  tm_p.tm_hour = (int) hour;
  tm_p.tm_mday = day;
  tm_p.tm_year = year - 1900;
  tm_p.tm_mon = month - 1;
  tm_p.tm_isdst = 0;

  t = mktime(&tm_p);

  printf("fday = %f\n", time2fday(t));

  height = 65.0;

  pomme_p = pomme_alloc_default();
  pomme_set_deg(16, pomme_p);

  pomme_calc(theta, lon, t, height, B, pomme_p);

  pomme_calc_int(theta, lon, t, height, B_int, pomme_p);
  pomme_calc_ext(theta, lon, t, height, B_ext, pomme_p);

  pomme_free(pomme_p);

  printf("calc 1: Bx = %.12e, By = %.12e, Bz = %.12e, Bf = %.12e\n",
         B[0], B[1], B[2], B[3]);

  printf("calc 2: Bx = %.12e, By = %.12e, Bz = %.12e, Bf = NC\n",
         B_int[0] + B_ext[0], B_int[1] + B_ext[1], B_int[2] + B_ext[2]);

  return 0;
} /* main() */
