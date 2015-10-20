/*
 * test.c
 * Patrick Alken
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "hwm.h"

static char tz_str[] = "TZ=GMT";

int
main(int argc, char *argv[])

{
  hwm_workspace *hwm_p;
  double lat,
         lon;
  double theta;
  int year,
      month,
      day;
  double hour; /* local time hour */
  double height, hstop, hstep;
  time_t t;
  size_t i;
  size_t nalt;
  struct tm tm_p;
  double merid[512], zonal[512];

  lat = -11.95 * M_PI / 180.0;
  lon = -76.87 * M_PI / 180.0;
  year = 2000;
  month = 8;
  day = 3;
  hour = 12.665199;

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

  t = 1144339320;

  height = 80.0;
  hstop = 500.0;
  nalt = 420;

  hstep = (hstop - height) / (nalt - 1);

  hwm_p = hwm_alloc(F107_IDX_FILE);

  hwm_f107_override(175.0, 175.0, hwm_p);
  /*hwm_set_error_scale(2.0, hwm_p);*/

  theta = M_PI / 2.0 - lat;

  hwm_call(theta,
           lon,
           t,
           height,
           hstep,
           nalt,
           merid,
           zonal,
           hwm_p);

  hwm_free(hwm_p);

  for (i = 0; i < nalt; ++i)
    {
      double h = height + i * hstep;

      printf("%f %e %e\n", h, merid[i], zonal[i]);
    }

  return 0;
} /* main() */
