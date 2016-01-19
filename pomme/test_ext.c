/*
 * test.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "pomme.h"

static char tz_str[] = "TZ=GMT";

int
main(int argc, char *argv[])

{
  pomme_workspace *pomme_p;
  double lat,
         lon;
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
  year = 2012;
  month = 10;
  day = 3;
  hour = 12.16;

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
  t = 1134086400;
  /*t = 1369935449;*/

  printf("fday = %f\n", time2fday(t));

  height = 355.1;

  pomme_p = pomme_alloc_default();
  pomme_set_deg(16, pomme_p);

  /*lon = 60.569813255592663;*/
  lon = -60.8;

  /*for (lon = 0.0; lon < 360.0; lon += 5.0)*/
  {
  for (lat = -90.0; lat < 90.0; lat += 2.0)
    {
      double dF_ext = 0.0;
      size_t k;
      double theta = M_PI / 2.0 - lat * M_PI / 180.0;
      double phi = lon * M_PI / 180.0;

      pomme_calc(theta, phi, t, height, B, pomme_p);

      pomme_calc_int(theta, phi, t, height, B_int, pomme_p);
      pomme_calc_ext(theta, phi, t, height, B_ext, pomme_p);

      for (k = 0; k < 3; ++k)
        dF_ext += B_ext[k] * B_int[k] / B_int[3];

      printf("%f %f %.12e %.12e\n",
             lat,
             lon,
             B_int[3] * 1.0e9,
             dF_ext * 1.0e9);
    }
    printf("\n");
  }


  pomme_free(pomme_p);

  return 0;
} /* main() */
