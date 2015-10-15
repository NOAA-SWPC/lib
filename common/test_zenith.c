/*
 * test_zenith.c
 * Patrick Alken
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "solarpos.h"
#include "solpos00.h"

int
main(int argc, char *argv[])
{
  struct posdata pd;
  long retval;
  int c;
  double lat, lon;
  time_t ts;
  struct tm tms;
  int year, hour, min, sec;
  int mon, mday;
  solarpos_workspace *spw;
  double zenith;
  char *tz_str = "TZ=GMT";

  lat = 45.0;
  lon = -120.0;
  year = 2013;
  hour = 15;
  min = 50;
  sec = 0;
  mon = 8;
  mday = 1;

  spw = solarpos_alloc();

  while ((c = getopt(argc, argv, "y:h:m:s:a:o:n:d:")) != (-1))
    {
      switch (c)
        {
          case 'y':
            year = atoi(optarg);
            break;

          case 'h':
            hour = atoi(optarg);
            break;

          case 'm':
            min = atoi(optarg);
            break;

          case 's':
            sec = atoi(optarg);
            break;

          case 'a':
            lat = atoi(optarg);
            break;

          case 'o':
            lon = atoi(optarg);
            break;

          case 'n':
            mon = atoi(optarg);
            break;

          case 'd':
            mday = atoi(optarg);
            break;

          case '?':
          default:
            printf("usage: %s [-h hour] [-m min] [-s sec] [-n month (1-12)] [-d day of month] [-y year] [-a lat] [-o lon]\n", argv[0]);
        }
    }

  tms.tm_sec = sec;
  tms.tm_min = min;
  tms.tm_hour = hour;
  tms.tm_mday = mday;
  tms.tm_mon = mon - 1;
  tms.tm_year = year - 1900;
  tms.tm_isdst = 0;

  putenv(tz_str);
  ts = mktime(&tms);

  S_init(&pd);

  pd.longitude = (float) lon;
  pd.latitude = (float) lat;
  pd.timezone = 0.0;
  pd.year = year;
  pd.daynum = tms.tm_yday + 1;
  pd.hour = hour;
  pd.minute = min;
  pd.second = sec;

  pd.temp = 10.0;
  pd.press = 1013.0;
  pd.aspect = 180.0;

  retval = S_solpos(&pd);
  S_decode(retval, &pd);

  solarpos_calc_zenith(ts,
                       lat * M_PI / 180.0,
                       lon * M_PI / 180.0,
                       &zenith,
                       spw);
  zenith *= 180.0 / M_PI;

  printf("zenith = %f\n", pd.zenetr);
  printf("zenith_me = %f\n", zenith);
  printf("diff = %e\n", pd.zenetr - zenith);

  solarpos_free(spw);

  return 0;
} /* main() */
