#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "lunar.h"

int
main()
{
  int year, month, day;
  lunardata data;
  time_t t;
  struct tm tms;
  char *tz_str = "TZ=GMT";
  double phase;

  year = 2008;
  month = 1;
  day = 22;

  tms.tm_sec = 0;
  tms.tm_min = 0;
  tms.tm_hour = 12;
  tms.tm_mday = day;
  tms.tm_mon = month - 1;
  tms.tm_year = year - 1900;
  tms.tm_isdst = 0;

  putenv(tz_str);
  t = mktime(&tms);

  lunarcalc(year, month, day, &data);

  phase = lunarphase(t);
  printf("p = %f\n", phase * 180.0 / M_PI);

  phase = time2lunar(t);
  printf("p2 = %f\n", phase * 180.0 / M_PI);

  return 0;
}
