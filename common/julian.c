/*
 * julian.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "julian.h"

/*
date2julian()
  Convert a given date to julian day

Inputs: year -
        month - [1,12]
        day   - [1,31]
        hour  - UT hour [0,23]
        min   - UT min [0,59]
        sec   - UT sec [0.0,59.999]

Return: julian day (number of days from 1 Jan 4713 BC)
*/

double
date2julian(const int year, const int month, const int day,
            const int hour, const int min, const double sec)
{
  double jd;

  jd = 367.0 * year -
       floor((7 * (year + floor((month + 9) / 12.0))) * 0.25) +
       floor( 275 * month / 9.0 ) +
       day + 1721013.5 +
       ((sec / 60.0 + min) / 60.0 + hour) / 24.0;

  return jd;
} /* date2julian() */

double
timet2julian(const time_t t)
{
  double jd;

#if 0

  int year, month, day, hour, min;
  double sec;
  struct tm *tmp;

  tmp = gmtime(&t);

  year = tmp->tm_year + 1900;
  month = tmp->tm_mon + 1;
  day = tmp->tm_mday;
  hour = tmp->tm_hour;
  min = tmp->tm_min;
  sec = (double) tmp->tm_sec;

  jd = date2julian(year, month, day, hour, min, sec);

#else

  jd = (t / 86400.0) + 2440587.5;

#endif

  return jd;
} /* timet2julian() */

/*
julian2GMST()
  Convert Julian day to Greenwich Mean Sidereal Time

Inputs: jd - julian day

Return: mean sidereal time in radians

Reference:
http://aa.usno.navy.mil/faq/docs/GAST.php
*/

double
julian2GMST(double jd)
{
  double jdmin, jdmax;
  double jd0 = 0.0;
  double H, D, D0, T;
  double GMST;

  jdmin = floor(jd) - 0.5;
  jdmax = floor(jd) + 0.5;

  if (jd > jdmin)
    jd0 = jdmin;
  if (jd > jdmax)
    jd0 = jdmax;

  H = (jd - jd0) * 24.0;
  D = jd - 2451545.0;
  D0 = jd0 - 2451545.0;
  T = D / 36525.0;

  /* calculate GMST in hours */
  GMST = 6.697374558 +
         0.06570982441908 * D0 +
         1.00273790935 * H +
         0.000026 * T * T;

  /* convert to radians */
  GMST *= M_PI / 12.0;

  return GMST;
} /* julian2GMST() */

/*
julian2GAST()
  Convert Julian day to Greenwich Apparent Sidreal Time

Inputs: jd - julian day

Return: Greenwich apparent sidereal time in radians
*/

double
julian2GAST(double jd)
{
  double T;
  double THETA, THETAm;
  double L, dL, omega;
  double epsm, dpsi, deps;

  T = (jd - 2451545.0) / 36525.0;

  /* mean sidereal time in radians */
  THETAm = julian2GMST(jd);

  /* compute nutations in obliquity and longitude */
  L = 280.4665 + 36000.7698 * T;
  dL = 218.3165 + 481267.8813 * T;
  omega = 125.04452 - 1934.136261 * T;

  /* convert to radians */
  L *= M_PI / 180.0;
  dL *= M_PI / 180.0;
  omega *= M_PI / 180.0;

  /* compute epsilon_m in degrees */
  epsm = 23.439291 -
         0.0130111 * T -
         1.64e-7 * T * T +
         5.04e-7 * T * T * T;

  /* convert to radians */
  epsm *= M_PI / 180.0;

  dpsi = -17.2 * sin(omega) -
          1.32 * sin(2.0 * L) -
          0.23 * sin(2.0 * dL) +
          0.21 * sin(2.0 * omega);
  deps = 9.2 * cos(omega) +
         0.57 * cos(2.0 * L) +
         0.10 * cos(2.0 * dL) -
         0.09 * cos(2.0 * omega);

  /* convert from arc-seconds to radians */
  dpsi *= M_PI / 180.0 / 3600.0;
  deps *= M_PI / 180.0 / 3600.0;

  THETA = THETAm + dpsi * cos(epsm + deps);

  return wrap2pi(THETA);
} /* julian2GAST() */
