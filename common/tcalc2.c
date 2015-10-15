/* tcalc.c
 * 
 * Copyright (C) 2006, 2007 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "tcalc2.h"

static double tcalc_tzone(double lon);
static double tcalc_geom_mean_long_sun(double jc);
static double tcalc_geom_mean_anomaly_sun(double jc);
static double tcalc_eccentricity_earth_orbit(double jc);
static double tcalc_obliquity_correction(double jc);
static double tcalc_sun_eq_of_center(double jc);
static double tcalc_sun_true_long(double jc);
static double tcalc_sun_apparent_long(double jc);
static double tcalc_eqtime(double jc);
static double tcalc_solar_decl(double jc);
static double tcalc_julian_day(int year, int month, double day);
static double tcalc_julian_cent(double jd);
static int tcalc_is_leap_year(int year);
static void tcalc_get_date(time_t t, int *year, int *month, double *day);

static int months[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

/*
tcalc_sunrise()
  Compute local time of sunrise on 'ts' (in minutes)

Inputs: ts   - timestamp
        lat  - latitude in radians
        long - longitude in radians
*/

double
tcalc_sunrise2(time_t ts, double lat, double lon)

{
  int year, month;
  double day;
  double jd, jc;
  double eqtime, decl;
  double ha; /* hour angle */
  double t;
  double zone; /* utc offset in hours */

  tcalc_get_date(ts, &year, &month, &day);

  jd = tcalc_julian_day(year, month, day);
  jc = tcalc_julian_cent(jd);

  eqtime = tcalc_eqtime(jc);
  decl = tcalc_solar_decl(jc);

  ha = acos(cos(90.833*M_PI/180.0) / cos(lat) / cos(decl) - tan(lat) * tan(decl));
  ha *= 180.0 / M_PI;

  t = 720.0 + 4.0 * (lon*180.0/M_PI - ha) - eqtime;

  /* UTC offset in hours */
  zone = tcalc_tzone(lon);

  t = t - (60 * zone);

  return (t);
} /* tcalc_sunrise() */

/*
tcalc_sunset()
  Compute local time of sunset on 'ts' (in minutes)

Inputs: ts   - timestamp
        lat  - latitude in radians
        long - longitude in radians
*/

double
tcalc_sunset2(time_t ts, double lat, double lon)

{
  int year, month;
  double day;
  double jd, jc;
  double eqtime, decl;
  double ha; /* hour angle */
  double t;
  double zone; /* utc offset in hours */

  tcalc_get_date(ts, &year, &month, &day);

  jd = tcalc_julian_day(year, month, day);
  jc = tcalc_julian_cent(jd);

  eqtime = tcalc_eqtime(jc);
  decl = tcalc_solar_decl(jc);

  ha = -acos(cos(90.833*M_PI/180.0) / cos(lat) / cos(decl) - tan(lat) * tan(decl));
  ha *= 180.0 / M_PI;

  t = 720.0 + 4.0 * (lon*180.0/M_PI - ha) - eqtime;

  /* UTC offset in hours */
  zone = tcalc_tzone(lon);

  t = t - (60 * zone);

  return (t);
} /* tcalc_sunset() */

/************************************************
 *          INTERNAL ROUTINES                   *
 ************************************************/

static double
tcalc_tzone(double lon)

{
  double lond = lon * 180.0 / M_PI;

  return (lond / 15.0);
} /* tcalc_tzone() */

static double
tcalc_geom_mean_long_sun(double jc)

{
  double L0 = 280.46646 + jc * (36000.76983 + 0.0003032 * jc);

  while (L0 > 360.0)
    L0 -= 360.0;

  while (L0 < 0.0)
    L0 += 360.0;

  return (L0); /* in degrees */
} /* tcalc_geom_mean_long_sun */

static double
tcalc_geom_mean_anomaly_sun(double jc)

{
  /* in degrees */
  return (357.52911 + jc * (35999.05029 - 0.0001537 * jc));
} /* tcalc_geom_mean_anomaly_sun */

static double
tcalc_eccentricity_earth_orbit(double jc)

{
  /* unitless */
  return (0.016708634 - jc * (0.000042037 + 0.0000001267 * jc));
}

static double
tcalc_obliquity_correction(double jc)

{
  double e0; /* mean obliquity of ecliptic */
  double omega;
  double seconds;
  double e;

  seconds = 21.448 - jc * (46.8150 + jc*(0.00059 - jc*(0.001813)));
  e0 = 23.0 + (26.0 + (seconds/60.0))/60.0;

  omega = 125.04 - 1934.136 * jc;

  e = e0 + 0.00256 * cos(omega * M_PI / 180.0);

  return (e); /* in degrees */
} /* tcalc_obliquity_correction() */

static double
tcalc_sun_eq_of_center(double jc)

{
  double m = tcalc_geom_mean_anomaly_sun(jc);
  double mrad = m * M_PI / 180.0;
  double sinm = sin(mrad);
  double sin2m = sin(2.0 * mrad);
  double sin3m = sin(3.0 * mrad);
  double C;
  
  C = sinm * (1.914602 - jc * (0.004817 + 0.000014 * jc)) +
      sin2m * (0.019993 - 0.000101 * jc) + sin3m * 0.000289;

  return (C); /* in degrees */
} /* tcalc_sun_eq_of_center() */

static double
tcalc_sun_true_long(double jc)

{
  double l0 = tcalc_geom_mean_long_sun(jc);
  double c = tcalc_sun_eq_of_center(jc);
  double O = l0 + c;

  return (O); /* in degrees */
} /* tcalc_sun_true_long() */

static double
tcalc_sun_apparent_long(double jc)

{
  double o = tcalc_sun_true_long(jc);
  double omega = 125.04 - 1934.136 * jc;
  double lambda = o - 0.00569 - 0.00478 * sin(omega * M_PI / 180.0);

  return (lambda);
} /* tcalc_sun_apparent_long() */

static double
tcalc_eqtime(double jc)

{
  double epsilon = tcalc_obliquity_correction(jc);
  double l0 = tcalc_geom_mean_long_sun(jc);
  double e = tcalc_eccentricity_earth_orbit(jc);
  double m = tcalc_geom_mean_anomaly_sun(jc);

  double y = tan(epsilon * M_PI / 180.0 / 2.0);
  double sin2l0 = sin(2.0 * l0 * M_PI / 180.0);
  double sinm = sin(m * M_PI / 180.0);
  double cos2l0 = cos(2.0 * l0 * M_PI / 180.0);
  double sin4l0 = sin(4.0 * l0 * M_PI / 180.0);
  double sin2m = sin(2.0 * m * M_PI / 180.0);
  double E_time;

  y *= y;

  E_time = y * sin2l0 - 2.0 * e * sinm + 4.0 * e * y * sinm * cos2l0
           - 0.5 * y * y * sin4l0 - 1.25 * e * e * sin2m;

  /* in minutes of time */
  return (4.0 * 180.0 / M_PI * E_time);
} /* tcalc_eqtime() */

static double
tcalc_solar_decl(double jc)

{
  double e = tcalc_obliquity_correction(jc);
  double lambda = tcalc_sun_apparent_long(jc);
  double sint = sin(e * M_PI / 180.0) * sin(lambda * M_PI / 180.0);
  double theta = asin(sint);

  /* in radians */
  return (theta);
} /* tcalc_solar_decl() */

static double
tcalc_julian_day(int year, int month, double day)

{
  double a, b;
  double jd;

  if (month <= 2)
    {
      year -= 1;
      month += 12;
    }

  a = floor(year / 100.0);
  b = 2.0 - a + floor(a / 4.0);

  jd = floor(365.25*(year + 4716)) + floor(30.6001*(month+1)) + day + b - 1524.5;

  return (jd);
} /* tcalc_julian_day() */

static double
tcalc_julian_cent(double jd)

{
  return ((jd - 2451545.0)/36525.0);
} /* tcalc_julian_cent() */

static int
tcalc_is_leap_year(int year)

{
  return ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0);
}

static void
tcalc_get_date(time_t t, int *year, int *month, double *day)

{
  struct tm *tm_p;

  tm_p = gmtime(&t);

  *year = tm_p->tm_year + 1900;
  *month = tm_p->tm_mon + 1;
  *day = tm_p->tm_mday;
}
