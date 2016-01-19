/* common.c
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
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <sys/time.h>

#include "common.h"

#include <gsl/gsl_math.h>

/*
 * time stamp of 01-Jan-2000 00:00:00 UTC
 * command used: date -d "1/1/2000 00:00:00 UTC" +"%s"
 */
static time_t t0_fday = 946684800;

/*
 * time stamp of 22-May-1960 12:00:00 UTC
 * command used: date -d "5/22/2000 12:00:00 UTC" +"%s"
 */
static time_t t0_lunar = -303307200;

static char tz_str[] = "TZ=GMT";

time_t
fday2timet(double fday)
{
  time_t t;

  /* find time in seconds since 1-Jan-2000 00:00:00 UTC */
  t = (time_t) (fday * 86400.0);

  /* add time to the epoch 1-Jan-1970 00:00:00 UTC */
  t += t0_fday;

  return (t);
} /* fday2timet() */

double
time2fday(time_t t)
{
  double fday;

  fday = (double) (t - t0_fday);
  fday /= 86400.0;

  return (fday);
} /* time2fday() */

/*
date2timet()
  sec \in [0,60]
  min \in [0,59]
  hour \in [0,23]
  day \in [1,31]
  month \in [1,12]
  year \in [*,*]
*/

time_t
date2timet(int sec, int min, int hour, int day, int month, int year)
{
  struct tm tmp;
  time_t t;

  tmp.tm_sec = sec;
  tmp.tm_min = min;
  tmp.tm_hour = hour;
  tmp.tm_mday = day;
  tmp.tm_mon = month - 1;
  tmp.tm_year = year - 1900;
  tmp.tm_isdst = 0;

  /* use GMT time */
  putenv(tz_str);

  t = mktime(&tmp);

  return t;
} /* date2timet() */

/*
time2lunar
  Convert a timestamp to lunar phase

Inputs: t - timestamp

Return: lunar phase in radians
*/

double
time2lunar(time_t t)
{
  static double phi0 = 1.93131118802;        /* 110.65598 deg */
  static double w0 = 2.0 * M_PI / 2380752.0; /* 2pi / 27.555 days */
  double phi;

  phi = w0 * (t - t0_lunar) + phi0;

#if 1
  while (phi > M_PI)
    phi -= M_PI;
  while (phi < 0.0)
    phi += M_PI;
#endif

  return phi;
} /* time2lunar() */

/*
get_season()
  Determine the day of year for a given fday

Inputs: t - timestamp in seconds since 1-Jan-1970 00:00:00 UTC

Return: floating point days since Jan 1, in the range 0 to 365
*/

double
get_season(time_t t)
{
  struct tm tm_p;
  double days;

  gmtime_r(&t, &tm_p);

  days = (double) tm_p.tm_yday;
  days += tm_p.tm_hour / 24.0 +
          tm_p.tm_min / 1440.0 +
          tm_p.tm_sec / 86400.0;

  return (days);
} /* get_season() */

/*
get_year()
  Determine the year for a given timestamp

Inputs: t - timestamp

Return: floating point year
*/

double
get_year(time_t t)
{
  struct tm tm_p;
  double year;
  double days;

  gmtime_r(&t, &tm_p);

  days = (double) tm_p.tm_yday;
  days += tm_p.tm_hour / 24.0 +
          tm_p.tm_min / 1440.0 +
          tm_p.tm_sec / 86400.0;

  year = (double) tm_p.tm_year + 1900.0 + days / 365.25;

  return (year);
} /* get_year() */

int
get_doy(time_t t)
{
  struct tm tm_p;

  gmtime_r(&t, &tm_p);

  return (tm_p.tm_yday + 1);
} /* get_doy() */

/*
get_localtime()
  Calculate the localtime for a given timestamp

Inputs: t         - timestamp
        longitude - longitude in radians of where to compute the
                    local time

Return: local time in floating point hours
*/

double
get_localtime(time_t t, double longitude)
{
  struct tm tm_p;
  double lt;

  /* convert utc to local time by adding longitude offset */
  lt = (double) t + longitude / (15.0 * M_PI / 180.0 / 3600.0);

  t = (time_t) lt;

  gmtime_r(&t, &tm_p);

  lt = (double) tm_p.tm_hour +
       (double) tm_p.tm_min / 60.0 +
       (double) tm_p.tm_sec / 3600.0;

  return (lt);
} /* get_localtime() */

double
get_ut(time_t t)
{
  struct tm tm_p;
  double ut;

  gmtime_r(&t, &tm_p);

  ut = (double) tm_p.tm_hour +
       (double) tm_p.tm_min / 60.0 +
       (double) tm_p.tm_sec / 3600.0;

  return (ut);
}

/*
get_fday()
  Convert parameterized time to fday

Inputs: year       - year
        longitude  - longitude in radians
        local_time - local time in hours
        season     - day of year (0 - 365)
*/

double
get_fday(int year, double longitude, double local_time, int season)
{
  struct tm tminfo;
  time_t offset;
  double fday;
  double utc;
  double lt = local_time * 3600.0;
  double s = season * 86400.0;

  /* use GMT time */
  putenv(tz_str);

  /*
   * get time offset since epoch of 00:00:00 UTC of the specified year
   */
  tminfo.tm_sec = 0;
  tminfo.tm_min = 0;
  tminfo.tm_hour = 0;
  tminfo.tm_mday = 1;
  tminfo.tm_mon = 0;
  tminfo.tm_year = year - 1900;
  tminfo.tm_isdst = 0;

  offset = mktime(&tminfo);

  /* get offset since 01-Jan-2000 00:00:00 UTC and convert to days */
  fday = (double) (offset - t0_fday);

  /* convert local time to utc by adding longitude offset */
  utc = lt - longitude / (15.0 * M_PI / 180.0 / 3600.0);

  /* add day of year and time */
  fday += s;
  fday += utc;

  fday /= 86400.0;
  
  return (fday);
} /* get_fday() */

/*
doy2md()
  Convert a day of year to month/day, accounting for leap years
*/

int
doy2md(int year, int doy, int *month, int *day)
{
  int monthtable[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  int mm = 0;
  int monthdays = 0;
  int newday;

  if (doy > 366)
    {
      fprintf(stderr, "doy2md: error: doy = %d\n", doy);
      return 1;
    }

  if ((year % 4 == 0) &&
      ((!(year % 100 == 0)) || (year % 400 == 0)))
    monthtable[1] = 29;
  else
    monthtable[1] = 28;

  while (doy > monthdays)
    monthdays += monthtable[mm++];

   newday = doy - monthdays + monthtable[mm-1];

   *month = mm;
   *day = newday;

   return 0;
} /* doy2md() */

int
is_leap_year(int year)
{
  if ((year % 4 == 0) &&
      ((!(year % 100 == 0)) || (year % 400 == 0)))
    return 1;

  return 0;
} /* is_leap_year() */

double
lon_fix(double phi)
{
  while (phi > M_PI)
    phi -= 2.0 * M_PI;

  while (phi < -M_PI)
    phi += 2.0 * M_PI;

  return phi;
} /* lon_fix() */

/* wrap x to the interval [0,360] */
double
wrap360(double x)
{
  double y = fmod(x, 360.0);

  if (y < 0.0)
    y += 360.0;

  return y;
} /* wrap360() */

/* wrap x to the interval [-180,180] */
double
wrap180(double x)
{
  double y = fmod(x, 360.0);

  if (y < -180.0)
    y += 360.0;
  if (y > 180.0)
    y -= 360.0;

  return y;
} /* wrap180() */

/* wrap x to the interval [0,2pi] */
double
wrap2pi(double x)
{
  double y = fmod(x, 2.0 * M_PI);

  if (y < 0.0)
    y += 2.0 * M_PI;

  return y;
} /* wrap2pi() */

/* wrap x to the interval [-pi,pi] */
double
wrappi(double x)
{
  double y = fmod(x, 2.0 * M_PI);

  if (y < -M_PI)
    y += 2.0 * M_PI;
  if (y > M_PI)
    y -= 2.0 * M_PI;

  return y;
} /* wrappi() */

/*
Given f(a) and f(b), find f(x) where x \in [a,b]
*/

double
linint(double x, double a, double b, double fa, double fb)
{
  double slope, intcp;

  slope = (fa - fb) / (a - b);
  intcp = fa - slope*a;

  return (slope*x + intcp);
}

double **
array2d_alloc(size_t size_x, size_t size_y)
{
  double **a;
  size_t i;

  a = malloc(size_x * sizeof(double *));
  if (!a)
    {
      fprintf(stderr, "array2d_alloc: malloc failed: %s\n",
              strerror(errno));
      return 0;
    }

  for (i = 0; i < size_x; ++i)
    {
      a[i] = malloc(size_y * sizeof(double));
      if (!a[i])
        {
          fprintf(stderr, "array2d_alloc: malloc failed: %s\n",
                  strerror(errno));
          return 0;
        }
    }

  return a;
} /* alloc2d() */

void
array2d_free(size_t size_x, double **a)
{
  size_t i;

  for (i = 0; i < size_x; ++i)
    free(a[i]);

  free(a);
}

/* C = A x B
idx 0: r component
idx 1: theta component
idx 2: phi component
*/
int
sphcross(double A[3], double B[3], double C[3])
{
  C[0] = (A[1]*B[2] - A[2]*B[1]);
  C[1] = (A[2]*B[0] - A[0]*B[2]);
  C[2] = (A[0]*B[1] - A[1]*B[0]);

  return 0;
} /* sphcross() */

/* compute Cartesian components of spherical basis vectors */
int
sph_basis(const double theta, const double phi,
          double rhat[3], double that[3], double phat[3])
{
  rhat[0] = sin(theta) * cos(phi);
  rhat[1] = sin(theta) * sin(phi);
  rhat[2] = cos(theta);

  that[0] = cos(theta) * cos(phi);
  that[1] = cos(theta) * sin(phi);
  that[2] = -sin(theta);

  phat[0] = -sin(phi);
  phat[1] = cos(phi);
  phat[2] = 0.0;

  return 0;
} /* sph_basis() */

double
vec_dot(const double a[3], const double b[3])
{
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

double
vec_norm(const double v[3])
{
  return gsl_hypot3(v[0], v[1], v[2]);
}

int
vec_unit(const double v[3], double unit[3])
{
  size_t i;
  double norm = vec_norm(v);

  for (i = 0; i < 3; ++i)
    unit[i] = v[i] / norm;

  return 0;
}

double
time_diff(struct timeval a, struct timeval b)
{
  double sec1, sec2;

  sec1 = (double) a.tv_sec + a.tv_usec * 1.0e-6;
  sec2 = (double) b.tv_sec + b.tv_usec * 1.0e-6;

  return sec2 - sec1;
}
