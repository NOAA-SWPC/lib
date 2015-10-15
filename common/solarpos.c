/* solarpos.c
 * 
 * Copyright (C) 2007 Patrick Alken
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

#include <gsl/gsl_math.h>

#include "solarpos.h"

static void solarpos_docalc(time_t t, double latitude, double longitude,
                            solarpos_workspace *w);

solarpos_workspace *
solarpos_alloc(void)
{
  solarpos_workspace *w;

  w = (solarpos_workspace *) calloc(1, sizeof(solarpos_workspace));
  if (!w)
    return 0;

  solarpos_init(w);

  return (w);
} /* solarpos_alloc() */

void
solarpos_free(solarpos_workspace *w)
{
  if (!w)
    return;

  free(w);
} /* solarpos_free() */

void
solarpos_init(solarpos_workspace *w)
{
  S_init(&(w->pd));

  w->pd.temp = 27.0;
  w->pd.press = 1006.0;
} /* solarpos_init() */

/*
solarpos_calc_zenith()
  Calculate the solar zenith angle

Inputs: t         - timestamp in UTC
        latitude  - latitude (in radians)
        longitude - longitude (in radians)
        result    - (output) where to store zenith angle (in radians)
        w         - solarpos workspace

Return: 0 on success, 1 on failure
*/

int
solarpos_calc_zenith(time_t t, double latitude, double longitude,
                     double *result, solarpos_workspace *w)
{
  w->pd.timezone = 0.0;
  solarpos_docalc(t, latitude, longitude, w);

  *result = ((double)w->pd.zenetr * M_PI / 180.0);

  return 0;
} /* solarpos_calc_zenith() */

/*
solarpos_calc_sunrs()
  Calculate sunrise/sunset local times

Inputs: t         - timestamp in UTC
        latitude  - latitude (in radians)
        longitude - longitude (in radians)
        sunrise   - (output) where to store local sunrise (minutes)
        sunset    - (output) where to store local sunset (minutes)
        w         - solarpos workspace

Notes: The output of this function has been tested against the
       tcalc_xxx routines and the results agree to within less
       than 4.5 minutes.
*/

int
solarpos_calc_sunrs(time_t t, double latitude, double longitude,
                    double *sunrise, double *sunset,
                    solarpos_workspace *w)
{
  double sr, ss;

  w->pd.timezone = longitude * 180.0 / M_PI / 15.0;
  solarpos_docalc(t, latitude, longitude, w);

  /*
   * the computed sunrise/sunset times are in UTC, so convert
   * them to the appropriate time zone:
   *
   * LT = UT + phi/(15 deg/hr) = UT + phi / (pi/12/60 rad/min)
   */
  /*sr = (double) w->pd.sretr + longitude * 12.0 / M_PI * 60.0;
  ss = (double) w->pd.ssetr + longitude * 12.0 / M_PI * 60.0;*/
  sr = (double) w->pd.sretr;
  ss = (double) w->pd.ssetr;

  *sunrise = sr;
  *sunset = ss;

  return 0;
} /* solarpos_calc_sunrs() */

/*************************************
 *       INTERNAL ROUTINES           *
 *************************************/

/*
solarpos_docalc()
  Call the solar position calculation routine to compute
all solar quantities

Notes: w->pd.timezone must be initialized before calling this
       function
*/

static void
solarpos_docalc(time_t t, double latitude, double longitude,
                solarpos_workspace *w)
{
  long ret;
  struct tm *tmptr;
  time_t offset;

  tmptr = gmtime(&t);

  w->pd.longitude = (float) (longitude * 180.0 / M_PI);
  w->pd.latitude = (float) (latitude * 180.0 / M_PI);
  w->pd.year = tmptr->tm_year + 1900;
  w->pd.daynum = tmptr->tm_yday + 1;
  w->pd.hour = tmptr->tm_hour;
  w->pd.minute = tmptr->tm_min;
  w->pd.second = tmptr->tm_sec;

  ret = S_solpos(&(w->pd));
  S_decode(ret, &(w->pd));
} /* solarpos_docalc() */
