/*
 * hwm.c
 *
 * This module provides a C interface to the HWM model (fortran)
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_math.h>

#include <indices/indices.h>

#include "hwm.h"

/*
hwm_alloc()
  Allocate hwm workspace
*/

hwm_workspace *
hwm_alloc(const char *f107_datafile)
{
  hwm_workspace *w;

  w = calloc(1, sizeof(hwm_workspace));
  if (!w)
    {
      fprintf(stderr, "hwm_alloc: malloc failed: %s\n", strerror(errno));
      return 0;
    }

  w->f107_workspace_p = f107_alloc(f107_datafile);
  w->kp_workspace_p = kp_alloc(KP_IDX_FILE);

  w->f107_override = -1.0;
  w->f107a_override = -1.0;

  w->scale = 1.0;

  return w;
} /* hwm_alloc() */

void
hwm_free(hwm_workspace *w)
{
  if (w->f107_workspace_p)
    f107_free(w->f107_workspace_p);

  if (w->kp_workspace_p)
    kp_free(w->kp_workspace_p);

  free(w);
} /* hwm_free() */

void
hwm_f107_override(double f107, double f107a, hwm_workspace *w)
{
  w->f107_override = f107;
  w->f107a_override = f107a;
} /* hwm_f107_override() */

/*
hwm_set_error_scale()
  Set scaling factor for winds for error analysis
*/

int
hwm_set_error_scale(double scale, hwm_workspace *w)
{
  w->scale = scale;
  return 0;
}

/*
hwm_call()
  Call the HWM model

Inputs: theta     - geographic colatitude (radians)
        longitude - geographic longitude (radians)
        t         - timestamp (UTC)
        hstart    - minimum altitude in km
        hstep     - step size in km
        nalt      - number of altitudes to compute
        merid     - (output) where to store meridional winds
        zonal     - (output) where to store zonal winds
        w         - workspace

Return: 0 if error occurs, number of altitudes computed if success

Notes: units on output winds are m/s
*/

size_t
hwm_call(double theta, double longitude, time_t t,
         double hstart, double hstep, size_t nalt,
         double *merid, double *zonal, hwm_workspace *w)
{
  float f107, f107a;
  float lat, lon;
  int iyd;
  float ut, lt;
  float utsec;
  float ap[2], wind[2];
  double tmp;
  size_t i;
  struct tm *tm_p;

  lat = 90.0 - theta * 180.0 / M_PI;
  lon = longitude * 180.0 / M_PI;

  tm_p = gmtime(&t);

  iyd = 90000 + tm_p->tm_yday + 1;
  ut = (float) tm_p->tm_hour +
       (float) tm_p->tm_min / 60.0 +
       (float) tm_p->tm_sec / 3600.0;
  lt = ut + lon / 15.0;

  if (lt > 24.0)
    lt -= 24.0;

  utsec = ut * 3600.0;

  if (w->f107_override > 0.0)
    f107 = (float) w->f107_override;
  else
    {
      f107_get(t, &tmp, w->f107_workspace_p);
      f107 = (float) tmp;
    }

  if (w->f107a_override > 0.0)
    f107a = (float) w->f107a_override;
  else
    {
      f107a_get(t, &tmp, w->f107_workspace_p);
      f107a = (float) tmp;
    }

  ap_get(t, &tmp, w->kp_workspace_p);
  ap[0] = (float) tmp;
  ap_get(t - 3*3600, &tmp, w->kp_workspace_p);
  ap[1] = (float) tmp;

  for (i = 0; i < nalt; ++i)
    {
      float height = (float) (hstart + i * hstep);

#if 0
      fprintf(stderr, "IYD = %d, UT = %f, ALT = %f, LAT = %f, LON = %f, STL = %f, AP = %f\n", 
              iyd, ut, height, lat, lon, lt, ap[1]);
#endif

      hwm14_(&iyd,
             &utsec,
             &height,
             &lat,
             &lon,
             &lt,
             &f107a,
             &f107,
             ap,
             wind);

      if (merid)
        merid[i] = wind[0] * w->scale;
      if (zonal)
        zonal[i] = wind[1] * w->scale;
    }

  return nalt;
} /* hwm_call() */
