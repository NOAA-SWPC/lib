/*
 * msis.c
 * Patrick Alken
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_math.h>

#include <indices/indices.h>

#include "msis.h"

static int msis_call(float latitude, float longitude, int year,
                     float hour, float height, float f107,
                     float f107a, msis_workspace *w);

msis_workspace *
msis_alloc(size_t nalt, const char *f107_datadir)
{
  msis_workspace *w;

  w = calloc(1, sizeof(msis_workspace));
  if (!w)
    {
      fprintf(stderr, "msis_alloc: calloc failed: %s\n", strerror(errno));
      return 0;
    }

  w->msis_results = malloc(nalt * sizeof(msis_result));
  if (!w->msis_results)
    {
      fprintf(stderr, "msis_alloc: malloc failed: %s\n", strerror(errno));
      msis_free(w);
      return 0;
    }

  w->f107_workspace_p = f107_alloc(f107_datadir);
  if (!w->f107_workspace_p)
    {
      msis_free(w);
      return 0;
    }

  w->kp_workspace_p = kp_alloc(KP_IDX_FILE);
  if (!w->kp_workspace_p)
    {
      msis_free(w);
      return 0;
    }

  w->nalt = nalt;

  w->f107_override = -1.0;
  w->f107a_override = -1.0;

  return w;
} /* msis_alloc() */

void
msis_free(msis_workspace *w)
{
  if (w->msis_results)
    free(w->msis_results);

  if (w->f107_workspace_p)
    f107_free(w->f107_workspace_p);

  if (w->kp_workspace_p)
    kp_free(w->kp_workspace_p);

  free(w);
} /* msis_free() */

void
msis_f107_override(double f107, double f107a, msis_workspace *w)
{
  w->f107_override = f107;
  w->f107a_override = f107a;
}

/*
msis_calc()
  Calculate MSIS parameters

Inputs: theta  - geographic colatitude (radians)
        phi    - geographic longitude (radians)
        t      - timestamp
        altmin - minimum altitude in km
        altstp - altitude step size in km
        nalt   - number of altitude steps
        w      - msis workspace

Return: success or error

Notes: MSIS parameters are stored in w->msis_results
*/

int
msis_calc(double theta, double phi, time_t t, double altmin,
          double altstp, size_t nalt, msis_workspace *w)
{
  int s = 0;
  size_t i;
  int year, month, day;
  float hour;
  float flat = 90.0 - theta * 180.0 / M_PI,
        flon = phi * 180.0 / M_PI;
  float f107, f107a;
  double f107d, f107ad;
  double ap;
  struct tm *tm_p;
  float *T = w->T;
  float *D = w->D;

  if (nalt > w->nalt)
    {
      fprintf(stderr, "msis_calc: specified nalt exceeds previous value (%zu,%zu)\n",
              nalt, w->nalt);
      return 1;
    }

  if (w->f107_override > 0.0)
    f107d = w->f107_override;
  else
    f107_get(t, &f107d, w->f107_workspace_p);

  if (w->f107a_override > 0.0)
    f107ad = w->f107a_override;
  else
    f107a_get(t, &f107ad, w->f107_workspace_p);

  f107 = (float) f107d;
  f107a = (float) f107ad;

  ap_get(t, &ap, w->kp_workspace_p);
  w->ap[0] = ap;

  putenv("TZ=GMT");
  tm_p = gmtime(&t);

  year = tm_p->tm_year + 1900;
  month = tm_p->tm_mon + 1;
  day = tm_p->tm_mday;
  hour = (float) tm_p->tm_hour +
         (float) tm_p->tm_min / 60.0 +
         (float) tm_p->tm_sec / 3600.0;

  w->day_of_year = tm_p->tm_yday + 1;

  for (i = 0; i < nalt; ++i)
    {
      float height = altmin + i * altstp;

      msis_call(flat,
                flon,
                year,
                hour,
                height,
                f107,
                f107a,
                w);

      /* save results into w->msis_results */

      /* convert to m-3 */
      w->msis_results[i].n_He = D[0] * 1.0e6;
      w->msis_results[i].n_O = D[1] * 1.0e6;
      w->msis_results[i].n_N2 = D[2] * 1.0e6;
      w->msis_results[i].n_O2 = D[3] * 1.0e6;
      w->msis_results[i].n_Ar = D[4] * 1.0e6;
      w->msis_results[i].n_H = D[6] * 1.0e6;
      w->msis_results[i].n_N = D[7] * 1.0e6;

      w->msis_results[i].n_n = w->msis_results[i].n_He +
                               w->msis_results[i].n_O +
                               w->msis_results[i].n_N2 +
                               w->msis_results[i].n_O2 +
                               w->msis_results[i].n_Ar +
                               w->msis_results[i].n_H +
                               w->msis_results[i].n_N;

      w->msis_results[i].T_n = T[1];
    }

  return s;
} /* msis_calc() */

/*
msis_get_result()
*/

msis_result *
msis_get_result(size_t idx, msis_workspace *w)
{
  return (&(w->msis_results[idx]));
} /* msis_get_result() */

/*
msis_call()

Inputs: latitude  - latitude in degrees
        longitude - longitude in degrees
        year      - year
        hour      - decimal hour UT
        height    - altitude in km
        f107      - F10.7
        f107a     - F10.7A
        w         - msis workspace

Notes: w->day_of_year must be set to the correct day of the year
*/

static int
msis_call(float latitude, float longitude, int year,
          float hour, float height, float f107, float f107a,
          msis_workspace *w)
{
  int s = 0;
  size_t i;
  int iyd;
  int yy;
  int day_of_year;
  float sec;      /* UT in sec */
  float stl;      /* local time in hours */
  float ap[7];
  int mass;      /* mass number */

  day_of_year = w->day_of_year;

  /* year ignored in current model */

  /* get last 2 digits of year */
  yy = year - (year/100) * 100;
  iyd = yy * 1000 + day_of_year;

  sec = hour * 3600.0;
  stl = hour + longitude / 15.0;

  for (i = 0; i < 7; ++i)
    ap[i] = (float) w->ap[i];

  /* all gases */
  mass = 48;

  gtd7_(&iyd,
        &sec,
        &height,
        &latitude,
        &longitude,
        &stl,
        &f107a,
        &f107,
        ap,
        &mass,
        w->D,
        w->T);

  return s;
} /* msis_call() */
