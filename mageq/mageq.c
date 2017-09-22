/*
 * mageq.c
 *
 * Calculate latitude of magnetic equator at a given geocentric radius
 * and longitude
 *
 * This is done by finding the theta which minimizes:
 *
 * |B_r (r,theta,phi)|
 *
 * holding r and phi fixed.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_deriv.h>

#include <msynth/msynth.h>

#include "mageq.h"

static double mageq_func_Br(double x, void *params);
static double mageq_func_lat(double phi, void *params);

mageq_workspace *
mageq_alloc()
{
  mageq_workspace *w;

  w = calloc(1, sizeof(mageq_workspace));

  w->s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent); 
  w->msynth_workspace_p = msynth_igrf_read(MSYNTH_IGRF_FILE);

  return w;
} /* mageq_alloc() */

void
mageq_free(mageq_workspace *w)
{
  if (w->s)
    gsl_min_fminimizer_free(w->s);

  if (w->msynth_workspace_p)
    msynth_free(w->msynth_workspace_p);

  free(w);
} /* mageq_free() */

/*
mageq_calc()
  Calculate latitude of magnetic equator

Inputs: longitude - geographic longitude (radians)
        r         - geocentric radius (km)
        t         - decimal year
        w         - mageq workspace

Return: geocentric latitude in radians
*/

double
mageq_calc(double longitude, double r, double t, mageq_workspace *w)
{
  int status;
  int iter = 0, max_iter = 100;
  double m = M_PI / 2.0;
  double a = M_PI / 2.0 - 70.0 * M_PI / 180.0;
  double b = M_PI / 2.0 + 70.0 * M_PI / 180.0;
  gsl_function F;
  mageq_params params;

  params.r = r;
  params.phi = longitude;
  params.t = t;
  params.w = w;

  F.function = &mageq_func_Br;
  F.params = &params;

  gsl_min_fminimizer_set (w->s, &F, m, a, b);

  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (w->s);

      m = gsl_min_fminimizer_x_minimum (w->s);
      a = gsl_min_fminimizer_x_lower (w->s);
      b = gsl_min_fminimizer_x_upper (w->s);

      status = gsl_min_test_interval (a, b, 1.0e-5, 0.0);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  /* m is theta, convert to latitude */
  return (M_PI / 2.0 - m);
} /* mageq_calc() */

/*
mageq_angle()
  Calculate the angle between the dip equator and the horizontal
(eastward) direction at a given longitude

Inputs: longitude - geographic longitude in radians
        r         - geocentric radius in km
        t         - decimal year
        w         - mageq workspace

Return: angle in radians
*/

double
mageq_angle(double longitude, double r, double t, mageq_workspace *w)
{
  gsl_function F;
  mageq_params params;
  double result, abserr;

  params.r = r;
  params.phi = longitude;
  params.t = t;
  params.w = w;

  F.function = &mageq_func_lat;
  F.params = &params;

  gsl_deriv_central(&F, longitude, 1e-2, &result, &abserr);

  return (atan(result));
} /* mageq_angle() */

static double
mageq_func_Br(double theta, void *params)
{
  mageq_params *p = (mageq_params *) params;
  mageq_workspace *w = p->w;
  double B[4];

  msynth_eval(p->t, p->r, theta, p->phi, B, w->msynth_workspace_p);

  /* return |B_z| = |B_r| */
  return fabs(B[2]);
} /* mageq_func_Br() */

static double
mageq_func_lat(double phi, void *params)
{
  mageq_params *p = (mageq_params *) params;
  double lat;

  lat = mageq_calc(phi, p->r, p->t, p->w);

  return lat;
} /* mageq_func_lat() */
