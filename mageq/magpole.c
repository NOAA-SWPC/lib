/*
 * magpole.c
 *
 * Calculate latitude and longitude of magnetic pole
 *
 * This is done by finding the theta and phi which minimize:
 *
 * B_H^2(t,r,theta,phi)
 *
 * holding r fixed, where:
 *
 * B_H = sqrt( B_theta^2 + B_phi^2 )
 *
 * The method used is Monte Carlo sampling
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <msynth/msynth.h>

#include "magpole.h"

static double func_f(const gsl_vector * x, void *params);

/*
magpole_alloc()
  Allocate magpole workspace

Inputs: N - number of Monte Carlo samples
*/

magpole_workspace *
magpole_alloc(const size_t N)
{
  magpole_workspace *w;

  gsl_rng_env_setup();

  w = calloc(1, sizeof(magpole_workspace));

  w->N = N;
  w->igrf_workspace_p = msynth_igrf_read(MSYNTH_IGRF_FILE);
  w->rng_p = gsl_rng_alloc(gsl_rng_default);

  return w;
}

void
magpole_free(magpole_workspace *w)
{
  if (w->igrf_workspace_p)
    msynth_free(w->igrf_workspace_p);

  if (w->rng_p)
    gsl_rng_free(w->rng_p);

  free(w);
}

/*
magpole_calc()
  Calculate latitude and longitude of magnetic pole

Inputs: longitude - geographic longitude (radians)
        r         - geocentric radius (km)
        t         - decimal year
        w         - magpole workspace

Return: geocentric latitude in radians
*/

int
magpole_calc(const double t, const double r, magpole_workspace *w)
{
  const size_t N = w->N;
  const double theta_max = M_PI / 6.0;
  double x_data[2];
  gsl_vector_view x = gsl_vector_view_array(x_data, 2);
  magpole_params params;
  double BH_min = 1.0e8;
  double theta_min = 0.0;
  double phi_min = 0.0;
  size_t i;

  params.r = r;
  params.t = t;
  params.w = w;

  for (i = 0; i < N; ++i)
    {
      double BH;

      x_data[0] = gsl_rng_uniform(w->rng_p) * theta_max;         /* theta \in [0, theta_max] */
      x_data[1] = 2.0 * M_PI * gsl_rng_uniform(w->rng_p) - M_PI; /* phi \in [-pi, pi] */

      BH = func_f(&x.vector, &params);
      if (BH < BH_min)
        {
          BH_min = BH;
          theta_min = x_data[0];
          phi_min = x_data[1];
        }
    }

  w->theta_pole = theta_min;
  w->phi_pole = phi_min;

  return GSL_SUCCESS;
}

static double
func_f(const gsl_vector * x, void *params)
{
  magpole_params *p = (magpole_params *) params;
  magpole_workspace *w = p->w;
  const double theta = gsl_vector_get(x, 0);
  const double phi = gsl_vector_get(x, 1);
  double B[4], BH;

  msynth_eval(p->t, p->r, theta, phi, B, w->igrf_workspace_p);

  BH = gsl_hypot(B[0], B[1]);

  return (BH * BH);
}
