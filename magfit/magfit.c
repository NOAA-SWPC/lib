/*
 * magfit.c
 *
 * High-level routines for fitting models to magnetic vector
 * data
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <satdata/satdata.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

#include "common.h"
#include "interp.h"
#include "magfit.h"
#include "track.h"

/*
magfit_alloc()
  Allocate secs 1d workspace

Inputs: flags        - SECS1D_FLG_xxx
        lmax         - maximum degree for Legendre functions in expansion
        R_iono       - radius of ionosphere (km)
        pole_spacing - along-orbit latitude spacing of SECS poles (degrees)

Return: pointer to workspace
*/

magfit_workspace *
magfit_alloc(const magfit_type * T, const magfit_parameters * params)
{
  magfit_workspace *w;

  w = calloc(1, sizeof(magfit_workspace));
  if (!w)
    return 0;

  w->state = (T->alloc)(params);
  if (w->state == NULL)
    {
      magfit_free(w);
      return 0;
    }

  w->type = T;
  w->params = *params;

  return w;
}

void
magfit_free(magfit_workspace *w)
{
  if (w->state)
    (w->type->free)(w->state);

  free(w);
}

magfit_parameters
magfit_default_parameters(void)
{
  magfit_parameters params;

  params.lat_spacing = 0.5;
  params.R = R_EARTH_KM + 110.0;
  params.lmax = MAGFIT_SECS_LMAX;
  params.secs_flags = MAGFIT_SECS_FLG_FIT_DF;

  return params;
}

int
magfit_reset(magfit_workspace *w)
{
  int status = (w->type->reset) (w->state);
  return status;
}

/*
magfit_add_track()
  Add satellite data from a single track to LS system

Inputs: tptr - track pointer
        data - satellite data
        w    - workspace

Return: number of data added
*/

size_t
magfit_add_track(const track_data *tptr, const satdata_mag *data,
                 magfit_workspace *w)
{
  size_t n = (w->type->add_track)(tptr, data, w->state);
  return n;
}

/*
magfit_fit()
  Fit magnetic model to previously added tracks

Inputs: w - workspace

Return: success/error

Notes:
1) Data must be added to workspace via magfit_add_track()
*/

int
magfit_fit(magfit_workspace *w)
{
  int status = (w->type->fit)(w->state);
  return status;
}

/*
magfit_eval_B()
  Evaluate magnetic field at a given (r,theta) using
previously computed coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        B     - (output) magnetic field vector (nT)
        w     - workspace
*/

int
magfit_eval_B(const double r, const double theta,
              double B[3], magfit_workspace *w)
{
  int status = (w->type->eval_B)(r, theta, B, w->state);
  return status;
}

/*
magfit_eval_J()
  Evaluate current density at a given (r,theta) using
previously computed coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        J     - (output) current density vector [A/km]
        w     - workspace
*/

int
magfit_eval_J(const double r, const double theta,
              double J[3], magfit_workspace *w)
{
  int status = (w->type->eval_J)(r, theta, J, w->state);
  return status;
}

/*
magfit_print_track()
  Print a track with magnetic model to a file
*/

int
magfit_print_track(const int header, FILE *fp, const track_data *tptr,
                   const satdata_mag *data, magfit_workspace *w)
{
  size_t i;

  if (header)
    {
      i = 1;
      fprintf(fp, "# Field %zu: timestamp\n", i++);
      fprintf(fp, "# Field %zu: radius (km)\n", i++);
      fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: X measurement (nT)\n", i++);
      fprintf(fp, "# Field %zu: Y measurement (nT)\n", i++);
      fprintf(fp, "# Field %zu: Z measurement (nT)\n", i++);
      fprintf(fp, "# Field %zu: model B_X (nT)\n", i++);
      fprintf(fp, "# Field %zu: model B_Y (nT)\n", i++);
      fprintf(fp, "# Field %zu: model B_Z (nT)\n", i++);
      fprintf(fp, "# Field %zu: model J_X\n", i++);
      fprintf(fp, "# Field %zu: model J_Y\n", i++);
      fprintf(fp, "# Field %zu: model J_Z\n", i++);
      return 0;
    }

  for (i = 0; i < tptr->n; ++i)
    {
      size_t didx = i + tptr->start_idx;
      double r = data->altitude[didx] + data->R;
      double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
      double B_model[3], J_model[3];

      if (!SATDATA_AvailableData(data->flags[didx]))
        continue;

      /* fit only low-latitude data */
      if (fabs(data->qdlat[didx]) > MAGFIT_QDMAX)
        continue;

      /* compute SECS model prediction */
      magfit_eval_B(r, theta, B_model, w);
      magfit_eval_J(r, theta, J_model, w);

      fprintf(fp, "%ld %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
              satdata_epoch2timet(data->t[didx]),
              r,
              data->longitude[didx],
              data->latitude[didx],
              data->qdlat[didx],
              tptr->Bx[i],
              tptr->By[i],
              tptr->Bz[i],
              B_model[0],
              B_model[1],
              B_model[2],
              J_model[0],
              J_model[1],
              J_model[2]);
    }

  fprintf(fp, "\n\n");

  return 0;
}
