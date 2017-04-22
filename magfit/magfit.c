/*
 * magfit.c
 *
 * High-level routines for fitting models to magnetic vector
 * data.
 *
 * Models supported:
 * 1. 1D SECS
 * 2. 2D SECS (DF only)
 * 3. PCA
 * 4. Internal potential field (gaussint)
 * 5. Simple ring current model (rc)
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
  Allocate magfit workspace

Inputs: T      - model type
        params - parameters

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

  params.lat_spacing1d = 0.5;
  params.lat_spacing2d = 2.0;
  params.lat_min = -60.0;
  params.lat_max = 60.0;
  params.lon_spacing = 3.0;
  params.lon_min = 0.0; /* these should be set by user */
  params.lon_max = 0.0;
  params.R = R_EARTH_KM + 110.0;
  params.lmax = MAGFIT_SECS_LMAX;
  params.secs_flags = MAGFIT_SECS_FLG_FIT_DF;

  params.pca_modes = 16;

  params.qdmax = 40.0;

  params.rc_p = 2;
  params.rc_fit_Y = 1;
  params.rc_subtract_crust = 0;

  return params;
}

int
magfit_reset(magfit_workspace *w)
{
  int status = (w->type->reset) (w->state);
  return status;
}

/*
magfit_add_datum()
  Add single datum to magfit workspace

Inputs: t     - timestamp (CDF_EPOCH)
        r     - radius (km)
        theta - colatitude (radians)
        phi   - colatitude (radians)
        qdlat - QD latitude (degrees)
        B     - magnetic field NEC vector (nT)
        w     - workspace
*/

int
magfit_add_datum(const double t, const double r, const double theta, const double phi,
                 const double qdlat, double B[3], magfit_workspace *w)
{
  int s;

  s = (w->type->add_datum)(t, r, theta, phi, qdlat, B, w->state);
  if (s)
    {
      fprintf(stderr, "magfit_add_datum: error adding data: %d\n", s);
      return GSL_FAILURE;
    }

  return GSL_SUCCESS;
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
  int s;
  size_t n = 0;
  size_t i;

  for (i = 0; i < tptr->n; ++i)
    {
      size_t didx = i + tptr->start_idx;
      double t = data->t[didx];
      double r = data->r[didx];
      double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
      double phi = data->longitude[didx] * M_PI / 180.0;
      double qdlat = data->qdlat[didx];
      double B[3];

      if (!SATDATA_AvailableData(data->flags[didx]))
        continue;

      /* fit only low-latitude data */
      if (fabs(qdlat) > w->params.qdmax)
        continue;

      B[0] = tptr->Bx[i];
      B[1] = tptr->By[i];
      B[2] = tptr->Bz[i];

      s = (w->type->add_datum)(t, r, theta, phi, qdlat, B, w->state);
      if (s)
        {
          fprintf(stderr, "magfit_add_track: error adding data\n");
          break;
        }

      ++n;
    }

  return n;
}

/*
magfit_fit()
  Fit magnetic model to previously added tracks

Inputs: rnorm - residual norm || y - A x ||
        snorm - solution norm || L x ||
        w     - workspace

Return: success/error

Notes:
1) Data must be added to workspace via magfit_add_track()
*/

int
magfit_fit(double * rnorm, double * snorm, magfit_workspace *w)
{
  int status = (w->type->fit)(rnorm, snorm, w->state);
  return status;
}

/*
magfit_apply_track()
  Subtract previously fitted model from a track

Inputs: tptr - track pointer
        data - satellite data
        w    - workspace

Return: number of data added
*/

int
magfit_apply_track(track_data *tptr, satdata_mag *data, magfit_workspace *w)
{
  int s = 0;
  size_t i;

  for (i = 0; i < tptr->n; ++i)
    {
      size_t didx = i + tptr->start_idx;
      double t = data->t[didx];
      double r = data->r[didx];
      double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
      double phi = data->longitude[didx] * M_PI / 180.0;
      double B[3];

      if (!SATDATA_AvailableData(data->flags[didx]))
        continue;

      magfit_eval_B(t, r, theta, phi, B, w);

      /* add RC correction model to data->B_ext */
      SATDATA_VEC_X(data->B_ext, didx) += B[0];
      SATDATA_VEC_Y(data->B_ext, didx) += B[1];
      SATDATA_VEC_Z(data->B_ext, didx) += B[2];
    }

  /* recalculate track residuals */
  track_calc_residuals(tptr, (const satdata_mag *) data);

  return s;
}

/*
magfit_eval_B()
  Evaluate magnetic field at a given (t,r,theta,phi) using
previously computed coefficients

Inputs: t     - timestamp (CDF_EPOCH)
        r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        B     - (output) magnetic field vector (nT)
        w     - workspace
*/

int
magfit_eval_B(const double t, const double r, const double theta, const double phi,
              double B[3], magfit_workspace *w)
{
  int status = (w->type->eval_B)(t, r, theta, phi, B, w->state);
  return status;
}

/*
magfit_eval_J()
  Evaluate current density at a given (r,theta) using
previously computed coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        J     - (output) current density vector [A/km]
        w     - workspace
*/

int
magfit_eval_J(const double r, const double theta, const double phi,
              double J[3], magfit_workspace *w)
{
  int status = (w->type->eval_J)(r, theta, phi, J, w->state);
  return status;
}

/*
magfit_print_map()
  Print lat/lon map of current flow and magnetic field components

Inputs: fp - output file
        r  - radius of shell for magnetic field maps (km)
        w  - workspace
*/

int
magfit_print_map(FILE *fp, const double r, magfit_workspace *w)
{
  int s = 0;
  const magfit_parameters *params = &(w->params);
  const size_t p = (w->type->ncoeff)(w->state);
  double (*eval_chi)(const double theta, const double phi, void * state) = w->type->eval_chi;
  void *chi_state = w->state;
  magfit_workspace *gaussint_p = NULL;
  double lat, lon;
  size_t i;

  if (eval_chi == NULL)
    {
      /*
       * For 1D and 2D SECS, generate magnetic field values on a spherical shell and
       * invert for internal gauss coefficients. Then we can use standard formulas to
       * determine stream function chi
       */
      double rnorm, snorm;

      gaussint_p = magfit_alloc(magfit_gaussint, &(w->params));

      fprintf(stderr, "magfit_print_map: building dataset for for gauss coefficient inversion...");

      for (lon = -180.0; lon <= 180.0; lon += 2.5)
        {
          double phi = lon * M_PI / 180.0;

          for (lat = -89.9; lat <= 89.9; lat += 1.0)
            {
              double theta = M_PI / 2.0 - lat * M_PI / 180.0;
              double B[3];

              /* compute magnetic field vector */
              s += (w->type->eval_B)(0.0, r, theta, phi, B, w->state);

              /* add B to gaussint workspace */
              (gaussint_p->type->add_datum)(0.0, r, theta, phi, 0.0, B, gaussint_p->state);
            }
        }

      fprintf(stderr, "done\n");

      /* fit internal gauss coefficients to magnetic field data */
      fprintf(stderr, "magfit_print_map: inverting for gauss coefficients...");
      magfit_fit(&rnorm, &snorm, gaussint_p);
      fprintf(stderr, "done\n");

      eval_chi = gaussint_p->type->eval_chi;
      chi_state = gaussint_p->state;
    }

  i = 1;
  fprintf(fp, "# Number of basis functions: %zu\n", p);
  fprintf(fp, "# Current shell radius: %g [km]\n", params->R);
  fprintf(fp, "# Magnetic field shell radius: %g [km]\n", r);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: chi (kA / nT)\n", i++);
  fprintf(fp, "# Field %zu: B_x (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_y (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_z (nT)\n", i++);

  for (lon = -180.0; lon <= 180.0; lon += 1.0)
    {
      double phi = lon * M_PI / 180.0;

      for (lat = -89.9; lat <= 89.9; lat += 1.0)
        {
          double theta = M_PI / 2.0 - lat * M_PI / 180.0;
          double B[3];
          double chi;

          /* compute current stream function */
          chi = (eval_chi)(theta, phi, chi_state);

          /* compute magnetic field vector */
          s += (w->type->eval_B)(0.0, r, theta, phi, B, w->state);

          fprintf(fp, "%f %f %f %f %f %f\n",
                  lon,
                  lat,
                  chi,
                  B[0],
                  B[1],
                  B[2]);
        }

      fprintf(fp, "\n");
    }

  fprintf(fp, "\n\n");

  fflush(fp);

  if (gaussint_p)
    magfit_free(gaussint_p);

  return s;
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
      fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
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
      fprintf(fp, "# Field %zu: model J_X (A/km)\n", i++);
      fprintf(fp, "# Field %zu: model J_Y (A/km)\n", i++);
      fprintf(fp, "# Field %zu: model J_Z (A/km)\n", i++);
      return 0;
    }

  for (i = 0; i < tptr->n; ++i)
    {
      size_t didx = i + tptr->start_idx;
      double t = data->t[didx];
      double r = data->r[didx];
      double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
      double phi = data->longitude[didx] * M_PI / 180.0;
      double B_model[3], J_model[3];

      if (!SATDATA_AvailableData(data->flags[didx]))
        continue;

      /* fit only low-latitude data */
      if (fabs(data->qdlat[didx]) > w->params.qdmax)
        continue;

      /* compute model prediction */
      magfit_eval_B(t, r, theta, phi, B_model, w);
      magfit_eval_J(r, theta, phi, J_model, w);

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

/*
magfit_print_rms()
  Print rms between a track and magnetic model to a file

Inputs: header - print header?
        fp     - file pointer
        lon0   - reference longitude (degrees)
        tptr   - track data
        data   - satellite data
        w      - workspace
*/

int
magfit_print_rms(const int header, FILE *fp, const double lon0, const track_data *tptr,
                 const satdata_mag *data, magfit_workspace *w)
{
  size_t i;
  double rms[3] = { 0.0, 0.0, 0.0 };
  size_t n = 0;

  if (header)
    {
      i = 1;
      fprintf(fp, "# Field %zu: timestamp of equator crossing\n", i++);
      fprintf(fp, "# Field %zu: dlon (degrees)\n", i++);
      fprintf(fp, "# Field %zu: X rms (nT)\n", i++);
      fprintf(fp, "# Field %zu: Y rms (nT)\n", i++);
      fprintf(fp, "# Field %zu: Z rms (nT)\n", i++);
      return 0;
    }

  for (i = 0; i < tptr->n; ++i)
    {
      size_t didx = i + tptr->start_idx;
      double t = data->t[didx];
      double r = data->r[didx];
      double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
      double phi = data->longitude[didx] * M_PI / 180.0;
      double B_model[3];

      if (!SATDATA_AvailableData(data->flags[didx]))
        continue;

      /* use only low-latitude data */
      if (fabs(data->qdlat[didx]) > w->params.qdmax)
        continue;

      /* compute model prediction */
      magfit_eval_B(t, r, theta, phi, B_model, w);

      rms[0] += pow(tptr->Bx[i] - B_model[0], 2.0);
      rms[1] += pow(tptr->By[i] - B_model[1], 2.0);
      rms[2] += pow(tptr->Bz[i] - B_model[2], 2.0);
      ++n;
    }

  if (n > 0)
    {
      for (i = 0; i < 3; ++i)
        rms[i] = sqrt(rms[i] / (double)n);
    }

  fprintf(fp, "%ld %.4f %.2f %.2f %.2f\n",
          satdata_epoch2timet(tptr->t_eq),
          wrap180(tptr->lon_eq - lon0),
          rms[0],
          rms[1],
          rms[2]);
  fflush(fp);

  return 0;
}
