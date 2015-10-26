/*
 * mag_log.c
 *
 * Contains routines which log various steps of the inversion
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "interp.h"
#include "mag.h"
#include "pde.h"

#include "pde_common.c"

/*
mag_log_profile()
  Log final EEF value to output file

Inputs: header - 1 = print header, 0 = don't
        ntrack - track number
        kp     - KP index
        dir    - satellite flight direction
        w      - workspace
*/

int
mag_log_profile(const int header, const size_t ntrack,
                const double kp, const int dir, const mag_workspace *w)
{
  int s = 0;
  size_t i;
  mag_track *track = (mag_track *) &(w->track);
  time_t t = satdata_epoch2timet(track->t_eq);
  double lt = get_localtime(t, track->phi_eq);
  double doy = get_season(t);
  double F2_peak;

  if (header)
    {
      /* print header information */
      i = 1;
      log_proc(w->log_profile, "# Field %zu: track number\n", i++);
      log_proc(w->log_profile, "# Field %zu: timestamp of equator crossing (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
      log_proc(w->log_profile, "# Field %zu: longitude of equator crossing (degrees)\n", i++);
      log_proc(w->log_profile, "# Field %zu: local time (hours)\n", i++);
      log_proc(w->log_profile, "# Field %zu: season (day of year in [0,365])\n", i++);
      log_proc(w->log_profile, "# Field %zu: peak F^(2) value at equator (nT)\n", i++);
      log_proc(w->log_profile, "# Field %zu: peak height-integrated current density (mA/m)\n", i++);
      log_proc(w->log_profile, "# Field %zu: kp\n", i++);
      log_proc(w->log_profile, "# Field %zu: satellite direction (+1 north, -1 south)\n", i++);
      return s;
    }

  F2_peak = interp_xy(track->qdlat, track->F2, track->n, 0.0);

  log_proc(w->log_profile, "%5zu %ld %9.4f %6.3f %7.3f %6.2f %9.6f %3.1f %2d\n",
           ntrack,
           t,
           track->phi_eq * 180.0 / M_PI,
           lt,
           doy,
           F2_peak,
           w->EEJ[w->ncurr / 2] * 1000.0,
           kp,
           dir);

  return s;
} /* mag_log_profile() */

/*
mag_log_F2()
  Log computed F^(2) residuals to a file

Inputs: header - 1 = print header, 0 = don't
        w      - workspace

Notes
1) track is stored in w->track
*/

int
mag_log_F2(const int header, const mag_workspace *w)
{
  int s = 0;
  const size_t downsample = 5; /* downsample to keep file size reasonable */
  size_t i;
  mag_track *track = (mag_track *) &(w->track);

  if (header)
    {
      /* print header information */
      i = 1;
      log_proc(w->log_F2, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
      log_proc(w->log_F2, "# Field %zu: local time (hours)\n", i++);
      log_proc(w->log_F2, "# Field %zu: season (day of year in [0,365])\n", i++);
      log_proc(w->log_F2, "# Field %zu: geocentric radius (km)\n", i++);
      log_proc(w->log_F2, "# Field %zu: longitude (degrees)\n", i++);
      log_proc(w->log_F2, "# Field %zu: geocentric latitude (degrees)\n", i++);
      log_proc(w->log_F2, "# Field %zu: QD latitude (degrees)\n", i++);
      log_proc(w->log_F2, "# Field %zu: F_sat (nT)\n", i++);
      log_proc(w->log_F2, "# Field %zu: F_internal (nT)\n", i++);
      log_proc(w->log_F2, "# Field %zu: dF_external (nT)\n", i++);
      log_proc(w->log_F2, "# Field %zu: F^(1) (nT)\n", i++);
      log_proc(w->log_F2, "# Field %zu: internal Sq model (nT)\n", i++);
      log_proc(w->log_F2, "# Field %zu: external Sq model (nT)\n", i++);
      log_proc(w->log_F2, "# Field %zu: F^(2) = F^(1) - Sq_model (nT)\n", i++);
      log_proc(w->log_F2, "# Field %zu: F^(2) fit from EEJ model (nT)\n", i++);
      return s;
    }

  for (i = 0; i < track->n; i += downsample)
    {
      time_t unix_time = satdata_epoch2timet(track->t[i]);
      double lt = get_localtime(unix_time, track->phi[i]);

      log_proc(w->log_F2, "%ld %6.3f %6.2f %9.4f %9.4f %8.4f %8.4f %8.2f %8.2f %6.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
               unix_time,
               lt,
               get_season(unix_time),
               track->r[i],
               track->phi[i] * 180.0 / M_PI,
               track->lat_deg[i],
               track->qdlat[i],
               track->F[i],
               track->F_int[i],
               track->dF_ext[i],
               track->F1[i],
               track->Sq_int[i],
               track->Sq_ext[i],
               track->F2[i],
               track->F2_fit[i]);
    }

  log_proc(w->log_F2, "\n\n");

  return s;
} /* mag_log_F2() */

/*
mag_log_Sq_Lcurve()
  Log Sq filter L-curves to a file

Inputs: header - 1 = print header, 0 = don't
        w      - workspace
*/

int
mag_log_Sq_Lcurve(const int header, const mag_workspace *w)
{
  int s = 0;
  const size_t downsample = 1;
  size_t i;
  mag_sqfilt_workspace *sqfilt_p = w->sqfilt_workspace_p;

  if (header)
    {
      /* print header information */
      i = 1;
      log_proc(w->log_Sq_Lcurve, "# Field %zu: regularization parameter\n", i++);
      log_proc(w->log_Sq_Lcurve, "# Field %zu: residual norm ||y - A x||\n", i++);
      log_proc(w->log_Sq_Lcurve, "# Field %zu: solution norm ||L x||\n", i++);
      return s;
    }

  for (i = 0; i < sqfilt_p->nreg; i += downsample)
    {
      log_proc(w->log_Sq_Lcurve, "%.6e %.12e %.12e\n",
               gsl_vector_get(sqfilt_p->reg_param, i),
               gsl_vector_get(sqfilt_p->rho, i),
               gsl_vector_get(sqfilt_p->eta, i));
    }

  log_proc(w->log_Sq_Lcurve, "\n\n");

  return s;
} /* mag_log_Sq_Lcurve() */

/*
mag_log_Sq_Lcorner()
  Log Sq filter L-corners to a file

Inputs: header - 1 = print header, 0 = don't
        w      - workspace
*/

int
mag_log_Sq_Lcorner(const int header, const mag_workspace *w)
{
  int s = 0;
  size_t i;
  mag_sqfilt_workspace *sqfilt_p = w->sqfilt_workspace_p;

  if (header)
    {
      /* print header information */
      i = 1;
      log_proc(w->log_Sq_Lcorner, "# Field %zu: L-corner regularization parameter\n", i++);
      log_proc(w->log_Sq_Lcorner, "# Field %zu: L-corner residual norm ||y - A x||\n", i++);
      log_proc(w->log_Sq_Lcorner, "# Field %zu: L-corner solution norm ||L x||\n", i++);
      return s;
    }

  log_proc(w->log_Sq_Lcorner, "%.6e %.12e %.12e\n",
           gsl_vector_get(sqfilt_p->reg_param, sqfilt_p->reg_idx),
           gsl_vector_get(sqfilt_p->rho, sqfilt_p->reg_idx),
           gsl_vector_get(sqfilt_p->eta, sqfilt_p->reg_idx));

  log_proc(w->log_Sq_Lcorner, "\n\n");

  return s;
} /* mag_log_Sq_Lcorner() */

/*
mag_log_LC()
  Log computed line currents to a file, 1 row per profile

Inputs: header - 1 = print header, 0 = don't
        w      - workspace
*/

int
mag_log_LC(const int header, const mag_workspace *w)
{
  int s = 0;
  size_t i;
  mag_track *track = (mag_track *) &(w->track);
  time_t t = satdata_epoch2timet(track->t_eq);

  if (header)
    {
      /* print header information */
      i = 1;
      log_proc(w->log_LC, "# Number of line currents: %zu\n", w->ncurr);
      log_proc(w->log_LC, "# Field %zu: timestamp of equator crossing (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
      log_proc(w->log_LC, "# Field %zu: longitude of equator crossing (degrees)\n", i++);
      log_proc(w->log_LC, "# Fields %zu-%zu: line current profile (mA/m)\n", i, i + w->ncurr - 1);
      ++i;
      return s;
    }

  log_proc(w->log_LC, "%zu %9.4f", t, track->phi_eq * 180.0 / M_PI);

  for (i = 0; i < w->ncurr; ++i)
    log_proc(w->log_LC, " %9.6f", w->EEJ[i] * 1000.0);

  log_proc(w->log_LC, "\n");

  return s;
} /* mag_log_LC() */

/*
mag_log_EEJ()
  Log computed F^(2) residuals to a file

Inputs: header - 1 = print header, 0 = don't
        w      - workspace

Notes
1) interpolated track is stored in 'w' arrays
*/

int
mag_log_EEJ(const int header, const mag_workspace *w)
{
  int s = 0;
  size_t i;
  mag_eej_workspace *eej_p = w->eej_workspace_p;

  if (header)
    {
      /* print header information */
      i = 1;
      log_proc(w->log_EEJ, "# Field %zu: QD latitude (degrees)\n", i++);
      log_proc(w->log_EEJ, "# Field %zu: EEJ height-integrated current (mA/m)\n", i++);
      return s;
    }

  for (i = 0; i < w->ncurr; ++i)
    {
      double qdlat = -eej_p->qdlat_max + i * eej_p->dqdlat;

      log_proc(w->log_EEJ, "%8.4f %9.6f\n",
               qdlat,
               w->EEJ[i] * 1000.0);
    }

  log_proc(w->log_EEJ, "\n\n");

  return s;
} /* mag_log_EEJ() */

/*
mag_log_EEJ_Lcurve()
  Log EEJ L-curves to a file

Inputs: header - 1 = print header, 0 = don't
        w      - workspace
*/

int
mag_log_EEJ_Lcurve(const int header, const mag_workspace *w)
{
  int s = 0;
  const size_t downsample = 1;
  size_t i;
  mag_eej_workspace *eej_p = w->eej_workspace_p;

  if (header)
    {
      /* print header information */
      i = 1;
      log_proc(w->log_EEJ_Lcurve, "# Field %zu: regularization parameter\n", i++);
      log_proc(w->log_EEJ_Lcurve, "# Field %zu: residual norm ||y - A x||\n", i++);
      log_proc(w->log_EEJ_Lcurve, "# Field %zu: solution norm ||L x||\n", i++);
      return s;
    }

  for (i = 0; i < eej_p->nreg; i += downsample)
    {
      log_proc(w->log_EEJ_Lcurve, "%.6e %.12e %.12e\n",
               gsl_vector_get(eej_p->reg_param, i),
               gsl_vector_get(eej_p->rho, i),
               gsl_vector_get(eej_p->eta, i));
    }

  log_proc(w->log_EEJ_Lcurve, "\n\n");

  return s;
} /* mag_log_EEJ_Lcurve() */

/*
mag_log_EEJ_Lcorner()
  Log EEJ L-corners to a file

Inputs: header - 1 = print header, 0 = don't
        w      - workspace
*/

int
mag_log_EEJ_Lcorner(const int header, const mag_workspace *w)
{
  int s = 0;
  size_t i;
  mag_eej_workspace *eej_p = w->eej_workspace_p;

  if (header)
    {
      /* print header information */
      i = 1;
      log_proc(w->log_EEJ_Lcorner, "# Field %zu: L-corner regularization parameter\n", i++);
      log_proc(w->log_EEJ_Lcorner, "# Field %zu: L-corner residual norm ||y - A x||\n", i++);
      log_proc(w->log_EEJ_Lcorner, "# Field %zu: L-corner solution norm ||L x||\n", i++);
      return s;
    }

  log_proc(w->log_EEJ_Lcorner, "%.6e %.12e %.12e\n",
           gsl_vector_get(eej_p->reg_param, eej_p->reg_idx),
           gsl_vector_get(eej_p->rho, eej_p->reg_idx),
           gsl_vector_get(eej_p->eta, eej_p->reg_idx));

  log_proc(w->log_EEJ_Lcorner, "\n\n");

  return s;
} /* mag_log_EEJ_Lcorner() */

/*
mag_log_PDE()
  Log computed PDE solutions residuals to a file

Inputs: header - 1 = print header, 0 = don't
        w      - workspace
*/

int
mag_log_PDE(const int header, const mag_workspace *w)
{
  int s = 0;
  size_t i;
  pde_workspace *pde_p = w->pde_workspace_p;

  if (header)
    {
      /* print header information */
      i = 1;
      log_proc(w->log_PDE, "# Field %zu: QD latitude (degrees)\n", i++);
      log_proc(w->log_PDE, "# Field %zu: height-integrated J(E, u = 0) (A/m)\n", i++);
      log_proc(w->log_PDE, "# Field %zu: height-integrated J(E = 0, u) (A/m)\n", i++);
      return s;
    }

  for (i = 0; i < pde_p->ntheta; ++i)
    {
      log_proc(w->log_PDE, "%f %.12e %.12e\n",
               90.0 - pde_theta(i, pde_p) * 180.0 / M_PI,
               gsl_vector_get(pde_p->J_lat_E, i),
               gsl_vector_get(pde_p->J_lat_u, i));
    }

  log_proc(w->log_PDE, "\n\n");

  return s;
} /* mag_log_PDE() */

/*
mag_log_model()
  Log modeled PDE solution to a file

Inputs: header - 1 = print header, 0 = don't
        w      - workspace
*/

int
mag_log_model(const int header, const mag_workspace *w)
{
  int s = 0;
  size_t i;
  inverteef_workspace *inveef_p = w->inverteef_workspace_p;
  mag_eej_workspace *eej_p = w->eej_workspace_p;

  if (header)
    {
      /* print header information */
      i = 1;
      log_proc(w->log_model, "# Field %zu: QD latitude (degrees)\n", i++);
      log_proc(w->log_model, "# Field %zu: modeled current density (mA/m)\n", i++);
      log_proc(w->log_model, "# Field %zu: satellite current density (mA/m)\n", i++);
      return s;
    }

  for (i = 0; i < w->ncurr; ++i)
    {
      double qdlat = -eej_p->qdlat_max + i * eej_p->dqdlat;

      log_proc(w->log_model, "%f %.12e %.12e\n",
               qdlat,
               gsl_vector_get(inveef_p->J_pde, i) * 1000.0,
               w->EEJ[i] * 1000.0);
    }

  log_proc(w->log_model, "\n\n");

  return s;
} /* mag_log_model() */

/*
mag_log_EEF()
  Log final EEF value to output file

Inputs: header - 1 = print header, 0 = don't
        t      - timestamp of equator crossing
        phi    - longitude of equator crossing (radians)
        kp     - kp index
        w      - workspace
*/

int
mag_log_EEF(const int header, const time_t t, const double phi,
            const double kp, const mag_workspace *w)
{
  int s = 0;
  size_t i;
  double lt = get_localtime(t, phi);
  double doy = get_season(t);

  if (header)
    {
      /* print header information */
      i = 1;
      log_proc(w->log_EEF, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
      log_proc(w->log_EEF, "# Field %zu: decimal year (UT)\n", i++);
      log_proc(w->log_EEF, "# Field %zu: longitude (degrees)\n", i++);
      log_proc(w->log_EEF, "# Field %zu: local time (hours)\n", i++);
      log_proc(w->log_EEF, "# Field %zu: season (day of year in [0,365])\n", i++);
      log_proc(w->log_EEF, "# Field %zu: eastward electric field (mV/m)\n", i++);
      log_proc(w->log_EEF, "# Field %zu: relative error between modeled and observed current profile\n", i++);
      log_proc(w->log_EEF, "# Field %zu: peak height-integrated current density (mA/m)\n", i++);
      log_proc(w->log_EEF, "# Field %zu: kp\n", i++);
      return s;
    }

  log_proc(w->log_EEF, "%ld %10.4f %10.4f %8.4f %8.2f %10.6f %10.6f %10.6f %5.1f\n",
           t,
           get_year(t),
           phi * 180.0 / M_PI,
           lt,
           doy,
           w->EEF,
           w->RelErr,
           w->EEJ[w->ncurr / 2] * 1000.0,
           kp);

  return s;
} /* mag_log_EEF() */
