/*
 * mag.c
 *
 * This module contains routines for inverting a magnetic profile
 * to an EEJ current profile
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <cdf.h>

#include <gsl/gsl_math.h>

#include <satdata/satdata.h>

#include "common.h"
#include "track.h"

#include "green.h"
#include "inverteef.h"
#include "log.h"
#include "mag.h"
#include "pde.h"

static int mag_track_datagap(const size_t sidx, const size_t eidx,
                             const satdata_mag *data);
static int mag_compute_F1(const double t_eq, const double phi_eq,
                          const size_t sidx, const size_t eidx,
                          const satdata_mag *data, mag_workspace *w);
static int mag_output(const int header, const mag_workspace *w);
static int mag_callback_season(const double doy, const void *params);

mag_workspace *
mag_alloc(mag_params *params)
{
  mag_workspace *w;
  
  w = calloc(1, sizeof(mag_workspace));
  if (!w)
    return 0;

  w->params = params;

  if (params->output_file)
    {
      w->fp_output = fopen(params->output_file, "a");
      if (!w->fp_output)
        {
          fprintf(stderr, "mag_alloc: unable to open %s: %s\n", 
                  params->output_file, strerror(errno));
          mag_free(w);
          return 0;
        }
    }

  w->ncurr = params->ncurr;
  w->EEJ = malloc(w->ncurr * sizeof(double));

  w->green_workspace_p = green_alloc(params->sq_nmax_int);

  if (params->use_vector)
    {
      w->sqfilt_vector_workspace_p = mag_sqfilt_vector_alloc(params->sq_nmax_int, params->sq_mmax_int,
                                                             params->sq_nmax_ext, params->sq_mmax_ext);
    }
  else
    {
      w->sqfilt_scalar_workspace_p = mag_sqfilt_scalar_alloc(params->sq_nmax_int, params->sq_mmax_int,
                                                             params->sq_nmax_ext, params->sq_mmax_ext);
    }

  w->eej_workspace_p = mag_eej_alloc(params->year, params->ncurr,
                                     params->curr_altitude, params->qdlat_max);

  {
    pde_parameters pde_params;
    inverteef_parameters inverteef_params;

    /*
     * I found ntheta = 200 causes problems in the sparse matrix
     * solver in the pde module, reducing to 100 improves things
     * but there are still occasional failures with the lis iterative
     * solver
     */

    pde_params.rmin = (R_EARTH_KM + 65.0) * 1.0e3;
    pde_params.rmax = (R_EARTH_KM + 500.0) * 1.0e3;
    pde_params.nr = 200;
    pde_params.theta_min = 65.0 * M_PI / 180.0;
    pde_params.theta_max = 115.0 * M_PI / 180.0;
    pde_params.ntheta = 100;

    w->pde_workspace_p = pde_alloc(&pde_params);

    inverteef_params.ntheta = pde_params.ntheta;
    inverteef_params.theta_min = pde_params.theta_min;
    inverteef_params.theta_max = pde_params.theta_max;
    inverteef_params.ncurr = params->ncurr;

    w->inverteef_workspace_p = inverteef_alloc(&inverteef_params);
  }

  w->kp_workspace_p = kp_alloc(KP_IDX_FILE);

  w->log_general = log_alloc(LOG_APPEND|LOG_TIMESTAMP, "%s/invert.log", params->log_dir);
  w->log_profile = log_alloc(LOG_WRITE, "%s/profile.dat", params->log_dir);
  w->log_F2 = log_alloc(LOG_WRITE, "%s/F2.dat", params->log_dir);
  w->log_B2 = log_alloc(LOG_WRITE, "%s/B2.dat", params->log_dir);
  w->log_Sq_Lcurve = log_alloc(LOG_WRITE, "%s/Sq_lcurve.dat", params->log_dir);
  w->log_Sq_Lcorner = log_alloc(LOG_WRITE, "%s/Sq_lcorner.dat", params->log_dir);
  w->log_EEJ_Lcurve = log_alloc(LOG_WRITE, "%s/EEJ_lcurve.dat", params->log_dir);
  w->log_EEJ_Lcorner = log_alloc(LOG_WRITE, "%s/EEJ_lcorner.dat", params->log_dir);
  w->log_EEJ = log_alloc(LOG_WRITE, "%s/EEJ.dat", params->log_dir);
  w->log_LC = log_alloc(LOG_WRITE, "%s/LC.dat", params->log_dir);
  w->log_PDE = log_alloc(LOG_WRITE, "%s/PDE.dat", params->log_dir);
  w->log_model = log_alloc(LOG_WRITE, "%s/model.dat", params->log_dir);
  w->log_EEF = log_alloc(LOG_WRITE, "%s/EEF.dat", params->log_dir);

  /* initialize headers in log files */
  mag_log_profile(1, 0, 0.0, 1, w);
  mag_log_F2(1, w);
  mag_log_B2(1, w);
  mag_log_Sq_Lcurve(1, w);
  mag_log_Sq_Lcorner(1, w);
  mag_log_LC(1, w);
  mag_log_EEJ(1, w);
  mag_log_EEJ_Lcurve(1, w);
  mag_log_EEJ_Lcorner(1, w);
  mag_log_PDE(1, w);
  mag_log_model(1, w);
  mag_log_EEF(1, 0, 0.0, 0.0, w);

  mag_output(1, w);

  return w;
} /* mag_alloc() */

void
mag_free(mag_workspace *w)
{
  if (w->EEJ)
    free(w->EEJ);

  if (w->green_workspace_p)
    green_free(w->green_workspace_p);

  if (w->sqfilt_vector_workspace_p)
    mag_sqfilt_vector_free(w->sqfilt_vector_workspace_p);

  if (w->sqfilt_scalar_workspace_p)
    mag_sqfilt_scalar_free(w->sqfilt_scalar_workspace_p);

  if (w->eej_workspace_p)
    mag_eej_free(w->eej_workspace_p);

  if (w->pde_workspace_p)
    pde_free(w->pde_workspace_p);

  if (w->inverteef_workspace_p)
    inverteef_free(w->inverteef_workspace_p);

  if (w->kp_workspace_p)
    kp_free(w->kp_workspace_p);

  if (w->log_general)
    log_free(w->log_general);

  if (w->log_profile)
    log_free(w->log_profile);

  if (w->log_F2)
    log_free(w->log_F2);

  if (w->log_B2)
    log_free(w->log_B2);

  if (w->log_Sq_Lcurve)
    log_free(w->log_Sq_Lcurve);

  if (w->log_Sq_Lcorner)
    log_free(w->log_Sq_Lcorner);

  if (w->log_LC)
    log_free(w->log_LC);

  if (w->log_EEJ)
    log_free(w->log_EEJ);

  if (w->log_EEJ_Lcurve)
    log_free(w->log_EEJ_Lcurve);

  if (w->log_EEJ_Lcorner)
    log_free(w->log_EEJ_Lcorner);

  if (w->log_PDE)
    log_free(w->log_PDE);

  if (w->log_model)
    log_free(w->log_model);

  if (w->log_EEF)
    log_free(w->log_EEF);

  if (w->fp_output)
    fclose(w->fp_output);

  free(w);
} /* mag_free() */

/*
mag_preproc()
  Preprocess satellite data according to given parameters

Inputs: data - satellite data
        w    - workspace
*/

int
mag_preproc(const mag_params *params, track_workspace *track_p,
            satdata_mag *data, mag_workspace *w)
{
  int status = 0;

  /* flag tracks with high rms */
  {
    /*
     * rms thresholds; set threshold high to process data
     * during strong storms: 17 March 2015 storm has scalar rms of up to 120 nT.
     *
     * The 22-23 June 2015 storm has scalar rms up to 140 nT
     */
    /*const double thresh[] = { -1.0, -1.0, -1.0, 150.0 };*/
    const double thresh[] = { 210.0, 170.0, 150.0, 160.0 };
    size_t nrms = track_flag_rms("rms.dat", thresh, NULL, data, track_p);

    log_proc(w->log_general, "mag_preproc: flagged %zu/%zu (%.1f%%) tracks due to high rms\n",
            nrms, track_p->n, (double) nrms / (double) track_p->n * 100.0);
  }

  /* flag data outside [lt_min,lt_max] */
  {
    size_t nlt = track_flag_lt(params->lt_min, params->lt_max, NULL, data, track_p);

    log_proc(w->log_general, "mag_preproc: flagged %zu/%zu (%.1f%%) tracks due to LT window [%g,%g]\n",
            nlt, track_p->n, (double)nlt / (double)track_p->n * 100.0,
            params->lt_min, params->lt_max);
  }

  /* flag data outside [lon_min,lon_max] */
  {
    size_t nlon = track_flag_lon(params->lon_min, params->lon_max, NULL, data, track_p);

    log_proc(w->log_general, "mag_preproc: flagged %zu/%zu (%.1f%%) tracks due to longitude window [%g,%g]\n",
            nlon, track_p->n, (double)nlon / (double)track_p->n * 100.0,
            params->lon_min, params->lon_max);
  }

  /* flag data according to season */
  {
    size_t nseas = track_flag_season(mag_callback_season, params, data, track_p);

    log_proc(w->log_general, "mag_preproc: flagged %zu/%zu (%.1f%%) data due to seasonal windows [%g,%g],[%g,%g]\n",
            nseas, data->n, (double)nseas / (double)data->n * 100.0,
            params->season_min, params->season_max,
            params->season_min2, params->season_max2);
  }

  /* flag high kp data */
  {
    const double kp_min = 0.0;
    const double kp_max = MAG_MAX_KP;
    size_t nkp = track_flag_kp(kp_min, kp_max, data, track_p);

    log_proc(w->log_general, "mag_preproc: flagged %zu/%zu (%.1f%%) data due to kp [%g,%g]\n",
            nkp, data->n, (double)nkp / (double)data->n * 100.0,
            kp_min, kp_max);
  }

  /* last check: flag tracks with very few good data points left */
  {
    size_t nflag = track_flag_n(2000, data, track_p);

    log_proc(w->log_general, "mag_preproc: flagged data due to low data points: %zu/%zu (%.1f%%) data flagged)\n",
            nflag, data->n, (double)nflag / (double)data->n * 100.0);
  }

  {
    size_t nflagged = satdata_nflagged(data);
    size_t nleft = data->n - nflagged;
    size_t nflagged_track = track_nflagged(track_p);
    size_t nleft_track = track_p->n - nflagged_track;

    log_proc(w->log_general, "mag_preproc: total flagged data: %zu/%zu (%.1f%%)\n",
             nflagged, data->n, (double)nflagged / (double)data->n * 100.0);
    log_proc(w->log_general, "mag_preproc: total remaining data: %zu/%zu (%.1f%%)\n",
             nleft, data->n, (double)nleft / (double)data->n * 100.0);

    log_proc(w->log_general, "mag_preproc: total flagged tracks: %zu/%zu (%.1f%%)\n",
             nflagged_track, track_p->n, (double)nflagged_track / (double)track_p->n * 100.0);
    log_proc(w->log_general, "mag_preproc: total remaining tracks: %zu/%zu (%.1f%%)\n",
             nleft_track, track_p->n, (double)nleft_track / (double)track_p->n * 100.0);
  }

  return status;
} /* mag_preproc() */

/*
mag_proc()
  Main processing loop for satellite data. Read in
tracks which cross the geographic equator and send them
to other modules for processing.

Inputs: data - satellite data
        w    - workspace
*/

int
mag_proc(const mag_params *params, track_workspace *track_p,
         satdata_mag *data, mag_workspace *w)
{
  int status = 0;
  size_t i = 0;
  size_t ntrack = 0;     /* number of tracks processed */
  size_t nrejgap = 0,    /* number of rejections due to data gap */
         nrejlat = 0,    /* number of rejections due to latitude */
         nrejflag = 0;   /* number of rejections due to various flags */

  if (data->n < 1000)
    {
      log_proc(w->log_general, "mag_proc: too few data for processing (%zu)\n",
               data->n);
      return -1;
    }

  log_proc(w->log_general, "mag_proc: processing %zu satellite data (%zu tracks)\n",
           data->n, track_p->n);

  /* preprocess tracks using given parameters */
  mag_preproc(params, track_p, data, w);

  for (i = 0; i < track_p->n; ++i)
    {
      int s;
      track_data *tptr = &(track_p->tracks[i]);
      size_t sidx, eidx;
      double lon_eq, t_eq, lt_eq, kp;
      time_t unix_time;
      int dir;
      char buf[2048];

      /* discard flagged tracks */
      if (tptr->flags)
        {
          ++nrejflag;
          continue;
        }

      sidx = tptr->start_idx;
      eidx = tptr->end_idx;
      lon_eq = tptr->lon_eq;
      t_eq = tptr->t_eq;
      lt_eq = tptr->lt_eq;
      unix_time = satdata_epoch2timet(tptr->t_eq);

      kp_get(unix_time, &kp, w->kp_workspace_p);

      /* reject tracks with insufficient latitude coverage */
      if (fabs(data->qdlat[sidx]) < 35.0 || fabs(data->qdlat[eidx]) < 35.0)
        {
          ++nrejlat;
          continue;
        }

      /* check for data gaps */
      s = mag_track_datagap(sidx, eidx, data);
      if (s)
        {
          ++nrejgap;
          continue;
        }

      ++ntrack;
      dir = tptr->satdir;

      sprintf(buf, "%s", ctime(&unix_time));
      buf[strlen(buf) - 1] = '\0';

      fprintf(stderr, "mag_proc: found track %zu, %s, lon = %g, lt = %g, kp = %g, dir = %d\n",
              ntrack, buf, lon_eq, lt_eq, kp, dir);

      /*
       * store track in workspace, compute magnetic coordinates, and
       * F^(1) residuals
       */ 
      s = mag_compute_F1(t_eq, lon_eq * M_PI / 180.0, sidx, eidx, data, w);
      if (s)
        return s; /* error occurred */

      /* filter out Sq and compute B^(2) or F^(2) residuals */
      if (params->use_vector)
        s = mag_sqfilt_vector(w, w->sqfilt_vector_workspace_p);
      else
        s = mag_sqfilt_scalar(w, w->sqfilt_scalar_workspace_p);

      if (s)
        return s;

      /* invert for EEJ height-integrated current density */
      s = mag_eej_proc(&(w->track), w->EEJ, w->eej_workspace_p);
      if (s)
        return s;

      /* log information for this profile (track number, t_eq, lon_eq) */
      mag_log_profile(0, ntrack, kp, dir, w);

      /* print F^(2) or B^(2) residuals to log file */
      if (params->use_vector)
        mag_log_B2(0, w);
      else
        mag_log_F2(0, w);

      /* output Sq L-curve and corners */
      mag_log_Sq_Lcurve(0, w);
      mag_log_Sq_Lcorner(0, w);

      /* print line current profiles to log file */
      mag_log_LC(0, w);

      /* print EEJ current density to log file */
      mag_log_EEJ(0, w);

      /* output EEJ L-curve and corners */
      mag_log_EEJ_Lcurve(0, w);
      mag_log_EEJ_Lcorner(0, w);

      /* output EEJ current density to output directory */
      mag_output(0, w);

      /* stop processing here if we only want profiles */
      if (params->profiles_only)
        continue;

      /* solve PDE */
      s = pde_proc(unix_time, lon_eq * M_PI / 180.0, w->pde_workspace_p);

      /*
       * print PDE solution J(E,u=0) and J(E=0,u) to log file - solution
       * is printed even if pde solver fails to keep indexing
       * consistent in log files
       */
      mag_log_PDE(0, w);

      if (s)
        {
          log_proc(w->log_general, "mag_proc: pde_proc failed on profile %zu\n",
                   ntrack);
        }

      /* invert for EEF */
      {
        gsl_vector_view J_sat = gsl_vector_view_array(w->EEJ, w->ncurr);

        s += inverteef_calc(&J_sat.vector,
                            w->pde_workspace_p->J_lat_E,
                            w->pde_workspace_p->J_lat_u,
                            w->inverteef_workspace_p);

        /* compute final EEF value in mV/m */
        w->EEF = w->inverteef_workspace_p->E_scale * EEF_PHI_0 * 1.0e3;

        /* store relative error between modeled and satellite profiles */
        w->RelErr = w->inverteef_workspace_p->RelErr;

        log_proc(w->log_general,
                 "mag_proc: profile %zu: [t,LT,lon,EEF,RelErr,kp,dir] = [%d,%.2f,%7.2f,%6.3f,%.2f,%.1f,%d]\n",
                 ntrack,
                 unix_time,
                 lt_eq,
                 lon_eq,
                 w->EEF,
                 w->RelErr,
                 kp,
                 dir);
        fprintf(stderr, "mag_proc: electric field = %f mV/m\n", w->EEF);
        fprintf(stderr, "mag_proc: relative error = %f\n", w->RelErr);

        /* log modeled profile */
        mag_log_model(0, w);

        /* if no errors, log EEF value to output file; since EEF is only printed
         * upon success, this file could be desynced from the other log files
         * which are all synced together for each profile, regardless of failure
         */
        if (s == GSL_SUCCESS)
          mag_log_EEF(0, unix_time, lon_eq * M_PI / 180.0, kp, w);
      }
    }

  log_proc(w->log_general,
           "mag_proc: %zu/%zu tracks rejected due to insufficient latitude coverage\n",
           nrejlat, track_p->n);
  log_proc(w->log_general,
           "mag_proc: %zu/%zu tracks rejected due to latitude gap\n",
           nrejgap, track_p->n);
  log_proc(w->log_general,
           "mag_proc: %zu/%zu tracks rejected due to flags (LT, kp, rms)\n",
           nrejflag, track_p->n);

  return status;
} /* mag_proc() */

/*
mag_track_datagap()
  Check satellite track for data gaps
*/

static int
mag_track_datagap(const size_t sidx, const size_t eidx,
                  const satdata_mag *data)
{
  size_t i;

  for (i = sidx + 1; i <= eidx; ++i)
    {
      double prevlat = data->latitude[i - 1];
      double lat = data->latitude[i];

      if (fabs(lat - prevlat) > MAG_DLAT)
        return -1;
    }

  return 0;
} /* mag_track_datagap() */

/*
mag_compute_F1()
  When a track is found for processing, compute magnetic
latitudes and copy the track into mag workspace for the
desired latitudes. Then call a main field model for the internal
and external field values along the track and compute F^(1)
residuals

Inputs: t_eq   - time of equator crossing (CDF_EPOCH)
        phi_eq - longitude of equator crossing (radians)
        sidx   - start index of track in 'data'
        eidx   - end index of track in 'data'
        data   - satellite data
        w      - workspace

Notes:
1) Upon exit, w->track contains all the information for this
track including the main field model values, QD latitudes and F^(1)
residuals
*/

static int
mag_compute_F1(const double t_eq, const double phi_eq, const size_t sidx,
               const size_t eidx, const satdata_mag *data, mag_workspace *w)
{
  int s = 0;
  size_t i, j;
  size_t idx = 0;
  mag_track *track = &(w->track);

  for (i = sidx; i <= eidx; ++i)
    {
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double qdlat = data->qdlat[i];
      double B_int[4], B_crust[4], B_ext[4], B_tot[4];

      if (fabs(qdlat) > MAG_MAX_QD_LATITUDE)
        continue;

      B_int[0] = SATDATA_VEC_X(data->B_main, i);
      B_int[1] = SATDATA_VEC_Y(data->B_main, i);
      B_int[2] = SATDATA_VEC_Z(data->B_main, i);

      B_crust[0] = SATDATA_VEC_X(data->B_crust, i);
      B_crust[1] = SATDATA_VEC_Y(data->B_crust, i);
      B_crust[2] = SATDATA_VEC_Z(data->B_crust, i);

      B_ext[0] = SATDATA_VEC_X(data->B_ext, i);
      B_ext[1] = SATDATA_VEC_Y(data->B_ext, i);
      B_ext[2] = SATDATA_VEC_Z(data->B_ext, i);

      /* add crustal field to main field */
      for (j = 0; j < 3; ++j)
        {
          B_tot[j] = B_int[j] + B_crust[j] + B_ext[j];
          B_int[j] += B_crust[j];
        }

      B_int[3] = gsl_hypot3(B_int[0], B_int[1], B_int[2]);
      B_tot[3] = gsl_hypot3(B_tot[0], B_tot[1], B_tot[2]);

      /* store this track */
      track->t_eq = t_eq;
      track->phi_eq = phi_eq;
      track->t[idx] = data->t[i];
      track->theta[idx] = theta;
      track->lat_deg[idx] = 90.0 - theta*180.0/M_PI;
      track->phi[idx] = phi;
      track->r[idx] = data->R + data->altitude[i];
      track->thetaq[idx] = M_PI / 2.0 - qdlat * M_PI / 180.0;
      track->qdlat[idx] = qdlat;
      track->F[idx] = data->F[i];
      track->X[idx] = SATDATA_VEC_X(data->B, i);
      track->Y[idx] = SATDATA_VEC_Y(data->B, i);
      track->Z[idx] = SATDATA_VEC_Z(data->B, i);
      track->Bx_int[idx] = B_int[0];
      track->By_int[idx] = B_int[1];
      track->Bz_int[idx] = B_int[2];
      track->F_int[idx] = B_int[3];
      track->dF_ext[idx] = vec_dot(B_int, B_ext) / B_int[3];

      /* compute F^(1) residuals */
      track->F1[idx] = data->F[i] - B_tot[3];
      track->X1[idx] = SATDATA_VEC_X(data->B, i) - B_tot[0];
      track->Y1[idx] = SATDATA_VEC_Y(data->B, i) - B_tot[1];
      track->Z1[idx] = SATDATA_VEC_Z(data->B, i) - B_tot[2];

      if (++idx >= MAG_MAX_TRACK)
        {
          fprintf(stderr, "mag_compute_F1: MAG_MAX_TRACK too small\n");
          return -1;
        }
    }

  track->n = idx;

  return s;
} /* mag_compute_F1() */

/*
mag_output()
  After a track has been successfully inverted for the EEJ,
write the EEJ current density to an output file with the format:

EEJ-YYYYMMDD-HHMMSS.txt
*/

static int
mag_output(const int header, const mag_workspace *w)
{
  FILE *fp = w->fp_output;
  char str[EPOCH_STRING_LEN + 1];
  double t = w->track.t_eq;
  size_t i;
  long year, month, day, hour, min, sec, msec;

  if (w->params->output_file == NULL)
    return 0; /* no output desired */

  if (header)
    {
      /* print header */
      fprintf(fp, "# Field 1: date (YYYY-MM-DD)\n");
      fprintf(fp, "# Field 2: universal time of equator crossing (HH:MM:SS)\n");
      fprintf(fp, "# Field 3: longitude of equator crossing (degrees)\n");
      fprintf(fp, "# Field 4-64: EEJ height-integrated current density (A/m)\n");
      return 0;
    }

  EPOCHbreakdown(t, &year, &month, &day, &hour, &min, &sec, &msec);

  sprintf(str, "%04ld-%02ld-%02ld %02ld:%02ld:%02ld",
          year,
          month,
          day,
          hour,
          min,
          sec);

  fprintf(fp, "%s %.3f ", str, w->track.phi_eq * 180.0 / M_PI);

  for (i = 0; i < w->ncurr; ++i)
    fprintf(fp, "%f ", w->EEJ[i]);

  fprintf(fp, "\n");
  fflush(fp);

  return GSL_SUCCESS;
} /* mag_output() */

static int
mag_callback_season(const double doy, const void *params)
{
  const mag_params *p = (mag_params *) params;

  if (doy >= p->season_min && doy <= p->season_max)
    return 0;

  if (p->season_min2 >= 0.0 && p->season_max2 >= 0.0 &&
     (doy >= p->season_min2 && doy <= p->season_max2))
    return 0;

  return -1;
}
