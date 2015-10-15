/*
 * track_synth.c
 *
 * Synthesize main, crustal, external field values along
 * satellite track
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <indices/indices.h>
#include <satdata/satdata.h>
#include <gsl/gsl_math.h>

#include "apex.h"
#include "chaos.h"
#include "common.h"
#include "msynth.h"
#include "pomme.h"

#include "track.h"

/*
track_synth()
  Synthesize main, crustal and external fields along satellite track.
Main field is provided as function argument. Crustal field is MF7.
External field is POMME

Inputs: down_sample - number of samples to throw out (>= 1)
                      (ie: if this is 5, every 5th sample is kept and
                       the rest discarded)
        data_in     - satellite data input
        data_out    - satellite data output
        msynth_core - msynth core workspace (degrees 1 to 15)
*/

int
track_synth(const int down_sample, const satdata_mag *data_in,
            satdata_mag *data_out, msynth_workspace *msynth_core)
{
  int s = 0;
  size_t i, j;
  double B_core[4], B_crust[4], B_int[4], B_ext[4], F_pomme, dF_ext;
  msynth_workspace *msynth_crust = msynth_mf7_read(MSYNTH_MF7_FILE);
  pomme_workspace *pomme_p = pomme_alloc_default();
  int year = (int) satdata_epoch2year(data_in->t[0]);
  apex_workspace *apex_p = apex_alloc(year);
  size_t idx = 0; /* data_out index */
  int prev_satdir = 0;
  double E_st = 0.0, I_st = 0.0;
  double IMF_By = 0.0, Em = 0.5;
  double f107 = 120.0;
  estist_workspace *estist_workspace_p = estist_alloc(ESTIST_IDX_FILE);

  pomme_set_radius(R_EARTH_KM, pomme_p);

  msynth_set(1, 15, msynth_core);
  msynth_set(16, 133, msynth_crust);

  for (i = 0; i < data_in->n; i += down_sample)
    {
      time_t t = satdata_epoch2timet(data_in->t[i]);
      double tyr = satdata_epoch2year(data_in->t[i]);
      double theta = M_PI / 2.0 - data_in->latitude[i] * M_PI / 180.0;
      double phi = data_in->longitude[i] * M_PI / 180.0;
      double alt = data_in->altitude[i];
      int satdir = satdata_mag_satdir(i, data_in);
      double alon, alat, qdlat;

      if (satdir != prev_satdir)
        {
          /*
           * satellite changed flight direction, new track; re-compute
           * indices for this track; this helps ensure POMME external
           * model is smooth/continuous over each track
           */
          s = pomme_get_indices(t, &E_st, &I_st, &IMF_By, &Em, &f107, pomme_p);
          if (s)
            return s;

          prev_satdir = satdir;
        }

      /* copy the entire data record to output structure */
      satdata_mag_copy_record(i, data_in, idx, data_out);

      /* compute QD latitude */
      apex_transform(theta, phi, (data_in->R + alt) * 1.0e3,
                     &alon, &alat, &qdlat, NULL, NULL, NULL,
                     apex_p);

      /* compute internal and external fields */
      msynth_eval(tyr, alt + data_in->R, theta, phi, B_core, msynth_core);
      msynth_eval(tyr, alt + data_in->R, theta, phi, B_crust, msynth_crust);

      for (j = 0; j < 3; ++j)
        B_int[j] = B_core[j] + B_crust[j];

      B_int[3] = gsl_hypot3(B_int[0], B_int[1], B_int[2]);

      /* Est/Ist now interpolates so we can calculate it at each point */
      estist_get(t, &E_st, &I_st, estist_workspace_p);
      s = pomme_calc_ext_indices(theta, phi, t, data_in->altitude[i],
                                 E_st, I_st, IMF_By, Em, f107, B_ext, pomme_p);
      if (s)
        return s;

      /* convert to nT */
      for (j = 0; j < 3; ++j)
        B_ext[j] *= 1.0e9;

      /* compute dF_ext = B_ext . b_int */
      dF_ext = 0.0;
      for (j = 0; j < 3; ++j)
        dF_ext += B_int[j] * B_ext[j] / B_int[3];

      F_pomme = B_int[3] + dF_ext;

      /* store scalar field (internal + external) */
      data_out->F_main[idx] = F_pomme;

      /* store vector core field */
      SATDATA_VEC_X(data_out->B_main, idx) = B_core[0];
      SATDATA_VEC_Y(data_out->B_main, idx) = B_core[1];
      SATDATA_VEC_Z(data_out->B_main, idx) = B_core[2];

      /* store vector crustal field */
      SATDATA_VEC_X(data_out->B_crust, idx) = B_crust[0];
      SATDATA_VEC_Y(data_out->B_crust, idx) = B_crust[1];
      SATDATA_VEC_Z(data_out->B_crust, idx) = B_crust[2];

      /* store vector external field */
      SATDATA_VEC_X(data_out->B_ext, idx) = B_ext[0];
      SATDATA_VEC_Y(data_out->B_ext, idx) = B_ext[1];
      SATDATA_VEC_Z(data_out->B_ext, idx) = B_ext[2];

      /* store QD latitude */
      data_out->qdlat[idx] = qdlat;

      ++idx;
    }

  data_out->n = idx;
  data_out->R = data_in->R;

  pomme_free(pomme_p);
  msynth_free(msynth_crust);
  apex_free(apex_p);
  estist_free(estist_workspace_p);

  return s;
} /* track_synth() */

/*
track_synth_chaos()
  Synthesize main, crustal and external fields along satellite track.
Main field is provided as function argument. Crustal field is MF7.
External field is CHAOS

Inputs: down_sample - number of samples to throw out (>= 1)
                      (ie: if this is 5, every 5th sample is kept and
                       the rest discarded)
        data_in     - satellite data input
        data_out    - satellite data output
        msynth_core - msynth core workspace (degrees 1 to 15)
*/

int
track_synth_chaos(const int down_sample, const satdata_mag *data_in,
                  satdata_mag *data_out, msynth_workspace *msynth_core)
{
  int s = 0;
  const size_t n = data_in->n;
  size_t i, j;
  double B_core[4], B_crust[4], B_int[4], B_ext[4], F_pomme, dF_ext;
  msynth_workspace *msynth_crust = msynth_mf7_read(MSYNTH_MF7_FILE);
  int year = (int) satdata_epoch2year(data_in->t[0]);
  apex_workspace *apex_p = apex_alloc(year);
  size_t idx = 0; /* data_out index */
  rc_workspace *rc_workspace_p = rc_alloc(RC_IDX_FILE);

  msynth_set(1, 15, msynth_core);
  msynth_set(16, 133, msynth_crust);

  /* compute CHAOS external field first with single call to MATLAB */
  {
    gsl_vector_view tv = gsl_vector_view_array(data_in->t, n);
    gsl_vector_view altv = gsl_vector_view_array(data_in->altitude, n);
    gsl_vector_view latv = gsl_vector_view_array(data_in->latitude, n);
    gsl_vector_view lonv = gsl_vector_view_array(data_in->longitude, n);
    gsl_vector *RC_ext = gsl_vector_alloc(n);
    gsl_vector *RC_ind = gsl_vector_alloc(n);
    gsl_matrix_view B_extv = gsl_matrix_view_array(data_in->B_ext, n, 3);

    /* store RC indices */
    for (i = 0; i < n; ++i)
      {
        time_t t = satdata_epoch2timet(data_in->t[i]);
        double RC_e, RC_i;

        s = rc_get(t, &RC_e, &RC_i, rc_workspace_p);
        if (s)
          return s;

        gsl_vector_set(RC_ext, i, RC_e);
        gsl_vector_set(RC_ind, i, RC_i);
      }

    chaos_ext(&tv.vector, &altv.vector, &latv.vector, &lonv.vector, RC_ext, RC_ind,
              &B_extv.matrix);

    gsl_vector_free(RC_ext);
    gsl_vector_free(RC_ind);
  }

  for (i = 0; i < data_in->n; i += down_sample)
    {
      double tyr = satdata_epoch2year(data_in->t[i]);
      double r = data_in->R + data_in->altitude[i];
      double theta = M_PI / 2.0 - data_in->latitude[i] * M_PI / 180.0;
      double phi = data_in->longitude[i] * M_PI / 180.0;
      double alon, alat, qdlat;

      /* copy the entire data record to output structure */
      satdata_mag_copy_record(i, data_in, idx, data_out);

      /* compute QD latitude */
      apex_transform(theta, phi, r * 1.0e3,
                     &alon, &alat, &qdlat, NULL, NULL, NULL,
                     apex_p);

      /* compute internal and external fields */
      msynth_eval(tyr, r, theta, phi, B_core, msynth_core);
      msynth_eval(tyr, r, theta, phi, B_crust, msynth_crust);

      for (j = 0; j < 3; ++j)
        B_int[j] = B_core[j] + B_crust[j];

      B_int[3] = gsl_hypot3(B_int[0], B_int[1], B_int[2]);

      /* get external field from MATLAB */
      B_ext[0] = SATDATA_VEC_X(data_in->B_ext, i);
      B_ext[1] = SATDATA_VEC_Y(data_in->B_ext, i);
      B_ext[2] = SATDATA_VEC_Z(data_in->B_ext, i);

      /* compute dF_ext = B_ext . b_int */
      dF_ext = 0.0;
      for (j = 0; j < 3; ++j)
        dF_ext += B_int[j] * B_ext[j] / B_int[3];

      F_pomme = B_int[3] + dF_ext;

      /* store scalar field (internal + external) */
      data_out->F_main[idx] = F_pomme;

      /* store vector core field */
      SATDATA_VEC_X(data_out->B_main, idx) = B_core[0];
      SATDATA_VEC_Y(data_out->B_main, idx) = B_core[1];
      SATDATA_VEC_Z(data_out->B_main, idx) = B_core[2];

      /* store vector crustal field */
      SATDATA_VEC_X(data_out->B_crust, idx) = B_crust[0];
      SATDATA_VEC_Y(data_out->B_crust, idx) = B_crust[1];
      SATDATA_VEC_Z(data_out->B_crust, idx) = B_crust[2];

      /* store vector external field */
      SATDATA_VEC_X(data_out->B_ext, idx) = B_ext[0];
      SATDATA_VEC_Y(data_out->B_ext, idx) = B_ext[1];
      SATDATA_VEC_Z(data_out->B_ext, idx) = B_ext[2];

      /* store QD latitude */
      data_out->qdlat[idx] = qdlat;

      ++idx;
    }

  data_out->n = idx;
  data_out->R = data_in->R;

  msynth_free(msynth_crust);
  apex_free(apex_p);
  rc_free(rc_workspace_p);

  return s;
} /* track_synth_chaos() */
