/*
 * track_synth.c
 *
 * Synthesize main, crustal, external field values along
 * satellite track
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

#include <indices/indices.h>
#include <satdata/satdata.h>
#include <gsl/gsl_math.h>

#include "apex.h"
#include "common.h"
#include "msynth.h"
#include "pomme.h"

#include "track.h"

/*
track_synth_int()
  Synthesize main and crustal fields along satellite track. Also
computes QD latitudes.

Inputs: data           - satellite data output
        msynth_core_p  - msynth core workspace (degrees 1 to 15)
        msynth_crust_p - msynth crustal field workspace
*/

int
track_synth_int(satdata_mag *data, msynth_workspace *msynth_core_p, msynth_workspace *msynth_crust_p)
{
  int s = 0;
  size_t i;
  const size_t max_threads = (size_t) omp_get_max_threads();
  msynth_workspace **crust_p = malloc(max_threads * sizeof(msynth_workspace *));
  msynth_workspace **core_p = malloc(max_threads * sizeof(msynth_workspace *));
  int year = (int) satdata_epoch2year(data->t[0]);
  apex_workspace *apex_p = apex_alloc(year);

  for (i = 0; i < max_threads; ++i)
    {
      core_p[i] = msynth_copy(msynth_core_p);
      crust_p[i] = msynth_copy(msynth_crust_p);

      msynth_set(1, 15, core_p[i]);
      msynth_set(16, msynth_crust_p->eval_nmax, crust_p[i]);
    }

#pragma omp parallel for private(i)
  for (i = 0; i < data->n; ++i)
    {
      int thread_id = omp_get_thread_num();
      double tyr = satdata_epoch2year(data->t[i]);
      double r = data->r[i];
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double B_core[4], B_crust[4];

      /* compute internal and external fields */
      msynth_eval(tyr, r, theta, phi, B_core, core_p[thread_id]);
      msynth_eval(tyr, r, theta, phi, B_crust, crust_p[thread_id]);

      /* store vector core field */
      SATDATA_VEC_X(data->B_main, i) = B_core[0];
      SATDATA_VEC_Y(data->B_main, i) = B_core[1];
      SATDATA_VEC_Z(data->B_main, i) = B_core[2];

      /* store vector crustal field */
      SATDATA_VEC_X(data->B_crust, i) = B_crust[0];
      SATDATA_VEC_Y(data->B_crust, i) = B_crust[1];
      SATDATA_VEC_Z(data->B_crust, i) = B_crust[2];
    }

  /* now compute QD latitudes (apex_transform is not thread safe) */
  for (i = 0; i < data->n; ++i)
    {
      double r = data->r[i];
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double alon, alat, qdlat;

      /* compute QD latitude (apex_transform is not thread-safe) */
      apex_transform(theta, phi, r * 1.0e3,
                     &alon, &alat, &qdlat, NULL, NULL, NULL, apex_p);

      /* store QD latitude */
      data->qdlat[i] = qdlat;
    }

  for (i = 0; i < max_threads; ++i)
    {
      msynth_free(core_p[i]);
      msynth_free(crust_p[i]);
    }

  free(crust_p);
  free(core_p);
  apex_free(apex_p);

  return s;
}

/*
track_synth_pomme()
  Synthesize POMME external field along satellite track.

Inputs: data - satellite data input/output

Notes:
1) On output, data->B_ext is filled with POMME values
*/

int
track_synth_pomme(satdata_mag *data)
{
  int s = 0;
  size_t i, j;
  const size_t max_threads = (size_t) omp_get_max_threads();
  pomme_workspace **ext_p = malloc(max_threads * sizeof(pomme_workspace *));
  size_t idx = 0; /* data_out index */
  estist_workspace *estist_workspace_p = estist_alloc(ESTIST_IDX_FILE);

  for (i = 0; i < max_threads; ++i)
    {
      ext_p[i] = pomme_alloc_default();
      pomme_set_radius(R_EARTH_KM, ext_p[i]);
    }

#pragma omp parallel for private(i)
  for (i = 0; i < data->n; ++i)
    {
      int thread_id = omp_get_thread_num();
      time_t t = satdata_epoch2timet(data->t[i]);
      double tyr = satdata_epoch2year(data->t[i]);
      double r = data->r[i];
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double B_ext[4];
      double E_st = 0.0, I_st = 0.0;
      double IMF_By = 0.0, Em = 0.5;
      double f107 = 120.0;

      /*
       * satellite changed flight direction, new track; re-compute
       * indices for this track; this helps ensure POMME external
       * model is smooth/continuous over each track
       */
      s = pomme_get_indices(0, t, &E_st, &I_st, &IMF_By, &Em, &f107, ext_p[thread_id]);
      if (s)
        {
          fprintf(stderr, "track_synth_pomme: pomme_get_indices failed: %s\n", s);
          continue;
        }

      /* Est/Ist now interpolates so we can calculate it at each point */
      estist_get(t, &E_st, &I_st, estist_workspace_p);
      s = pomme_calc_ext_indices(theta, phi, t, data->altitude[i],
                                 E_st, I_st, IMF_By, Em, f107, B_ext, ext_p[thread_id]);
      if (s)
        {
          fprintf(stderr, "track_synth_pomme: pomme_calc_ext_indices failed: %s\n", s);
          continue;
        }

      /* convert to nT */
      for (j = 0; j < 3; ++j)
        B_ext[j] *= 1.0e9;

      /* store vector external field */
      SATDATA_VEC_X(data->B_ext, i) = B_ext[0];
      SATDATA_VEC_Y(data->B_ext, i) = B_ext[1];
      SATDATA_VEC_Z(data->B_ext, i) = B_ext[2];
    }

  for (i = 0; i < max_threads; ++i)
    pomme_free(ext_p[i]);

  free(ext_p);
  estist_free(estist_workspace_p);

  return s;
}
