/*
 * mfield_preproc.c
 *
 * Pre-process satellite data and save in magdata format. This
 * program should be called with Stage1 data, in order to do track
 * rms test, and along-track gradient calculation
 *
 * Pre-processing steps are:
 * 1. Instrument flags (recommended CHAMP flags except 1 star camera allowed)
 * 2. Track rms test
 * 3. filter for WMM criteria
 * 4. Downsample by factor 20
 * 5. Compute and store along-track gradients
 *
 * The result is an output file in magdata format containing all data point to
 * be used in the modeling. All data points will have a MAGDATA_FLG_FIT_xxx flag
 * set, and other flags will vary.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <satdata/satdata.h>
#include <flow/flow.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>

#include "bin3d.h"
#include "common.h"
#include "euler.h"
#include "magdata.h"
#include "track.h"
#include "solarpos.h"

#include "mfield.h"
#include "mfield_data.h"

static size_t model_flags(const size_t magdata_flags, const double t,
                          const double theta, const double phi, const double qdlat);
size_t flag_local_time(const double lt_min, const double lt_max,
                       const double euler_lt_min, const double euler_lt_max,
                       satdata_mag *data);

/* local-time range for main field modeling */
#define MFIELD_LT_MIN              (5.0)
#define MFIELD_LT_MAX              (22.0)

/* local-time range for Euler angle fitting */
#define MFIELD_EULER_LT_MIN        (6.0)
#define MFIELD_EULER_LT_MAX        (18.0)

/* number of seconds for computing along-track differences */
#define MFIELD_GRAD_DT             (40.0)

/* maximum zenith angle in degrees at high latitudes */
#define MFIELD_MAX_ZENITH          (100.0)

/*
 * include apriori crustal field in modeling - set this if fitting only
 * a main field model, so that an apriori crustal field (MF7) is subtracted
 * from the satellite observations
 */
#define MFIELD_INC_CRUSTAL         0

/* maximum QD latitude for fitting vector data */
#define MFIELD_VECTOR_QDLAT        (55.0)

/* maximum QD latitude for fitting Euler angles */
#define MFIELD_EULER_QDLAT         MFIELD_VECTOR_QDLAT

/* geocentric latitude cutoff for using local time vs zenith selection */
#define MFIELD_HIGH_LATITUDE      (60.0)

/* Global */
solarpos_workspace *solarpos_workspace_p = NULL;

typedef struct
{
  int flag_rms;       /* use track rms test */
  size_t downsample;  /* downsampling factor */
} preprocess_parameters;

#if 0
/*
copy_track()
  Copy a single track to magdata structure

Inputs: track_idx - track index
        data      - satellite data
        track_p   - track data
        mdata     - (output) where to store new data
        ndata     - (output) updated to count scalar/vector/euler/gradient
                    data
                    ndata[0] = # scalar measurements for MF modeling
                    ndata[1] = # vector measurements for MF modeling
                    ndata[2] = # vector measurements for Euler angles
                    ndata[3] = # along-track gradient measurements
*/

size_t
copy_track(const size_t track_idx, const satdata_mag *data,
           const track_workspace *track_p, magdata *mdata,
           size_t npts[4])
{
  int s = 0;
  size_t i, k;
  track_data *tptr = &(track_p->tracks[track_idx]);
  const size_t start_idx = tptr->start_idx;
  const size_t end_idx = tptr->end_idx;
  const size_t grad_idx = (size_t) MFIELD_GRAD_DT;
  const double grad_dt_min = MFIELD_GRAD_DT - 2.0;
  const double grad_dt_max = MFIELD_GRAD_DT + 2.0;
  const size_t grad_mask = 0xffffffff & ~SATDATA_FLG_DOWNSAMPLE;
  magdata_datum datum;
  size_t ngrad = 0;
  int flagged_start = 0;

  for (i = start_idx; i <= end_idx; ++i)
    {
      size_t flags = 0;     /* final magdata flags */
      size_t fitting_flags;

      if (SATDATA_BadData(data->flags[i]))
        continue;

      if (data->flags[i] & SATDATA_FLG_DOWNSAMPLE)
        continue;

      /*
       * determine whether this data point will be used for MF coefficients,
       * Euler angles, or both
       */
      fitting_flags = model_flags(mdata->global_flags, i, data);
      if (fitting_flags == 0)
        continue; /* this point will not be used in the modeling */

      flags |= fitting_flags;

      /* initialize to 0 */
      magdata_datum_init(&datum);

      if (!flagged_start)
        {
          flags |= MAGDATA_FLG_TRACK_START;
          flagged_start = 1;
        }

      datum.B_nec[0] = SATDATA_VEC_X(data->B, i);
      datum.B_nec[1] = SATDATA_VEC_Y(data->B, i);
      datum.B_nec[2] = SATDATA_VEC_Z(data->B, i);
      datum.B_vfm[0] = SATDATA_VEC_X(data->B_VFM, i);
      datum.B_vfm[1] = SATDATA_VEC_Y(data->B_VFM, i);
      datum.B_vfm[2] = SATDATA_VEC_Z(data->B_VFM, i);
      datum.F = data->F[i];

      for (k = 0; k < 4; ++k)
        datum.q[k] = data->q[4 * i + k];

      /* subtract external field model from measurements */
      datum.B_model[0] = SATDATA_VEC_X(data->B_ext, i);
      datum.B_model[1] = SATDATA_VEC_Y(data->B_ext, i);
      datum.B_model[2] = SATDATA_VEC_Z(data->B_ext, i);

#if MFIELD_INC_CRUSTAL
      /* subtract crustal field also */
      datum.B_model[0] += SATDATA_VEC_X(data->B_crust, i);
      datum.B_model[1] += SATDATA_VEC_Y(data->B_crust, i);
      datum.B_model[2] += SATDATA_VEC_Z(data->B_crust, i);
#endif

      {
        size_t j = GSL_MIN(i + grad_idx, data->n - 1);
        double dt = (data->t[j] - data->t[i]) / 1000.0; /* in s */

        /*
         * if data j is flagged due to downsampling, that's ok, but reject
         * point if any other flags are set
         */
        if (!(data->flags[j] & grad_mask) &&
            (dt >= grad_dt_min && dt <= grad_dt_max))
          {
            assert(data->flags[j] == 0 || data->flags[j] == SATDATA_FLG_DOWNSAMPLE);

            datum.B_nec_ns[0] = SATDATA_VEC_X(data->B, j);
            datum.B_nec_ns[1] = SATDATA_VEC_Y(data->B, j);
            datum.B_nec_ns[2] = SATDATA_VEC_Z(data->B, j);
            datum.F_ns = data->F[j];

            /* subtract external field model from measurements */
            datum.B_model_ns[0] = SATDATA_VEC_X(data->B_ext, j);
            datum.B_model_ns[1] = SATDATA_VEC_Y(data->B_ext, j);
            datum.B_model_ns[2] = SATDATA_VEC_Z(data->B_ext, j);

#if MFIELD_INC_CRUSTAL
            datum.B_model_ns[0] += SATDATA_VEC_X(data->B_crust, j);
            datum.B_model_ns[1] += SATDATA_VEC_Y(data->B_crust, j);
            datum.B_model_ns[2] += SATDATA_VEC_Z(data->B_crust, j);
#endif

            datum.r_ns = data->altitude[j] + data->R;
            datum.theta_ns = M_PI / 2.0 - data->latitude[j] * M_PI / 180.0;
            datum.phi_ns = data->longitude[j] * M_PI / 180.0;
            datum.qdlat_ns = data->qdlat[j];

            /* set gradient flag to indicate gradient information available */
            flags |= MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS;
          }
      }

      /* scalar measurement available (only for MF fitting) */
      if (fitting_flags & MAGDATA_FLG_FIT_MF)
        flags |= MAGDATA_FLG_F;

      /* vector measurement available for both MF and Euler (if below MFIELD_VECTOR_QDLAT) */
      if (fabs(data->qdlat[i]) <= MFIELD_VECTOR_QDLAT)
        flags |= MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z;

      if (!(flags & (MAGDATA_FLG_F|MAGDATA_FLG_X|MAGDATA_FLG_Y|MAGDATA_FLG_Z)))
        {
          fprintf(stderr, "UH OH\n");
        }

      datum.t = data->t[i];
      datum.r = data->altitude[i] + data->R;
      datum.theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      datum.phi = data->longitude[i] * M_PI / 180.0;
      datum.qdlat = data->qdlat[i];
      datum.satdir = satdata_mag_satdir(i, data);
      datum.flags = flags;

      s = magdata_add(&datum, mdata);
      if (s == 0)
        {
          if (flags & MAGDATA_FLG_FIT_MF)
            {
              if (flags & MAGDATA_FLG_F)
                ++(npts[0]);

              if (flags & (MAGDATA_FLG_X|MAGDATA_FLG_Y|MAGDATA_FLG_Z))
                ++(npts[1]);
            }

          if (flags & MAGDATA_FLG_FIT_EULER)
            {
              if (flags & (MAGDATA_FLG_X|MAGDATA_FLG_Y|MAGDATA_FLG_Z))
                ++(npts[2]);
            }

          if (flags & MAGDATA_FLG_DZ_NS)
            ++npts[3];
        }

      if (s == 0 && flags & MAGDATA_FLG_DZ_NS)
        ++ngrad;
    }

  return ngrad;
} /* copy_track() */

#endif

magdata *
copy_data(const size_t magdata_flags, const satdata_mag *data, const track_workspace *track_p)
{
  const size_t nflagged = satdata_nflagged(data);
  size_t ndata = data->n - nflagged;
  magdata *mdata;
  magdata_params params;
  size_t i;
  size_t npts[4] = { 0, 0, 0, 0 };
  size_t nmodel[4] = { 0, 0, 0, 0 };

  params.grad_dt_ns = MFIELD_GRAD_DT;
  params.model_main = 0;
#if MFIELD_INC_CRUSTAL
  /* subtract MF7 crustal field from data prior to modeling */
  params.model_crust = 1;
#else
  params.model_crust = 0;
#endif
  /* subtract external field from data prior to modeling */
  params.model_ext = 1;

  /*
   * for some reason, the loop below adds more than 'ndata'
   * points into mdata, it shouldn't do this but I'm not sure why
   * so add 10000
   */
  mdata = magdata_alloc(ndata + 10000, data->R);
  if (!mdata)
    return 0;

  mdata->global_flags = magdata_flags;

  /* copy tracks into mdata structure */
  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);

      /* discard flagged tracks */
      if (tptr->flags)
        continue;

#if 0
      ngrad += copy_track(i, data, track_p, mdata, npts);
#else
      magdata_copy_track(&params, i, data, track_p, mdata, npts);
#endif
    }

  /*
   * now determine which points in mdata will be used for
   * MF modeling, Euler angle fitting, etc
   */
  for (i = 0; i < mdata->n; ++i)
    {
      size_t fitting_flags = model_flags(mdata->global_flags,
                                         mdata->t[i], mdata->theta[i],
                                         mdata->phi[i], mdata->qdlat[i]);

      if ((fitting_flags & MAGDATA_FLG_FIT_MF) &&
          MAGDATA_ExistScalar(mdata->flags[i]))
        ++(nmodel[0]);
      if ((fitting_flags & MAGDATA_FLG_FIT_MF) &&
          MAGDATA_ExistVector(mdata->flags[i]))
        ++(nmodel[1]);
      if ((fitting_flags & MAGDATA_FLG_FIT_EULER) &&
          MAGDATA_ExistVector(mdata->flags[i]))
        ++(nmodel[2]);

      mdata->flags[i] |= fitting_flags;
    }

#if 0
  fprintf(stderr, "\n");
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) scalar measurements for MF modeling\n",
          npts[0], mdata->n, (double) npts[0] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) vector measurements for MF modeling\n",
          npts[1], mdata->n, (double) npts[1] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) vector measurements for Euler angle modeling\n",
          npts[2], mdata->n, (double) npts[2] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) along-track gradient data available\n",
          ngrad, ndata, (double)ngrad / (double)ndata * 100.0);
#else
  fprintf(stderr, "\n");
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) scalar measurements available\n",
          npts[0], mdata->n, (double) npts[0] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) vector measurements available\n",
          npts[1], mdata->n, (double) npts[1] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) along-track scalar measurements available\n",
          npts[2], mdata->n, (double) npts[2] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) along-track vector measurements available\n",
          npts[3], mdata->n, (double) npts[3] / (double) mdata->n * 100.0);

  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) scalar measurements selected for MF modeling\n",
          nmodel[0], mdata->n, (double) nmodel[0] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) vector measurements selected for MF modeling\n",
          nmodel[1], mdata->n, (double) nmodel[1] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) vector measurements selected for Euler angle modeling\n",
          nmodel[2], mdata->n, (double) nmodel[2] / (double) mdata->n * 100.0);
#endif

  return mdata;
}

satdata_mag *
read_swarm(const char *filename, const int asmv)
{
  size_t nflag;
  satdata_mag *data;
  struct timeval tv0, tv1;

  fprintf(stderr, "read_swarm: reading %s...", optarg);
  gettimeofday(&tv0, NULL);

  if (asmv)
    data = satdata_swarm_asmv_read_idx(optarg, 0);
  else
    data = satdata_swarm_read_idx(optarg, 0);

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu data read, %g seconds)\n",
          data->n, time_diff(tv0, tv1));

  /* check for instrument flags since we use Stage1 data */

  fprintf(stderr, "read_swarm: filtering for instrument flags...");
  nflag = satdata_swarm_filter_instrument(1, data);
  fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
          nflag, data->n, (double)nflag / (double)data->n * 100.0);

  return data;
} /* read_swarm() */

satdata_mag *
read_champ(const char *filename)
{
  satdata_mag *data;
  struct timeval tv0, tv1;
  size_t nflag;

  fprintf(stderr, "read_champ: reading %s...", optarg);
  gettimeofday(&tv0, NULL);
  data = satdata_champ_read_idx(optarg, 0);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu data read, %g seconds)\n",
          data->n, time_diff(tv0, tv1));

  /* check for instrument flags since we use Stage1 data */
  fprintf(stderr, "read_champ: filtering for instrument flags...");
  nflag = satdata_champ_filter_instrument(1, 0, data);
  fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
          nflag, data->n, (double)nflag / (double)data->n * 100.0);

  return data;
} /* read_champ() */

/*
check_data_point()
  Check an individual data point to determine if it will be used to
fit various model parameters.

Inputs: magdata_flags - MAGDATA_GLOBFLG_xxx
        t             - CDF_EPOCH timestamp
        theta         - colatitude (radians)
        phi           - longitude (radians)
        qdlat         - QD latitude (degrees)

Return: flags indicating fit parameters (MAGDATA_FLG_FIT_xxx)

Notes:
1) If data point is below MFIELD_HIGH_LATITUDE and local time is
within [MFIELD_LT_MIN,MFIELD_LT_MAX], flag is set to MAGDATA_FLG_FIT_MF

2) If data point is higher than MFIELD_HIGH_LATITUDE, zenith angle
is computed and if the point is in darkness, MAGDATA_FLG_FIT_MF is set

3) If we are fitting Euler angles to this satellite, and the data
point satisfies the criteria, MAGDATA_FLG_FIT_EULER is set
*/

static size_t
model_flags(const size_t magdata_flags, const double t,
            const double theta, const double phi, const double qdlat)
{
  size_t flags = 0;
  const time_t unix_time = satdata_epoch2timet(t);
  const double lt = get_localtime(unix_time, phi);
  const double lat_deg = 90.0 - theta * 180.0 / M_PI;

  /* check if we should fit Euler angles to this data point */
  if (magdata_flags & MAGDATA_GLOBFLG_EULER &&
      fabs(qdlat) <= MFIELD_EULER_QDLAT &&
      (lt >= MFIELD_EULER_LT_MAX || lt <= MFIELD_EULER_LT_MIN))
    {
      flags |= MAGDATA_FLG_FIT_EULER;
    }

  /* check if we should fit main field model to this data point */
  if (fabs(lat_deg) <= MFIELD_HIGH_LATITUDE)
    {
      if (lt >= MFIELD_LT_MAX || lt <= MFIELD_LT_MIN)
        flags |= MAGDATA_FLG_FIT_MF;
    }
  else
    {
      double lat_rad = lat_deg * M_PI / 180.0;
      double zenith;

      solarpos_calc_zenith(unix_time, lat_rad, phi, &zenith, solarpos_workspace_p);
      zenith *= 180.0 / M_PI;
      assert(zenith >= 0.0);

      /* large zenith angle means darkness */
      if (zenith >= MFIELD_MAX_ZENITH)
        flags |= MAGDATA_FLG_FIT_MF;
    }

  return flags;
}

/*
flag_local_time()

Inputs: lt_min       - min LT for modeling
        lt_max       - max LT for modeling
        euler_lt_min - min LT for Euler angles
        euler_lt_max - max LT for Euler angles
        data         - satellite data
*/

size_t
flag_local_time(const double lt_min, const double lt_max,
                const double euler_lt_min, const double euler_lt_max,
                satdata_mag *data)
{
  size_t i;
  size_t cnt = 0;
  solarpos_workspace *work_sp = solarpos_alloc();

  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time;
      double lt;

      /* ignore already flagged data to improve speed */
      if (data->flags[i])
        continue;

      unix_time = satdata_epoch2timet(data->t[i]);
      lt = get_localtime(unix_time, data->longitude[i]*M_PI/180.0);

#if 0 /* XXX */
      /* check if data point should be used for Euler angles */
      if (lt <= euler_lt_min || lt >= euler_lt_max)
        data->flags[i] |= SATDATA_FLG_EULER;
#endif

      /*
       * at high latitudes, calculate solar zenith angle for a better
       * measure of darkness/sunlight than local time
       */
      if (fabs(data->qdlat[i]) > 60.0)
        {
          double lat = data->latitude[i] * M_PI / 180.0;
          double lon = wrappi(data->longitude[i] * M_PI / 180.0);
          double zenith;

          solarpos_calc_zenith(unix_time, lat, lon, &zenith, work_sp);
          zenith *= 180.0 / M_PI;
          assert(zenith >= 0.0);

          /* small zenith angle means sunlight so flag it */
          if (zenith < MFIELD_MAX_ZENITH)
            {
              ++cnt;
              data->flags[i] |= SATDATA_FLG_LT;
            }
        }
      else if (lt >= lt_min && lt <= lt_max)
        {
          ++cnt;
          data->flags[i] |= SATDATA_FLG_LT;
        }
    }

  solarpos_free(work_sp);

  return cnt;
}

int
print_track_stats(const satdata_mag *data, const track_workspace *track_p)
{
  size_t nflagged = satdata_nflagged(data);
  size_t nleft = data->n - nflagged;
  size_t nflagged_track = track_nflagged(track_p);
  size_t nleft_track = track_p->n - nflagged_track;

  fprintf(stderr, "preprocess_data: total flagged data: %zu/%zu (%.1f%%)\n",
          nflagged, data->n, (double)nflagged / (double)data->n * 100.0);
  fprintf(stderr, "preprocess_data: total remaining data: %zu/%zu (%.1f%%)\n",
          nleft, data->n, (double)nleft / (double)data->n * 100.0);

  fprintf(stderr, "preprocess_data: total flagged tracks: %zu/%zu (%.1f%%)\n",
          nflagged_track, track_p->n, (double)nflagged_track / (double)track_p->n * 100.0);
  fprintf(stderr, "preprocess_data: total remaining tracks: %zu/%zu (%.1f%%)\n",
          nleft_track, track_p->n, (double)nleft_track / (double)track_p->n * 100.0);

  return 0;
} /* print_track_stats() */

/*
preprocess_data()

Inputs: params - preprocess parameters
          downsample - downsampling factor
        data   - satellite data

Return: pointer to sorted track workspace (should be freed by caller)
*/

track_workspace *
preprocess_data(const preprocess_parameters *params, const size_t magdata_flags,
                satdata_mag *data)
{
  struct timeval tv0, tv1;
  track_workspace *track_p = track_alloc();

  fprintf(stderr, "preprocess_data: initializing tracks...");
  gettimeofday(&tv0, NULL);
  track_init(data, NULL, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* detect bad tracks with rms test */
  if (params->flag_rms)
    {
      const char *rmsfile = "satrms.dat";
      size_t nrms;
      double thresh[] = { 20.0, 25.0, 15.0, 15.0 };

      nrms = track_flag_rms(rmsfile, thresh, data, track_p);
      fprintf(stderr, "preprocess_data: flagged (%zu/%zu) (%.1f%%) points due to high rms\n",
              nrms, data->n, (double) nrms / (double) data->n * 100.0);
    }

  /* flag according to WMM criteria */
  satdata_filter_wmm(1, data);

  /* downsample data */
  {
    size_t i;

    fprintf(stderr, "preprocess_data: downsampling data by factor %zu...", params->downsample);

    for (i = 0; i < data->n; ++i)
      {
        if (i % params->downsample != 0)
          data->flags[i] |= SATDATA_FLG_DOWNSAMPLE;
      }

    fprintf(stderr, "done\n");
  }

  {
    size_t i;
    size_t nflag = 0;

    fprintf(stderr, "preprocess_data: flagging points due to modeling criteria...");
    for (i = 0; i < data->n; ++i)
      {
        double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
        double phi = data->longitude[i] * M_PI / 180.0;
        size_t fitting_flags = model_flags(magdata_flags,
                                           data->t[i], theta,
                                           phi, data->qdlat[i]);

        /* flag point if it will not be used in any modeling */
        if (fitting_flags == 0)
          {
            data->flags[i] |= SATDATA_FLG_OUTLIER;
            ++nflag;
          }
      }
    fprintf(stderr, "done (%zu/%zu (%.1f%%) points flagged)\n",
            nflag, data->n, (double)nflag / (double)data->n * 100.0);
  }

  /* print track statistics */
  {
    char *stat_file = "track_stats.dat";

    fprintf(stderr, "preprocess_data: printing track statistics to %s...", stat_file);
    track_print_stats(stat_file, track_p);
    fprintf(stderr, "done\n");
  }

  print_track_stats(data, track_p);

  return track_p;
} /* preprocess_data() */

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --champ_file      | -c champ_index_file       - CHAMP index file\n");
  fprintf(stderr, "\t --swarm_file      | -s swarm_index_file       - Swarm index file\n");
  fprintf(stderr, "\t --swarm_asmv_file | -a swarm_asmv_index_file  - Swarm ASM-V index file\n");
  fprintf(stderr, "\t --downsample      | -d downsample             - downsampling factor\n");
  fprintf(stderr, "\t --euler_file      | -e euler_file             - Euler angles file\n");
  fprintf(stderr, "\t --output_file     | -o output_file            - output file\n");
}

int
main(int argc, char *argv[])
{
  char *datamap_file = "datamap.dat";
  char *data_file = "data.dat";
  char *output_file = NULL;
  satdata_mag *data = NULL;
  magdata *mdata;
  euler_workspace *euler_p = NULL;
  struct timeval tv0, tv1;
  track_workspace *track_p;
  preprocess_parameters params;
  size_t magdata_flags = 0;

  /* defaults */
  params.flag_rms = 1;
  params.downsample = 20;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "swarm_file", required_argument, NULL, 's' },
          { "swarm_asmv_file", required_argument, NULL, 'a' },
          { "champ_file", required_argument, NULL, 'c' },
          { "downsample", required_argument, NULL, 'd' },
          { "output_file", required_argument, NULL, 'o' },
          { "euler_file", required_argument, NULL, 'e' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:c:d:e:o:s:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          /* Swarm ASM-V */
          case 'a':
            data = read_swarm(optarg, 1);
            magdata_flags = MAGDATA_GLOBFLG_EULER;
            params.flag_rms = 0; /* no NEC data for rms flagging */
            break;

          /* Swarm official */
          case 's':
            data = read_swarm(optarg, 0);
            magdata_flags = MAGDATA_GLOBFLG_EULER;
            break;

          case 'c':
            data = read_champ(optarg);
            break;

          case 'd':
            params.downsample = (size_t) atoi(optarg);
            break;

          case 'e':
              fprintf(stderr, "main: reading Euler angles from %s...", optarg);
              euler_p = euler_read(optarg);
              if (!euler_p)
                exit(1);
              fprintf(stderr, "done (%zu sets of angles read)\n", euler_p->n);
              break;

          case 'o':
            output_file = optarg;
            break;

          default:
            break;
        }
    }

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  solarpos_workspace_p = solarpos_alloc();

  fprintf(stderr, "main: LT minimum       = %.1f\n", MFIELD_LT_MIN);
  fprintf(stderr, "main: LT maximum       = %.1f\n", MFIELD_LT_MAX);
  fprintf(stderr, "main: Euler LT minimum = %.1f\n", MFIELD_EULER_LT_MIN);
  fprintf(stderr, "main: Euler LT maximum = %.1f\n", MFIELD_EULER_LT_MAX);

  if (euler_p)
    {
      fprintf(stderr, "main: rotating VFM measurements with new Euler angles...");
      euler_apply(data, euler_p);
      fprintf(stderr, "done\n");
    }

  track_p = preprocess_data(&params, magdata_flags, data);

  fprintf(stderr, "main: computing vector residuals...");
  gettimeofday(&tv0, NULL);
  mdata = copy_data(magdata_flags, data, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu residuals copied, %g seconds)\n",
          mdata->n, time_diff(tv0, tv1));

  /* free data after copying arrays to free up memory */
  satdata_mag_free(data);

  magdata_init(mdata);

  fprintf(stderr, "main: computing spatial weighting of data...");
  gettimeofday(&tv0, NULL);
  magdata_calc(mdata);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: writing data to %s...", data_file);
  magdata_print(data_file, mdata);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing data map to %s...", datamap_file);
  magdata_map(datamap_file, mdata);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: satellite rmin = %.1f (%.1f) [km]\n",
          mdata->rmin, mdata->rmin - mdata->R);
  fprintf(stderr, "main: satellite rmax = %.1f (%.1f) [km]\n",
          mdata->rmax, mdata->rmax - mdata->R);

  if (output_file)
    {
      fprintf(stderr, "main: writing data to %s...", output_file);
      magdata_write(output_file, mdata);
      fprintf(stderr, "done\n");
    }

  magdata_free(mdata);
  track_free(track_p);
  solarpos_free(solarpos_workspace_p);

  return 0;
} /* main() */
