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
void print_unflagged_data(const char *filename, const satdata_mag *data);

/* use SMDL index for filtering instead of WMM criteria */
#define MFIELD_FILTER_SMDL         0

/*
 * local-time range for main field modeling; the (mod 24) issue is handled
 * in check_lt(), so if you want to select data between 10pm and 5am, set
 * [min,max] = [22,5]
 */
#define MFIELD_LT_MIN              (0.0)
#define MFIELD_LT_MAX              (5.0)

/* local-time range for Euler angle fitting */
#define MFIELD_EULER_LT_MIN        (18.0)
#define MFIELD_EULER_LT_MAX        (6.0)

/* number of seconds for computing along-track differences */
#define MFIELD_GRAD_DT_NS          (40.0)

/* maximum allowed seconds between satellite measurements for computing east-west differences */
#define MFIELD_GRAD_DT_EW          (10.0)

/* maximum allowed longitudinal separation between satellites for computing east-west differences (degrees) */
#define MFIELD_GRAD_DPHI_MAX       (2.0)

/* maximum allowed latitudinal separation between satellites for computing east-west differences (degrees) */
#define MFIELD_GRAD_DLAT_MAX       (0.1)

/* maximum zenith angle in degrees at high latitudes */
#define MFIELD_MAX_ZENITH          (100.0)

/*
 * include apriori crustal field in modeling - set this if fitting only
 * a main field model, so that an apriori crustal field (MF7) is subtracted
 * from the satellite observations
 */
#define MFIELD_INC_CRUSTAL         0

/* QD latitude above which we do not fit X/Y data for MF modeling */
#define MFIELD_QDLAT_HIGH          (55.0)

/* maximum QD latitude for fitting Euler angles */
#define MFIELD_EULER_QDLAT         MFIELD_QDLAT_HIGH

/* geocentric latitude cutoff for using local time vs zenith selection */
#define MFIELD_HIGH_LATITUDE      (60.0)

/* define to fit Z component at high latitudes instead of F */
#define MFIELD_FIT_Z_HIGHLAT      1

#define MFIELD_IDX_X              0
#define MFIELD_IDX_Y              1
#define MFIELD_IDX_Z              2
#define MFIELD_IDX_F              3
#define MFIELD_IDX_DX_NS          4
#define MFIELD_IDX_DY_NS          5
#define MFIELD_IDX_DZ_NS          6
#define MFIELD_IDX_DF_NS          7
#define MFIELD_IDX_DX_EW          8
#define MFIELD_IDX_DY_EW          9
#define MFIELD_IDX_DZ_EW          10
#define MFIELD_IDX_DF_EW          11
#define MFIELD_IDX_B_EULER        12
#define MFIELD_IDX_END            13

/* Global */
solarpos_workspace *solarpos_workspace_p = NULL;

typedef struct
{
  int flag_rms;       /* use track rms test */
  size_t downsample;  /* downsampling factor */
} preprocess_parameters;

magdata *
copy_data(const size_t magdata_flags, const satdata_mag *data, const track_workspace *track_p,
          const size_t magdata_flags2, const satdata_mag *data2, const track_workspace *track_p2)
{
  const size_t nflagged = satdata_nflagged(data);
  size_t ndata = data->n - nflagged;
  magdata *mdata;
  magdata_params params;
  size_t i;
  size_t npts[6] = { 0, 0, 0, 0, 0, 0 };
  size_t nmodel[MFIELD_IDX_END];

  params.grad_dt_ns = MFIELD_GRAD_DT_NS;
  params.grad_dt_ew = MFIELD_GRAD_DT_EW;
  params.grad_dphi_max = MFIELD_GRAD_DPHI_MAX;
  params.grad_dlat_max = MFIELD_GRAD_DLAT_MAX;
  params.model_main = 0;
#if MFIELD_INC_CRUSTAL
  /* subtract MF7 crustal field from data prior to modeling */
  params.model_crust = 1;
#else
  params.model_crust = 0;
#endif
  /* subtract external field from data prior to modeling */
  params.model_ext = 1;

  /* initialize arrays */
  for (i = 0; i < MFIELD_IDX_END; ++i)
    nmodel[i] = 0;

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

      if (data2 == NULL)
        magdata_copy_track(&params, i, data, track_p, mdata, npts);
      else
        magdata_copy_track_EW(&params, i, data, track_p, data2, track_p2, mdata, npts);
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

      if (fitting_flags & MAGDATA_FLG_FIT_MF)
        {
          if (fabs(mdata->qdlat[i]) <= MFIELD_QDLAT_HIGH)
            {
              /* don't fit scalar data at low latitudes if vector is available */
              if (MAGDATA_ExistVector(mdata->flags[i]))
                mdata->flags[i] &= ~(MAGDATA_FLG_F | MAGDATA_FLG_DF_NS | MAGDATA_FLG_DF_EW);
            }
          else
            {
#if MFIELD_FIT_Z_HIGHLAT
              /* only fit Z data at high latitudes */
              mdata->flags[i] &= ~(MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_F);
              mdata->flags[i] &= ~(MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DF_NS);
              mdata->flags[i] &= ~(MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DF_EW);
#else
              /* only fit F data at high latitudes */
              mdata->flags[i] &= ~(MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z);
              mdata->flags[i] &= ~(MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS);
              mdata->flags[i] &= ~(MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW);
#endif
            }

#if MFIELD_FIT_Z_HIGHLAT
          /*XXX: sometimes there is a scalar measurement but no vector -
           * forcably discard all scalar measurements */
          mdata->flags[i] &= ~(MAGDATA_FLG_F | MAGDATA_FLG_DF_NS | MAGDATA_FLG_DF_EW);
#endif

          if (MAGDATA_ExistX(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_X]);

          if (MAGDATA_ExistY(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_Y]);

          if (MAGDATA_ExistZ(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_Z]);

          if (MAGDATA_ExistScalar(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_F]);

          /* count north-south gradient data */

          if (MAGDATA_ExistDX_NS(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DX_NS]);

          if (MAGDATA_ExistDY_NS(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DY_NS]);

          if (MAGDATA_ExistDZ_NS(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DZ_NS]);

          if (MAGDATA_ExistDF_NS(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DF_NS]);

          /* count east-west gradient data */

          if (MAGDATA_ExistDX_EW(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DX_EW]);

          if (MAGDATA_ExistDY_EW(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DY_EW]);

          if (MAGDATA_ExistDZ_EW(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DZ_EW]);

          if (MAGDATA_ExistDF_EW(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DF_EW]);
        }

      if ((fitting_flags & MAGDATA_FLG_FIT_EULER) &&
          MAGDATA_ExistVector(mdata->flags[i]))
        ++(nmodel[MFIELD_IDX_B_EULER]);

      mdata->flags[i] |= fitting_flags;
    }

  fprintf(stderr, "\n");
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) scalar measurements available\n",
          npts[0], mdata->n, (double) npts[0] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) vector measurements available\n",
          npts[1], mdata->n, (double) npts[1] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) along-track scalar measurements available\n",
          npts[2], mdata->n, (double) npts[2] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) along-track vector measurements available\n",
          npts[3], mdata->n, (double) npts[3] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) east-west scalar measurements available\n",
          npts[4], mdata->n, (double) npts[4] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) east-west vector measurements available\n",
          npts[5], mdata->n, (double) npts[5] / (double) mdata->n * 100.0);

  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) X vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_X], mdata->n, (double) nmodel[MFIELD_IDX_X] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) Y vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_Y], mdata->n, (double) nmodel[MFIELD_IDX_Y] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) Z vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_Z], mdata->n, (double) nmodel[MFIELD_IDX_Z] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) scalar measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_F], mdata->n, (double) nmodel[MFIELD_IDX_F] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) vector measurements selected for Euler angle modeling\n",
          nmodel[MFIELD_IDX_B_EULER], mdata->n, (double) nmodel[MFIELD_IDX_B_EULER] / (double) mdata->n * 100.0);

  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) North/South dX vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DX_NS], mdata->n, (double) nmodel[MFIELD_IDX_DX_NS] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) North/South dY vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DY_NS], mdata->n, (double) nmodel[MFIELD_IDX_DY_NS] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) North/South dZ vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DZ_NS], mdata->n, (double) nmodel[MFIELD_IDX_DZ_NS] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) North/South dF scalar measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DF_NS], mdata->n, (double) nmodel[MFIELD_IDX_DF_NS] / (double) mdata->n * 100.0);

  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) East/West dX vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DX_EW], mdata->n, (double) nmodel[MFIELD_IDX_DX_EW] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) East/West dY vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DY_EW], mdata->n, (double) nmodel[MFIELD_IDX_DY_EW] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) East/West dZ vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DZ_EW], mdata->n, (double) nmodel[MFIELD_IDX_DZ_EW] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) East/West dF scalar measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DF_EW], mdata->n, (double) nmodel[MFIELD_IDX_DF_EW] / (double) mdata->n * 100.0);

  return mdata;
}

satdata_mag *
read_swarm(const char *filename, const int asmv)
{
  size_t nflag;
  satdata_mag *data;
  struct timeval tv0, tv1;

  fprintf(stderr, "read_swarm: reading %s...", filename);
  gettimeofday(&tv0, NULL);

  if (asmv)
    data = satdata_swarm_asmv_read_idx(filename, 0);
  else
    data = satdata_swarm_read_idx(filename, 0);

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

  fprintf(stderr, "read_champ: reading %s...", filename);
  gettimeofday(&tv0, NULL);
  data = satdata_champ_read_idx(filename, 0);
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
check_lt()
  Check if a given LT is within [lt_min,lt_max] accounting for mod 24. So it
is possible to have input lt_min < lt_max in order to select data across midnight.

Example: [lt_min,lt_max] = [6,18] will select daytime data between 6am and 6pm
         [lt_min,lt_max] = [18,6] will select nighttime data between 6pm and 6am
         [lt_min,lt_max] = [22,5] will select nighttime data between 10pm and 5am
         [lt_min,lt_max] = [0,5] will select data between midnight and 5am
*/

static int
check_lt(const double lt, const double lt_min, const double lt_max)
{
  double a, b;

  b = fmod(lt_max - lt_min, 24.0);
  if (b < 0.0)
    b += 24.0;

  a = fmod(lt - lt_min, 24.0);
  if (a < 0.0)
    a += 24.0;

  if (a > b)
    return 0; /* invalid local time */

  /* valid local time */
  return 1;
}

/*
model_flags()
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
  int status;

  /* check if we should fit Euler angles to this data point */
  status = check_lt(lt, MFIELD_EULER_LT_MIN, MFIELD_EULER_LT_MAX);
  if ((status == 1) &&
      (magdata_flags & MAGDATA_GLOBFLG_EULER) &&
      (fabs(qdlat) <= MFIELD_EULER_QDLAT))
    {
      flags |= MAGDATA_FLG_FIT_EULER;
    }

  /* check if we should fit main field model to this data point */
  if (fabs(lat_deg) <= MFIELD_HIGH_LATITUDE)
    {
      status = check_lt(lt, MFIELD_LT_MIN, MFIELD_LT_MAX);

      if (status == 1)
        flags |= MAGDATA_FLG_FIT_MF;
    }
  else
    {
#if !MFIELD_FILTER_SMDL
      double lat_rad = lat_deg * M_PI / 180.0;
      double zenith;

      solarpos_calc_zenith(unix_time, lat_rad, phi, &zenith, solarpos_workspace_p);
      zenith *= 180.0 / M_PI;
      assert(zenith >= 0.0);

      /* large zenith angle means darkness */
      if (zenith >= MFIELD_MAX_ZENITH)
        flags |= MAGDATA_FLG_FIT_MF;
#else
      /* if using SMDL filtering, don't use zenith criteria, use all high
       * latitude points which pass the SMDL filter */
      flags |= MAGDATA_FLG_FIT_MF;
#endif
    }

  return flags;
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

#if 0
  print_unflagged_data("data_ts/data_ts.0", data);
#endif

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
#if 0
      double thresh[] = { 20.0, 25.0, 15.0, 15.0 };
#else
      /* 2 Jan 2017: since the main field model used for rms is out of date,
       * need to increase thresholds to prevent throwing out good data
       */
      double thresh[] = { 30.0, 30.0, 30.0, 30.0 };
#endif

      nrms = track_flag_rms(rmsfile, thresh, NULL, data, track_p);
      fprintf(stderr, "preprocess_data: flagged (%zu/%zu) (%.1f%%) tracks due to high rms\n",
              nrms, track_p->n, (double) nrms / (double) track_p->n * 100.0);
    }

#if 0
  print_unflagged_data("data_ts/data_ts.1", data);
#endif

#if MFIELD_FILTER_SMDL

  /* flag according to SMDL index */
  {
    size_t nsmdl;
    const double smdl_lo = 6.0;
    const double smdl_hi = 50.0;

    fprintf(stderr, "main: flagging according to SMDL thresholds (lo = %g, hi = %g)...",
            smdl_lo, smdl_hi);
    nsmdl = satdata_filter_smdl(smdl_lo, smdl_hi, data);
    fprintf(stderr, "done (%zu/%zu (%.1f%%) points flagged\n",
            nsmdl, data->n, (double)nsmdl / (double)data->n * 100.0);
  }

#else

  /* flag according to WMM criteria */
  satdata_filter_wmm(1, data);

#endif

#if 0
  print_unflagged_data("data_ts/data_ts.2", data);
#endif

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

#if 0
  print_unflagged_data("data_ts/data_ts.3", data);
#endif

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

#if 0
  print_unflagged_data("data_ts/data_ts.4", data);
#endif

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

/* print unflagged data points for data selection visualization */
void
print_unflagged_data(const char *filename, const satdata_mag *data)
{
  size_t i;
  FILE *fp;
  const double tmin = satdata_epoch2year(data->t[0]);
  const double tmax = satdata_epoch2year(data->t[data->n - 1]);
  const double dt = 10.0;            /* bin size in days */
  const double dt_yrs = dt / 365.25; /* convert to years */
  const size_t n = (size_t) ((tmax - tmin) / dt_yrs);
  gsl_histogram *h = gsl_histogram_alloc(n);

  gsl_histogram_set_ranges_uniform(h, tmin, tmax);

  fp = fopen(filename, "w");

  for (i = 0; i < data->n; ++i)
    {
      if (data->flags[i])
        continue;

      gsl_histogram_increment(h, satdata_epoch2year(data->t[i]));
    }

  gsl_histogram_fprintf(fp, h, "%g", "%g");

  fclose(fp);
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --champ_file      | -c champ_index_file       - CHAMP index file\n");
  fprintf(stderr, "\t --swarm_file      | -s swarm_index_file       - Swarm index file\n");
  fprintf(stderr, "\t --swarm_file2     | -t swarm_index_file2      - Swarm index file 2 (for E/W gradients)\n");
  fprintf(stderr, "\t --swarm_asmv_file | -a swarm_asmv_index_file  - Swarm ASM-V index file\n");
  fprintf(stderr, "\t --downsample      | -d downsample             - downsampling factor\n");
  fprintf(stderr, "\t --euler_file      | -e euler_file             - Euler angles file\n");
  fprintf(stderr, "\t --output_file     | -o output_file            - binary output data file (magdata format)\n");
}

int
main(int argc, char *argv[])
{
  char *datamap_file = "datamap.dat";
  char *data_file = "data.dat";
  char *output_file = NULL;
  satdata_mag *data = NULL;
  satdata_mag *data2 = NULL;
  magdata *mdata;
  euler_workspace *euler_p = NULL;
  struct timeval tv0, tv1;
  track_workspace *track_p = NULL;
  track_workspace *track_p2 = NULL;
  preprocess_parameters params;
  size_t magdata_flags = 0;
  size_t magdata_flags2 = 0;

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
          { "swarm_file2", required_argument, NULL, 't' },
          { "swarm_asmv_file", required_argument, NULL, 'a' },
          { "champ_file", required_argument, NULL, 'c' },
          { "downsample", required_argument, NULL, 'd' },
          { "output_file", required_argument, NULL, 'o' },
          { "euler_file", required_argument, NULL, 'e' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:c:d:e:o:s:t:", long_options, &option_index);
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

          /* For E/W gradients */
          case 't':
            data2 = read_swarm(optarg, 0);
            magdata_flags2 = MAGDATA_GLOBFLG_EULER;
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

  fprintf(stderr, "main: === PREPROCESSING SATELLITE 1 ===\n");
  track_p = preprocess_data(&params, magdata_flags, data);

  if (data2)
    {
      fprintf(stderr, "main: === PREPROCESSING SATELLITE 2 ===\n");
      track_p2 = preprocess_data(&params, magdata_flags2, data2);
    }

  fprintf(stderr, "main: computing vector residuals...");
  gettimeofday(&tv0, NULL);
  mdata = copy_data(magdata_flags, data, track_p, magdata_flags2, data2, track_p2);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu residuals copied, %g seconds)\n",
          mdata->n, time_diff(tv0, tv1));

  /* free data after copying arrays to free up memory */
  satdata_mag_free(data);

  if (data2)
    satdata_mag_free(data2);

  magdata_init(mdata);

  fprintf(stderr, "main: computing spatial weighting of data...");
  gettimeofday(&tv0, NULL);
  magdata_calc(mdata);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

#if 0
  fprintf(stderr, "main: writing data to %s...", data_file);
  magdata_print(data_file, mdata);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing data map to %s...", datamap_file);
  magdata_map(datamap_file, mdata);
  fprintf(stderr, "done\n");
#endif

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

  if (track_p2)
    track_free(track_p2);

  return 0;
} /* main() */
