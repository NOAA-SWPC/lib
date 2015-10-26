/*
 * track_flag.c
 * Track filtering for rms, local time, etc
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <satdata/satdata.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>

#include "common.h"

#include "track.h"

static int track_calc_rms(const double qdlat_min, const double qdlat_max,
                          const track_data *tptr,
                          size_t *nrms_scal, size_t *nrms_vec, double rms[4], satdata_mag *data);
static int track_search_jump(const double threshold, const double qdlat_min,
                             const double qdlat_max, const track_data *tptr,
                             satdata_mag *data);
static size_t track_flag_outliers(const size_t sidx, const size_t eidx,
                                  const double thresh_f, const double thresh_z,
                                  satdata_mag *data);

/*
track_flag_rms()
  Compute rms residuals of each satellite track for latitudes in
[-QDLAT_MAX,QDLAT_MAX], and flag tracks whose rms is higher than
a given threshold

Inputs: outfile  - file to write rms info
        thresh[] - rms thresholds (nT)
                   thresh[0] = X rms threshold (nT)
                   thresh[1] = Y rms threshold (nT)
                   thresh[2] = Z rms threshold (nT)
                   thresh[3] = F rms threshold (nT)
        data     - satellite data
        w        - track workspace

Return: number of data points flagged due to high rms

Notes:
1) Tracks with high rms are flagged with SATDATA_FLG_OUTLIER

2) To disable a certain rms component test, set thresh[i] to -1.0

3) On output, the 'rms' and 'nrms' fields of each track are filled in
*/

size_t
track_flag_rms(const char *outfile, const double thresh[4],
               satdata_mag *data, track_workspace *w)
{
  FILE *fp;
  size_t i, j;
  size_t nflagged = 0;        /* number of points flagged */
  size_t ntrack_flagged = 0;  /* number of entire tracks flagged for rms */

  if (data->n == 0)
    return 0;

  fp = fopen(outfile, "w");

  i = 1;
  fprintf(fp, "# Field %zu: timestamp of equator crossing (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
  fprintf(fp, "# Field %zu: time of equator crossing (years)\n", i++);
  fprintf(fp, "# Field %zu: longitude of equator crossing (degrees)\n", i++);
  fprintf(fp, "# Field %zu: local time of equator crossing (hours)\n", i++);
  fprintf(fp, "# Field %zu: track mean altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: number of data points in track\n", i++);
  fprintf(fp, "# Field %zu: track X rms (nT)\n", i++);
  fprintf(fp, "# Field %zu: track Y rms (nT)\n", i++);
  fprintf(fp, "# Field %zu: track Z rms (nT)\n", i++);
  fprintf(fp, "# Field %zu: track scalar rms (nT)\n", i++);
  fprintf(fp, "# Field %zu: satellite direction (+1 north -1 south)\n", i++);
  fprintf(fp, "# Field %zu: number of along-track data used for scalar rms\n", i++);
  fprintf(fp, "# Field %zu: number of along-track data used for vector rms\n", i++);
  fprintf(fp, "# Field %zu: flagged due to rms (1 or 0)\n", i++);

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      size_t sidx = tptr->start_idx;
      size_t eidx = tptr->end_idx;
      int s;
      int flag_rms = 0;

      /* compute rms of track residuals for low-latitudes for |qdlat| < QDLAT_MAX */
      track_calc_rms(0.0, TRACK_QDLAT_MAX, tptr, &(tptr->nrms_scal), &(tptr->nrms_vec), tptr->rms, data);

      /* compute rms of track residuals for polar latitudes for |qdlat| > QDLAT_MAX */
      track_calc_rms(TRACK_QDLAT_MAX, 100.0, tptr, &(tptr->nrmspol_scal), &(tptr->nrmspol_vec), tptr->rmspol, data);

      /* check if there are enough points for good rms estimate */
      s = 0;
      if (thresh[3] > 0.0 && tptr->nrms_scal < TRACK_RMS_NDATA)
        s = -1;

#if 0
      /* don't check vector nrms since some CHAMP local-times have 0 vector
       * data but still have good scalar data
       */
      if ((thresh[0] > 0.0 || thresh[1] > 0.0 || thresh[2] > 0.0) &&
          tptr->nrms_vec < TRACK_RMS_NDATA)
        s = -1;
#endif

      if (s != 0)
        {
          /* too few points for rms calculation */
#if TRACK_DEBUG
          fprintf(stderr, "track_flag_rms: failed rms, flagging pts [%zu, %zu]\n", sidx, eidx);
#endif
          flag_rms = 1;
        }
      else
        {
          /* check rms against thresholds */
          for (j = 0; j < 4; ++j)
            {
              if (thresh[j] > 0.0 && tptr->rms[j] > thresh[j])
                flag_rms = 1;
            }

          if (!flag_rms)
            {
              /*
               * otherwise, look for individual outliers which weren't part of
               * the rms computation (ie: at the poles)
               */
              nflagged += track_flag_outliers(sidx, eidx, 10.0*thresh[3], 10.0*thresh[2], data);
            }
        }

      if (flag_rms)
        {
          /* failed rms test */
          nflagged += track_flag_track(i, TRACK_FLG_RMS, data, w);
          ++ntrack_flagged;
        }

      fprintf(fp, "%ld %f %10.4f %6.2f %6.2f %6zu %10.4f %10.4f %10.4f %10.4f %2d %6zu %6zu %d\n",
              satdata_epoch2timet(tptr->t_eq),
              satdata_epoch2year(tptr->t_eq),
              wrap180(tptr->lon_eq),
              tptr->lt_eq,
              tptr->meanalt,
              tptr->n,
              tptr->rms[0],
              tptr->rms[1],
              tptr->rms[2],
              tptr->rms[3],
              satdata_mag_satdir(sidx, data),
              tptr->nrms_scal,
              tptr->nrms_vec,
              flag_rms);
      fflush(fp);
    }

  fclose(fp);

  fprintf(stderr, "track_flag_rms: flagged %zu/%zu (%.1f%%) tracks due to rms\n",
          ntrack_flagged, w->n, (double) ntrack_flagged / (double) w->n * 100.0);
  fprintf(stderr, "track_flag_rms: rms data written to %s\n", outfile);

  return nflagged;
} /* track_flag_rms() */

/*
track_flag_jumps()
  Check for large jumps in vector residuals between adjacent data points

Inputs: thresh   - jump threshold (nT)
        data     - satellite data
        w        - track workspace

Return: number of data points flagged due to jumps

Notes:
1) Tracks with jumps are flagged with SATDATA_FLG_OUTLIER
*/

size_t
track_flag_jumps(const double thresh, satdata_mag *data, track_workspace *w)
{
  size_t i;
  size_t nflagged = 0;        /* number of points flagged */
  size_t ntrack_flagged = 0;  /* number of entire tracks flagged for jumps */

  if (data->n == 0)
    return 0;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      int s;

      s = track_search_jump(thresh, 0.0, TRACK_QDLAT_MAX, tptr, data);
      if (s)
        {
          /* failed jump test */
          nflagged += track_flag_track(i, TRACK_FLG_JUMP, data, w);
          ++ntrack_flagged;
        }
    }

  fprintf(stderr, "track_flag_jumps: flagged %zu/%zu (%.1f%%) tracks due to jumps\n",
          ntrack_flagged, w->n, (double) ntrack_flagged / (double) w->n * 100.0);

  return nflagged;
} /* track_flag_jumps() */

/*
track_flag_satdir()
  Flag any tracks flying in a specified direction (north/south)

Inputs: satdir - +1: flag north-flying tracks
                 -1: flag south-flying tracks
        data   - satellite data
        w      - track workspace

Return: number of data flagged
*/

size_t
track_flag_satdir(const int satdir, satdata_mag *data, track_workspace *w)
{
  size_t i;
  size_t nflagged = 0;        /* number of points flagged */
  size_t ntrack_flagged = 0;  /* number of entire tracks flagged for satdir */

  if (data->n == 0)
    return 0;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);

      if (tptr->satdir == satdir)
        {
          nflagged += track_flag_track(i, TRACK_FLG_SATDIR, data, w);
          ++ntrack_flagged;
        }
    }

  fprintf(stderr, "track_flag_satdir: flagged %zu/%zu (%.1f%%) tracks due to satellite direction\n",
          ntrack_flagged, w->n, (double) ntrack_flagged / (double) w->n * 100.0);

  return nflagged;
} /* track_flag_satdir() */

/*
track_flag_ut()
  Flag any tracks outside of [ut_min,ut_max]. The UT for
comparison is the UT at the time of the equator crossing

Inputs: ut_min - minimum UT (hours)
        ut_max - maximum UT (hours)
        data   - satellite data
        w      - track workspace

Return: number of data flagged
*/

size_t
track_flag_ut(const double ut_min, const double ut_max, satdata_mag *data, track_workspace *w)
{
  size_t i;
  size_t nflagged = 0;        /* number of points flagged */
  size_t ntrack_flagged = 0;  /* number of entire tracks flagged for UT */

  if (data->n == 0)
    return 0;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      double ut = satdata_epoch2ut(tptr->t_eq);

      if (ut < ut_min || ut > ut_max)
        {
          nflagged += track_flag_track(i, TRACK_FLG_UT, data, w);
          ++ntrack_flagged;
        }
    }

  fprintf(stderr, "track_flag_ut: flagged %zu/%zu (%.1f%%) tracks due to UT\n",
          ntrack_flagged, w->n, (double) ntrack_flagged / (double) w->n * 100.0);

  return nflagged;
} /* track_flag_ut() */

/*
track_flag_lt()
  Flag any tracks outside of [lt_min,lt_max]. The LT for
comparison is the LT at the time of the equator crossing

Inputs: lt_min - minimum LT (hours)
        lt_max - maximum LT (hours)
        data   - satellite data
        w      - track workspace

Return: number of data flagged
*/

size_t
track_flag_lt(const double lt_min, const double lt_max, satdata_mag *data, track_workspace *w)
{
  size_t i;
  size_t nflagged = 0;        /* number of points flagged */
  size_t ntrack_flagged = 0;  /* number of entire tracks flagged for UT */
  double a, b;

  if (data->n == 0)
    return 0;

  b = fmod(lt_max - lt_min, 24.0);
  if (b < 0.0)
    b += 24.0;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      time_t t = satdata_epoch2timet(tptr->t_eq);
      double lt = get_localtime(t, tptr->lon_eq * M_PI / 180.0);

      a = fmod(lt - lt_min, 24.0);
      if (a < 0.0)
        a += 24.0;

      if (a > b)
        {
          nflagged += track_flag_track(i, TRACK_FLG_LT, data, w);
          ++ntrack_flagged;
        }
    }

  fprintf(stderr, "track_flag_lt: flagged %zu/%zu (%.1f%%) tracks due to LT\n",
          ntrack_flagged, w->n, (double) ntrack_flagged / (double) w->n * 100.0);

  return nflagged;
} /* track_flag_lt() */

/*
track_flag_season()
  Flag any tracks outside of [lt_min,lt_max]. The UT for
comparison is the UT at the time of the equator crossing

Inputs: callback - callback function to check whether a given
                   doy argument is acceptable; return values:
                   -1 - reject doy
                   0  - accept doy
        data     - satellite data
        w        - track workspace

Return: number of data flagged
*/

size_t
track_flag_season(int (*callback)(const double doy, const void *params),
                  const void *params, satdata_mag *data, track_workspace *w)
{
  size_t i;
  size_t nflagged = 0;        /* number of points flagged */
  size_t ntrack_flagged = 0;  /* number of entire tracks flagged for UT */

  if (data->n == 0)
    return 0;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      time_t t = satdata_epoch2timet(tptr->t_eq);
      double doy = get_season(t);
      int s = (*callback)(doy, params);

      if (s)
        {
          nflagged += track_flag_track(i, TRACK_FLG_DOY, data, w);
          ++ntrack_flagged;
        }
    }

  fprintf(stderr, "track_flag_season: flagged %zu/%zu (%.1f%%) tracks due to season\n",
          ntrack_flagged, w->n, (double) ntrack_flagged / (double) w->n * 100.0);

  return nflagged;
} /* track_flag_season() */

/*
track_flag_lon()
  Flag any tracks outside of [lon_min,lon_max] using
the longitude of the equator crossing

Inputs: lon_min  - minimum longitude in [-180,180] (degrees)
        lon_max  - maximum longitude in [-180,180] (degrees)
        data     - satellite data
        w        - track workspace

Return: number of data flagged
*/

size_t
track_flag_lon(const double lon_min, const double lon_max,
               satdata_mag *data, track_workspace *w)
{
  size_t i;
  size_t nflagged = 0;        /* number of points flagged */
  size_t ntrack_flagged = 0;  /* number of entire tracks flagged for UT */

  if (data->n == 0)
    return 0;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      double lon = wrap180(tptr->lon_eq);

      if (lon < lon_min || lon > lon_max)
        {
          nflagged += track_flag_track(i, TRACK_FLG_LONGITUDE, data, w);
          ++ntrack_flagged;
        }
    }

  fprintf(stderr, "track_flag_lon: flagged %zu/%zu (%.1f%%) tracks due to longitude\n",
          ntrack_flagged, w->n, (double) ntrack_flagged / (double) w->n * 100.0);

  return nflagged;
} /* track_flag_lon() */

/*
track_flag_kp()
  Flag any tracks with kp outside of [kp_min,kp_max]. The kp for
comparison is the kp at the time of the equator crossing. Since
kp is a 3-hour index, it shouldn't change over the entire track.

Inputs: kp_min - minimum kp
        kp_max - maximum kp
        data   - satellite data
        w      - track workspace

Return: number of data flagged
*/

size_t
track_flag_kp(const double kp_min, const double kp_max, satdata_mag *data, track_workspace *w)
{
  size_t i;
  size_t nflagged = 0;        /* number of points flagged */
  size_t ntrack_flagged = 0;  /* number of entire tracks flagged for UT */
  kp_workspace *kp_p;
  int s;

  if (data->n == 0)
    return 0;

  kp_p = kp_alloc(KP_IDX_FILE);

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      time_t t = satdata_epoch2timet(tptr->t_eq);
      double kp;

      s = kp_get(t, &kp, kp_p);
      if (s == 0 &&
          (kp < kp_min || kp > kp_max))
        {
          nflagged += track_flag_track(i, TRACK_FLG_KP, data, w);
          ++ntrack_flagged;
        }
    }

  fprintf(stderr, "track_flag_kp: flagged %zu/%zu (%.1f%%) tracks due to kp\n",
          ntrack_flagged, w->n, (double) ntrack_flagged / (double) w->n * 100.0);

  kp_free(kp_p);

  return nflagged;
} /* track_flag_kp() */

/*
track_flag_meanalt()
  Flag any tracks with mean altitude outside of [alt_min,alt_max].

Inputs: alt_min - minimum altitude (km)
        alt_max - maximum altitude (km)
        data    - satellite data
        w       - track workspace

Return: number of data flagged
*/

size_t
track_flag_meanalt(const double alt_min, const double alt_max, satdata_mag *data, track_workspace *w)
{
  size_t i;
  size_t nflagged = 0;        /* number of points flagged */
  size_t ntrack_flagged = 0;  /* number of entire tracks flagged for altitude */

  if (data->n == 0)
    return 0;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);

      if (tptr->meanalt < alt_min || tptr->meanalt > alt_max)
        {
          nflagged += track_flag_track(i, TRACK_FLG_ALTITUDE, data, w);
          ++ntrack_flagged;
        }
    }

  fprintf(stderr, "track_flag_meanalt: flagged %zu/%zu (%.1f%%) tracks due to altitude\n",
          ntrack_flagged, w->n, (double) ntrack_flagged / (double) w->n * 100.0);

  return nflagged;
} /* track_flag_meanalt() */

/*
track_flag_incomplete()
  Flag incomplete tracks with missing data. Tracks which do not
contain data in [qd_min,qd_max] are discarded.

Inputs: qd_min - minimum QD latitude (degrees)
        qd_max - maximum QD latitude (degrees)
        data   - satellite data
        w      - track workspace

Return: number of data flagged
*/

size_t
track_flag_incomplete(const double qd_min, const double qd_max, satdata_mag *data, track_workspace *w)
{
  size_t i;
  size_t nflagged = 0;        /* number of points flagged */
  size_t ntrack_flagged = 0;  /* number of entire tracks flagged for missing data */

  if (data->n == 0)
    return 0;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      double qd0 = GSL_MIN(data->qdlat[tptr->start_idx], data->qdlat[tptr->end_idx]);
      double qd1 = GSL_MAX(data->qdlat[tptr->start_idx], data->qdlat[tptr->end_idx]);

      if (qd0 > qd_min || qd1 < qd_max)
        {
          nflagged += track_flag_track(i, TRACK_FLG_INCOMPLETE, data, w);
          ++ntrack_flagged;
        }
    }

  fprintf(stderr, "track_flag_incomplete: flagged %zu/%zu (%.1f%%) tracks due to missing data\n",
          ntrack_flagged, w->n, (double) ntrack_flagged / (double) w->n * 100.0);

  return nflagged;
} /* track_flag_incomplete() */

/*
track_flag_n()
  Flag tracks with less than a given number of unflagged data
points. Only currently unflagged tracks are counted

Inputs: nmin - minimum number of unflagged data points
        data - satellite data
        w    - track workspace

Return: number of data flagged
*/

size_t
track_flag_n(const size_t nmin, satdata_mag *data, track_workspace *w)
{
  size_t i;
  size_t nflagged = 0;        /* number of points flagged */
  size_t ntrack_flagged = 0;  /* number of entire tracks flagged for missing data */

  if (data->n == 0)
    return 0;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      size_t ngood;

      /* don't count already flagged tracks */
      if (tptr->flags)
        continue;
      
      ngood = tptr->n - track_data_nflagged(tptr, data);

      if (ngood < nmin)
        {
          nflagged += track_flag_track(i, TRACK_FLG_INCOMPLETE, data, w);
          ++ntrack_flagged;
        }
    }

  fprintf(stderr, "track_flag_n: flagged %zu/%zu (%.1f%%) tracks due to missing data\n",
          ntrack_flagged, w->n, (double) ntrack_flagged / (double) w->n * 100.0);

  return nflagged;
} /* track_flag_n() */

/*****************************************************
 * INTERNAL ROUTINES                                 *
 *****************************************************/

/*
track_calc_rms()
  Compute rms of vector and scalar residuals along satellite track for
QD latitudes in a given range, ignoring previously flagged data

Inputs: qdlat_min - points which satisfy the following are added to rms
                    qdlat_min <= |qdlat| <= qdlat_max
        qdlat_max - see above
        start_idx - starting index of track
        end_idx   - ending index of track
        nrms_scal - (output) number of data points used in scalar rms calculation
        nrms_vec  - (output) number of data points used in vector rms calculation
        rms       - (output)
                    rms[0] = rms_X
                    rms[1] = rms_Y
                    rms[2] = rms_Z
                    rms[3] = rms_F
        data      - satellite data

Return: success or error
*/

static int
track_calc_rms(const double qdlat_min, const double qdlat_max,
               const track_data *tptr,
               size_t *nrms_scal, size_t *nrms_vec, double rms[4], satdata_mag *data)
{
  int s = 0;
  size_t nscal = 0;
  size_t nvec = 0;
  size_t i;

  for (i = 0; i < 4; ++i)
    rms[i] = 0.0;

  /* compute along-track rms */
  for (i = 0; i < tptr->n; ++i)
    {
      size_t didx = i + tptr->start_idx;

      /* ignore bad/missing data */
      if (SATDATA_BadData(data->flags[didx]))
        continue;

      if (fabs(data->qdlat[didx]) < qdlat_min ||
          fabs(data->qdlat[didx]) > qdlat_max)
        continue;

      if (SATDATA_ExistVector(data->flags[didx]))
        {
          rms[0] += pow(tptr->Bx[i], 2.0);
          rms[1] += pow(tptr->By[i], 2.0);
          rms[2] += pow(tptr->Bz[i], 2.0);
          ++nvec;
        }

      rms[3] += pow(tptr->Bf[i], 2.0);
      ++nscal;
    }

  *nrms_scal = nscal;
  *nrms_vec = nvec;

  if (nvec > 0)
    {
      for (i = 0; i < 3; ++i)
        rms[i] = sqrt(rms[i] / nvec);
    }

  if (nscal > 0)
    rms[3] = sqrt(rms[3] / nscal);

  return s;
} /* track_calc_rms() */

/*
track_search_jump()
  Search vector and scalar residuals for sudden jumps between adjacent
data points

Inputs: threshold - jump threshold (nT)
        qdlat_min - points which satisfy the following are checked for jumps:
                    qdlat_min <= |qdlat| <= qdlat_max
        qdlat_max - see above
        tptr      - track pointer
        data      - satellite data

Return: 0 if no jumps found, -1 if jump found
*/

static int
track_search_jump(const double threshold, const double qdlat_min,
                  const double qdlat_max, const track_data *tptr,
                  satdata_mag *data)
{
  int s = 0;
  size_t i;

  for (i = 0; i < tptr->n - 1; ++i)
    {
      size_t didx = tptr->start_idx + i;

      if (SATDATA_BadData(data->flags[didx]))
        continue;

      if (fabs(data->qdlat[didx]) < qdlat_min ||
          fabs(data->qdlat[didx]) > qdlat_max)
        continue;

      /* check vector data for jumps */
      if (SATDATA_ExistVector(data->flags[didx]))
        {
          if (fabs(tptr->Bx[i] - tptr->Bx[i + 1]) > threshold)
            return -1;
          if (fabs(tptr->By[i] - tptr->By[i + 1]) > threshold)
            return -1;
          if (fabs(tptr->Bz[i] - tptr->Bz[i + 1]) > threshold)
            return -1;
        }

      /* check scalar data for jumps */
      if (fabs(tptr->Bf[i] - tptr->Bf[i + 1]) > threshold)
        return -1;
    }

  return s;
} /* track_search_jump() */

/*
track_flag_outliers()
  Search for individual outliers and flag if found

Inputs: sidx     - starting index of track in data
        eidx     - ending index of track in data
        thresh_f - allowed threshold for F residual
        thresh_z - allowed threshold for Z residual
        data     - satellite data
*/

static size_t
track_flag_outliers(const size_t sidx, const size_t eidx,
                    const double thresh_f, const double thresh_z,
                    satdata_mag *data)
{
  size_t i;
  size_t n = 0;

  for (i = sidx; i <= eidx; ++i)
    {
      double resf, resz;

      /* don't count if already flagged */
      if (data->flags[i] & SATDATA_FLG_OUTLIER)
        continue;

      resz = SATDATA_VEC_Z(data->B, i) -
             SATDATA_VEC_Z(data->B_main, i);
      resf = data->F[i] - data->F_main[i];

      if ((thresh_f > 0.0 && fabs(resf) > thresh_f) ||
          (thresh_z > 0.0 && fabs(resz) > thresh_z))
        {
          data->flags[i] |= SATDATA_FLG_OUTLIER;
          ++n;
        }
    }

  return n;
} /* track_flag_outliers() */
