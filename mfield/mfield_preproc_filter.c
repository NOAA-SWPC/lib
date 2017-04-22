/*
 * mfield_preproc_filter.c
 *
 * This module contains routines for filtering data according
 * to various indices requirements
 */

static size_t mfield_flag_RC(const double RC_max, track_workspace *track_p, satdata_mag *data);
static size_t mfield_flag_dRC(const double dRC_max, track_workspace *track_p, satdata_mag *data);
static int mfield_preprocess_filter(const preprocess_parameters *params, track_workspace *track_p, satdata_mag *data);

/*
mfield_flag_RC()
  Flag tracks for RC criteria

Inputs: RC_max  - maximum allowed RC index (nT)
        track_p - track workspace
        data    - data

Return: number of tracks flagged
*/

static size_t
mfield_flag_RC(const double RC_max, track_workspace *track_p, satdata_mag *data)
{
  size_t ntrack_flagged = 0; /* number of tracks flagged */
  rc_workspace *rc_p;
  size_t i;
  int s;

  if (data->n == 0)
    return 0;

  rc_p = rc_alloc(RC_IDX_FILE);

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;
      time_t t1 = satdata_epoch2timet(data->t[start_idx]);
      time_t t2 = satdata_epoch2timet(tptr->t_eq);
      time_t t3 = satdata_epoch2timet(data->t[end_idx]);
      double RC1, RC2, RC3;

      s = rc_get_RC(t1, &RC1, rc_p);
      s += rc_get_RC(t2, &RC2, rc_p);
      s += rc_get_RC(t3, &RC3, rc_p);
      if (s)
        {
          fprintf(stderr, "preprocess_flag_RC: error: RC not available for track %zu\n", i);
          continue;
        }

      if ((fabs(RC1) > RC_max) || (fabs(RC2) > RC_max) || (fabs(RC3) > RC_max))
        {
          track_flag_track(i, TRACK_FLG_RC, data, track_p);
          ++ntrack_flagged;
        }
    }

  rc_free(rc_p);

  return ntrack_flagged;
}

/*
mfield_flag_dRC()
  Flag tracks for dRC/dt criteria

Reject track if:

1. |dRC/dt| > dRC_max at any time along the track
2. |dRC/dt| > dRC_max within 3 hours prior to track

Inputs: dRC_max - maximum allowed dRC/dt index (nT/hour)
        track_p - track workspace
        data    - data

Return: number of tracks flagged
*/

static size_t
mfield_flag_dRC(const double dRC_max, track_workspace *track_p, satdata_mag *data)
{
  size_t ntrack_flagged = 0; /* number of tracks flagged */
  rc_workspace *rc_p;
  size_t i;
  int s;

  if (data->n == 0)
    return 0;

  rc_p = rc_alloc(RC_IDX_FILE);

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;
      time_t t1 = satdata_epoch2timet(data->t[start_idx]);
      time_t t2 = satdata_epoch2timet(tptr->t_eq);
      time_t t3 = satdata_epoch2timet(data->t[end_idx]);
      double dRC1, dRC2, dRC3;
      int flag_track = 0;

      s = rc_deriv_get(t1, &dRC1, rc_p);
      s += rc_deriv_get(t2, &dRC2, rc_p);
      s += rc_deriv_get(t3, &dRC3, rc_p);
      if (s)
        {
          fprintf(stderr, "preprocess_flag_RC: error: dRC/dt not available for track %zu\n", i);
          flag_track = 1;
        }

      if ((fabs(dRC1) > dRC_max) || (fabs(dRC2) > dRC_max) || (fabs(dRC3) > dRC_max))
        {
          flag_track = 1;
        }

      /* now check |dRC/dt| for 3 hours prior to track */
      s = rc_deriv_get(t1 - 3600, &dRC1, rc_p);
      s += rc_deriv_get(t1 - 2*3600, &dRC2, rc_p);
      s += rc_deriv_get(t1 - 3*3600, &dRC3, rc_p);
      if (s)
        {
          fprintf(stderr, "preprocess_flag_RC: error: dRC/dt not available for track %zu\n", i);
          flag_track = 1;
        }

      if ((fabs(dRC1) > dRC_max) || (fabs(dRC2) > dRC_max) || (fabs(dRC3) > dRC_max))
        {
          flag_track = 1;
        }

      if (flag_track)
        {
          track_flag_track(i, TRACK_FLG_RC, data, track_p);
          ++ntrack_flagged;
        }
    }

  rc_free(rc_p);

  return ntrack_flagged;
}

/*
mfield_preprocess_filter()
  Select data for geomagnetically quiet periods
*/

static int
mfield_preprocess_filter(const preprocess_parameters *params, track_workspace *track_p, satdata_mag *data)
{
  int s = 0;
  size_t nflagged_kp,
         nflagged_RC,
         nflagged_dRC;

  nflagged_kp = track_flag_kp(0.0, params->max_kp, data, track_p);
  /*nflagged_RC = mfield_flag_RC(RC_max, track_p, data);*/
  nflagged_dRC = mfield_flag_dRC(params->max_dRC, track_p, data);

  fprintf(stderr, "\n");

  fprintf(stderr, "\t mfield_preprocess_filter: flagged %zu/%zu (%.1f%%) tracks due to kp [%.1f]\n",
          nflagged_kp, track_p->n, (double) nflagged_kp / (double) track_p->n * 100.0, params->max_kp);
  /*fprintf(stderr, "\t mfield_preprocess_filter: flagged %zu/%zu (%.1f%%) tracks due to RC\n",
          nflagged_RC, track_p->n, (double) nflagged_RC / (double) track_p->n * 100.0);*/
  fprintf(stderr, "\t mfield_preprocess_filter: flagged %zu/%zu (%.1f%%) tracks due to dRC/dt [%.1f nT/hour]\n",
          nflagged_dRC, track_p->n, (double) nflagged_dRC / (double) track_p->n * 100.0, params->max_dRC);

  return s;
}
