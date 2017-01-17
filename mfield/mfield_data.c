/*
 * mfield_data.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#include <satdata/satdata.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>

#include "mfield_data.h"

#include "common.h"

/*
mfield_data_alloc()
  Allocate mfield_data_workspace

Inputs: nsources - number of data sources (satellites,observatories)
        epoch    - model epoch (decimal years)

Return: pointer to workspace
*/

mfield_data_workspace *
mfield_data_alloc(const size_t nsources, const double epoch)
{
  mfield_data_workspace *w;

  w = calloc(1, sizeof(mfield_data_workspace));
  if (!w)
    return 0;

  w->nsources = nsources;
  w->epoch = epoch;

  w->t0 = malloc(nsources * sizeof(double));
  w->t1 = malloc(nsources * sizeof(double));
  if (!w->t0 || !w->t1)
    {
      mfield_data_free(w);
      return 0;
    }

  w->mdata = calloc(1, nsources * sizeof(magdata *));
  if (!w->mdata)
    {
      mfield_data_free(w);
      return 0;
    }

  w->rstat_workspace_p = gsl_rstat_alloc();
  if (!w->rstat_workspace_p)
    {
      mfield_data_free(w);
      return 0;
    }

  w->t_mu = -1.0;
  w->t_sigma = -1.0;
  w->t0_data = -1.0;
  w->t1_data = -1.0;

  return w;
} /* mfield_data_alloc() */

void
mfield_data_free(mfield_data_workspace *w)
{
  if (w->t0)
    free(w->t0);

  if (w->t1)
    free(w->t1);

  if (w->rstat_workspace_p)
    gsl_rstat_free(w->rstat_workspace_p);

  if (w->mdata)
    {
      size_t i;

      for (i = 0; i < w->nsources; ++i)
        {
          magdata *mdata = w->mdata[i];

          if (mdata)
            magdata_free(mdata);
        }

      free(w->mdata);
    }

  free(w);
}

/*
mfield_data_filter_time()
  Flag any data points outside of [tmin,tmax] with
MAGDATA_FLG_DISCARD

Inputs: tmin - minimum time (decimal year)
        tmax - maximum time (decimal year)
        w    - workspace

Return: number of data flagged

Notes:
1) tmin/tmax can be set to -1 to exclude them from the test
*/

size_t
mfield_data_filter_time(const double tmin, const double tmax,
                        mfield_data_workspace *w)
{
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
          double t = satdata_epoch2year(mptr->t[j]);

          if ((tmin > 0.0 && t < tmin) ||
              (tmax > 0.0 && t > tmax))
            {
              mptr->flags[j] |= MAGDATA_FLG_DISCARD;
              ++cnt;
            }
        }
    }

  return cnt;
} /* mfield_data_filter_time() */

/*
mfield_data_filter_euler()
  We are not fitting Euler angles (MFIELD_FIT_EULER is 0),
so discard any data which is marked Euler only

Inputs: w    - workspace

Return: number of data flagged
*/

size_t
mfield_data_filter_euler(mfield_data_workspace *w)
{
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
          if ((mptr->flags[j] & MAGDATA_FLG_FIT_EULER) &&
              !(mptr->flags[j] & MAGDATA_FLG_FIT_MF))
            {
              mptr->flags[j] |= MAGDATA_FLG_DISCARD;
              ++cnt;
            }
        }
    }

  return cnt;
} /* mfield_data_filter_euler() */

/*
mfield_data_filter_comp()
  Discard any data according to the MFIELD_FIT_xxx flags

Inputs: w    - workspace

Return: number of data flagged
*/

size_t
mfield_data_filter_comp(mfield_data_workspace *w)
{
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
#if !MFIELD_FIT_X
          mptr->flags[j] &= ~MAGDATA_FLG_X;
#endif

#if !MFIELD_FIT_Y
          mptr->flags[j] &= ~MAGDATA_FLG_Y;
#endif

#if !MFIELD_FIT_Z
          mptr->flags[j] &= ~MAGDATA_FLG_Z;
#endif

#if !MFIELD_FIT_DX_NS
          mptr->flags[j] &= ~MAGDATA_FLG_DX_NS;
#endif

#if !MFIELD_FIT_DY_NS
          mptr->flags[j] &= ~MAGDATA_FLG_DY_NS;
#endif

#if !MFIELD_FIT_DZ_NS
          mptr->flags[j] &= ~MAGDATA_FLG_DZ_NS;
#endif
        }
    }

  return cnt;
}

/*
mfield_data_init()
  Compute mean and stddev of timestamps minus epoch for later time scaling

w_i = t_i - epoch
t_mu = mean(w_i)
t_sigma = stddev(w_i)

Inputs: w - workspace

Return: success/error

Notes:
1) w->t_mu and w->t_sigma are updated with timestamp mean/stddev in years

2) w->t0_data is initialized to the timestamp of the first data point (CDF_EPOCH)

3) w->t1_data is initialized to the timestamp of the last data point (CDF_EPOCH)

3) w->t0 and w->t1 are initialized to the first/last timestamps of each satellite
*/

int
mfield_data_init(mfield_data_workspace *w)
{
  int s = 0;
  size_t i, j;

  gsl_rstat_reset(w->rstat_workspace_p);

  w->t0_data = 1.0e15;
  w->t1_data = -1.0e15;
  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      magdata_t(&(w->t0[i]), &(w->t1[i]), mptr);

      if (mptr->n != 0)
        {
          if (w->t0[i] > 0.0)
            w->t0_data = GSL_MIN(w->t0_data, w->t0[i]);

          if (w->t1[i] > 0.0)
            w->t1_data = GSL_MAX(w->t1_data, w->t1[i]);
        }

      for (j = 0; j < mptr->n; ++j)
        {
          double t;

          if (mptr->flags[j] & MAGDATA_FLG_DISCARD)
            continue;

          t = satdata_epoch2year(mptr->t[j]) - w->epoch;
          gsl_rstat_add(t, w->rstat_workspace_p);
        }
    }

  w->t_mu = gsl_rstat_mean(w->rstat_workspace_p);
  w->t_sigma = gsl_rstat_sd(w->rstat_workspace_p);

  if (w->t_sigma == 0.0)
    {
      /* this can happen for a fixed time grid like EMAG2 */
      w->t_mu = 0.0;
      w->t_sigma = 1.0;
    }

  return s;
} /* mfield_data_init() */

/*
mfield_data_epoch()
  Compute epoch of input data by averaging all timestamps
*/

double
mfield_data_epoch(mfield_data_workspace *w)
{
  /* initialize t_mu and t_sigma */
  mfield_data_init(w);

  return w->t_mu + w->epoch;
} /* mfield_data_epoch() */

int
mfield_data_map(const char *filename, const mfield_data_workspace *w)
{
  int s = 0;
  size_t i;
  char buf[2048];

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      sprintf(buf, "%s.%zu", filename, i);

      fprintf(stderr, "mfield_data_map: printing spatial coverage of satellite %zu to %s...",
              i, buf);
      magdata_map(buf, mptr);
      fprintf(stderr, "done\n");
    }

  return s;
}

magdata *
mfield_data_ptr(const size_t idx, const mfield_data_workspace *w)
{
  if (idx >= w->nsources)
    {
      fprintf(stderr, "mfield_data_ptr: invalid index: %zu\n", idx);
      return 0;
    }

  return w->mdata[idx];
}
