/*
 * magdata_list.c
 *
 * This is a high-level wrapper for multiple magdata sources
 * (multiple satellites, ground observatories, etc)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <assert.h>

#include <satdata/satdata.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <common/common.h>

#include "magdata_list.h"

/*
magdata_list_alloc()
  Allocate magdata_list

Inputs: nsources - number of data sources (satellites, observatories)

Return: pointer to workspace
*/

magdata_list *
magdata_list_alloc(const size_t nsources)
{
  magdata_list *w;

  w = calloc(1, sizeof(magdata_list));
  if (!w)
    return 0;

  w->n_tot = nsources;
  w->n = 0;

  w->mdata = calloc(1, nsources * sizeof(magdata *));
  if (!w->mdata)
    {
      magdata_list_free(w);
      return 0;
    }

  return w;
}

void
magdata_list_free(magdata_list *w)
{
  if (w->mdata)
    {
      size_t i;

      for (i = 0; i < w->n; ++i)
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
magdata_list_add
  Read a magdata file and add it to the list of sources

Inputs: filename - magdata file
        w        - workspace

Return: number of data read
*/

size_t
magdata_list_add(const char *filename, magdata_list *w)
{
  if (w->n >= w->n_tot)
    {
      fprintf(stderr, "magdata_list_add: error: only %zu data sources allocated\n", w->n_tot);
      return 0;
    }
  else
    {
      size_t ndata;

      w->mdata[w->n] = magdata_read(filename, NULL);
      if (w->mdata[w->n] == NULL)
        {
          fprintf(stderr, "magdata_list_add: error reading %s\n", filename);
          return 0;
        }

      ndata = w->mdata[w->n]->n;

      ++(w->n);

      return ndata;
    }
}

/*
magdata_list_filter_time()
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
magdata_list_filter_time(const double tmin, const double tmax,
                         magdata_list *w)
{
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, w);

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
}

/*
magdata_list_filter_euler()
  We are not fitting Euler angles (fit_euler is 0),
so discard any data which is marked Euler only

Inputs: w    - workspace

Return: number of data flagged
*/

size_t
magdata_list_filter_euler(magdata_list *w)
{
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, w);

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
}

int
magdata_list_map(const char *dir_prefix, const magdata_list *w)
{
  int s = 0;
  size_t i;
  char buf[2048];

  for (i = 0; i < w->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, w);

      sprintf(buf, "%s/map%zu", dir_prefix, i);
      magdata_map(buf, mptr);
    }

  return s;
}

/*
magdata_list_print()
  Print out all data which will be used for field modeling

Inputs: dir_prefix - directory prefix
        w          - magdata_list
*/

int
magdata_list_print(const char *dir_prefix, const magdata_list *w)
{
  int s = 0;
  size_t i;
  char buf[2048];

  for (i = 0; i < w->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, w);

      sprintf(buf, "%s/data%zu", dir_prefix, i);
      magdata_print(buf, mptr);
    }

  return s;
}

magdata *
magdata_list_ptr(const size_t idx, const magdata_list *w)
{
  if (idx >= w->n)
    {
      fprintf(stderr, "magdata_list_ptr: invalid index: %zu\n", idx);
      return 0;
    }

  return w->mdata[idx];
}

/*
magdata_list_rminmax()
  Compute minimum and maximum radius of all magdata sources

Inputs: list - magdata list
        rmin - (output) minimum radius (km)
        rmax - (output) maximum radius (km)

Return: success/error
*/

int
magdata_list_rminmax(const magdata_list * list, double * rmin, double * rmax)
{
  int s = 0;
  size_t i;

  *rmin = 1.0e9;
  *rmax = -1.0e9;

  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);

      *rmin = GSL_MIN(*rmin, mptr->rmin);
      *rmax = GSL_MAX(*rmax, mptr->rmax);
    }

  return s;
}

/*
magdata_list_index()
  Count number of data of each type and build index for each
data point

Inputs: list  - magdata list
        count - (output) counts of each data type, indexed by
                MAGDATA_LIST_IDX_xxx (size MAGDATA_LIST_IDX_END)

Return: success/error

Notes:
1) The mptr->index[j] array is updated to contain the index of
mptr->datum[j] (the jth point in mptr in [0,data_total - 1])
*/

int
magdata_list_index(magdata_list * list, size_t count[])
{
  int s = 0;
  size_t i, j, k;

  for (i = 0; i < MAGDATA_LIST_IDX_END; ++i)
    count[i] = 0;

  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);

      for (j = 0; j < mptr->n; ++j)
        {
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* store index of this data point accounting for all previous counts */
          mptr->index[j] = 0;
          for (k = 0; k <= MAGDATA_LIST_IDX_DF_EW; ++k)
            {
              mptr->index[j] += count[k];
            }

          if (MAGDATA_ExistX(mptr->flags[j]))
            ++count[MAGDATA_LIST_IDX_X];

          if (MAGDATA_ExistY(mptr->flags[j]))
            ++count[MAGDATA_LIST_IDX_Y];

          if (MAGDATA_ExistZ(mptr->flags[j]))
            ++count[MAGDATA_LIST_IDX_Z];

          if (MAGDATA_ExistScalar(mptr->flags[j]))
            ++count[MAGDATA_LIST_IDX_F];

          if (MAGDATA_ExistDX_NS(mptr->flags[j]))
            ++count[MAGDATA_LIST_IDX_DX_NS];

          if (MAGDATA_ExistDY_NS(mptr->flags[j]))
            ++count[MAGDATA_LIST_IDX_DY_NS];

          if (MAGDATA_ExistDZ_NS(mptr->flags[j]))
            ++count[MAGDATA_LIST_IDX_DZ_NS];

          if (MAGDATA_ExistDF_NS(mptr->flags[j]))
            ++count[MAGDATA_LIST_IDX_DF_NS];

          if (MAGDATA_ExistDX_EW(mptr->flags[j]))
            ++count[MAGDATA_LIST_IDX_DX_EW];

          if (MAGDATA_ExistDY_EW(mptr->flags[j]))
            ++count[MAGDATA_LIST_IDX_DY_EW];

          if (MAGDATA_ExistDZ_EW(mptr->flags[j]))
            ++count[MAGDATA_LIST_IDX_DZ_EW];

          if (MAGDATA_ExistDF_EW(mptr->flags[j]))
            ++count[MAGDATA_LIST_IDX_DF_EW];
        }
    }

  /* compute totals */
  for (i = 0; i < MAGDATA_LIST_IDX_DF_EW; ++i)
    count[MAGDATA_LIST_IDX_TOTAL] += count[i];

  count[MAGDATA_LIST_IDX_VECTOR] = count[MAGDATA_LIST_IDX_X] +
                                   count[MAGDATA_LIST_IDX_Y] +
                                   count[MAGDATA_LIST_IDX_Z];

  count[MAGDATA_LIST_IDX_VGRAD] = count[MAGDATA_LIST_IDX_DX_NS] +
                                  count[MAGDATA_LIST_IDX_DY_NS] +
                                  count[MAGDATA_LIST_IDX_DZ_NS] +
                                  count[MAGDATA_LIST_IDX_DX_EW] +
                                  count[MAGDATA_LIST_IDX_DY_EW] +
                                  count[MAGDATA_LIST_IDX_DZ_EW];

  return s;
}
