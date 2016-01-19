/*
 * dst.c
 *
 * This module reads data from the automatically updated Dst/Est/Ist
 * data file in /nfs/satmag/data/Indices/Dst/Est_Ist_index.pli
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "common.h"
#include "dst.h"
#include "interp.h"

/*
dst_read_data()
  Read Dst/Est/Ist data from a SPIDR data file
*/

static size_t
dst_read_data(const char *filename, dst_workspace *w)
{
  FILE *fp;
  size_t n;
  int ret;
  char buf[DST_MAX_BUFFER];
  double dst;    /* dst value */
  double est;    /* est value */
  double ist;    /* ist value */
  double fday;
  time_t t;
  struct tm *tm_p;
  int year;
  int yidx, midx, didx;

  fp = fopen(filename, "r");
  if (!fp)
    {
      perror("fopen");
      return 0;
    }

  n = 0;
  while (fgets(buf, DST_MAX_BUFFER, fp) != NULL)
    {
      /* ignore comments */
      if (*buf == '#')
        continue;

      /*
       * file format is:
       *
       * yyyy-mm-dd hh:mm value qualifier description
       *
       * with hourly sampling
       */
      ret = sscanf(buf, "%lf %lf %lf %lf\n",
                   &fday, &dst, &est, &ist);
      if (ret < 4)
        continue;

      t = fday2timet(fday);
      tm_p = gmtime(&t);
      year = tm_p->tm_year + 1900;

      if (year < DST_FIRST_YEAR)
        {
          fprintf(stderr, "dst_read_data: invalid year: %d\n", year);
          continue;
        }

      yidx = year - DST_FIRST_YEAR;
      midx = tm_p->tm_mon;
      didx = tm_p->tm_mday - 1;

      assert(midx >= 0 && midx < 12);
      assert(didx >= 0 && didx < 31);

      n = w->data[yidx][midx][didx].n;
      w->data[yidx][midx][didx].dst[n] = dst;
      w->data[yidx][midx][didx].est[n] = est;
      w->data[yidx][midx][didx].ist[n] = ist;
      w->data[yidx][midx][didx].t[n] = t;

      if (++n >= DST_MAX_DATA)
        {
          fprintf(stderr, "dst_read_data: DST_MAX_DATA not large enough\n");
          return n;
        }

      w->data[yidx][midx][didx].n = n;
    }

  fclose(fp);

  return n;
} /* dst_read_data() */

/*
dst_alloc()
  Allocate a dst workspace and read in data from filename
*/

dst_workspace *
dst_alloc(const char *filename)
{
  dst_workspace *w;
  size_t i, j, k, l;

  w = calloc(1, sizeof(dst_workspace));
  if (!w)
    return 0;

  for (i = 0; i < DST_MAX_YEAR; ++i)
    {
      for (j = 0; j < 12; ++j)
        {
          for (k = 0; k < 31; ++k)
            {
              w->data[i][j][k].dst = malloc(sizeof(double) * DST_MAX_DATA);
              w->data[i][j][k].est = malloc(sizeof(double) * DST_MAX_DATA);
              w->data[i][j][k].ist = malloc(sizeof(double) * DST_MAX_DATA);
              w->data[i][j][k].t = malloc(sizeof(time_t) * DST_MAX_DATA);

              for (l = 0; l < DST_MAX_DATA; ++l)
                {
                  w->data[i][j][k].dst[l] = 0.0;
                  w->data[i][j][k].est[l] = 0.0;
                  w->data[i][j][k].ist[l] = 0.0;
                  w->data[i][j][k].t[l] = 0;
                  w->data[i][j][k].n = 0;
                }
            }
        }
    }

  w->n = dst_read_data(filename, w);

  return (w);
} /* dst_alloc() */

void
dst_free(dst_workspace *w)
{
  size_t i, j, k;

  for (i = 0; i < DST_MAX_YEAR; ++i)
    {
      for (j = 0; j < 12; ++j)
        {
          for (k = 0; k < 31; ++k)
            {
              dst_data *dptr = &(w->data[i][j][k]);

              if (dptr->dst)
                free(dptr->dst);

              if (dptr->est)
                free(dptr->est);

              if (dptr->ist)
                free(dptr->ist);

              if (dptr->t)
                free(dptr->t);
            }
        }
    }

  free(w);
} /* dst_free() */

/*
dst_get()
  Obtain the DST value for a given timestamp

Inputs: t      - timestamp in UTC seconds since epoch
        result - where to store DST value (in nT)
        w      - dst workspace

Return: 0 on success, 1 on failure
*/

int
dst_get(time_t t, double *result, dst_workspace *w)
{
  int s = 0;
  struct tm *tm_p;
  int yidx, midx, didx;
  size_t i;
  double val = 999.0;
  size_t n;
  dst_data *dptr;

  tm_p = gmtime(&t);
  yidx = tm_p->tm_year + 1900 - DST_FIRST_YEAR;
  midx = tm_p->tm_mon;
  didx = tm_p->tm_mday - 1;

  dptr = &(w->data[yidx][midx][didx]);

  n = dptr->n;
  if (n == 0)
    {
      fprintf(stderr, "dst_get: no data available for given timestamp (%ld)\n", t);
      return 1;
    }

  if (t < dptr->t[0])
    {
      *result = dptr->dst[0];
      return 0;
    }
  else if (t >= dptr->t[n - 1])
    {
      *result = dptr->dst[n - 1];
      return 0;
    }

  for (i = 0; i < n - 1; ++i)
    {
      if (dptr->t[i] <= t && t < dptr->t[i + 1])
        {
          /*val = dptr->dst[i];*/
          val = interp1d(dptr->t[i], dptr->t[i + 1],
                         dptr->dst[i], dptr->dst[i + 1], t);
          break;
        }
    }

  if (val > 990.0)
    {
      fprintf(stderr, "dst_get: unable to find data for timestamp (%ld)\n", t);
      return 1;
    }

  *result = val;

  return s;
} /* dst_get() */

/*
est_get()
  Obtain the EST value for a given timestamp

Inputs: t      - timestamp in UTC seconds since epoch
        result - where to store DST value (in nT)
        w      - dst workspace

Return: 0 on success, 1 on failure
*/

int
est_get(time_t t, double *result, dst_workspace *w)
{
  int s = 0;
  struct tm *tm_p;
  int yidx, midx, didx;
  size_t i;
  double val = 999.0;
  size_t n;
  dst_data *dptr;

  tm_p = gmtime(&t);
  yidx = tm_p->tm_year + 1900 - DST_FIRST_YEAR;
  midx = tm_p->tm_mon;
  didx = tm_p->tm_mday - 1;

  dptr = &(w->data[yidx][midx][didx]);

  n = dptr->n;
  if (n == 0)
    {
      fprintf(stderr, "est_get: no data available for given timestamp (%ld)\n", t);
      return 1;
    }

  if (t < dptr->t[0])
    {
      *result = dptr->est[0];
      return 0;
    }
  else if (t >= dptr->t[n - 1])
    {
      *result = dptr->est[n - 1];
      return 0;
    }

  for (i = 0; i < n - 1; ++i)
    {
      if (dptr->t[i] <= t && t < dptr->t[i + 1])
        {
          /*val = dptr->est[i];*/
          val = interp1d(dptr->t[i], dptr->t[i + 1],
                         dptr->est[i], dptr->est[i + 1], t);
          break;
        }
    }

  if (val > 990.0)
    {
      fprintf(stderr, "est_get: unable to find data for timestamp (%ld)\n", t);
      return 1;
    }

  *result = val;

  return s;
} /* est_get() */

/*
ist_get()
  Obtain the IST value for a given timestamp

Inputs: t      - timestamp in UTC seconds since epoch
        result - where to store DST value (in nT)
        w      - dst workspace

Return: 0 on success, 1 on failure
*/

int
ist_get(time_t t, double *result, dst_workspace *w)
{
  int s = 0;
  struct tm *tm_p;
  int yidx, midx, didx;
  size_t i;
  double val = 999.0;
  size_t n;
  dst_data *dptr;

  tm_p = gmtime(&t);
  yidx = tm_p->tm_year + 1900 - DST_FIRST_YEAR;
  midx = tm_p->tm_mon;
  didx = tm_p->tm_mday - 1;

  dptr = &(w->data[yidx][midx][didx]);

  n = dptr->n;
  if (n == 0)
    {
      fprintf(stderr, "ist_get: no data available for given timestamp (%ld)\n", t);
      return 1;
    }

  if (t < dptr->t[0])
    {
      *result = dptr->ist[0];
      return 0;
    }
  else if (t >= dptr->t[n - 1])
    {
      *result = dptr->ist[n - 1];
      return 0;
    }

  for (i = 0; i < n - 1; ++i)
    {
      if (dptr->t[i] <= t && t < dptr->t[i + 1])
        {
          /*val = dptr->ist[i];*/
          val = interp1d(dptr->t[i], dptr->t[i + 1],
                         dptr->ist[i], dptr->ist[i + 1], t);
          break;
        }
    }

  if (val > 990.0)
    {
      fprintf(stderr, "ist_get: unable to find data for timestamp (%ld)\n", t);
      return 1;
    }

  *result = val;

  return s;
} /* ist_get() */
