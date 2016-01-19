/*
 * julia.c
 *
 * This module reads Jicamarca magnetometer data files
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "julia.h"

static int julia_compare_sort(const void *a, const void *b);
static int julia_compare_search(const void *a, const void *b);

/* XXX Global variable */
double JULIA_TIME_WINDOW = 0.0;

/*
julia_read_data()
  Read data from a Jicamarca data file

Inputs: filename - data file
        w        - kp workspace
*/

int
julia_read_data(const char *filename, julia_data *data)
{
  int s = 0;
  FILE *fp;
  char buf[JULIA_MAX_BUFFER];
  size_t n;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "julia_read_data: fopen: cannot open %s: %s\n",
              filename, strerror(errno));
      return 1;
    }

  /* use GMT time */
  putenv("TZ=GMT");

  n = data->n;

  while (fgets(buf, sizeof(buf), fp) != NULL)
    {
      int c;
      struct tm tmp;
      int dd, mm, year, hh, min, sec;
      double glat, glon, galt, snl;
      double vipe1, dvipe1, vipn1, dvipn1;

      c = sscanf(buf, "%d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf",
                 &year,
                 &mm,
                 &dd,
                 &hh,
                 &min,
                 &sec,
                 &glat,
                 &glon,
                 &galt,
                 &snl,
                 &vipe1,
                 &dvipe1,
                 &vipn1,
                 &dvipn1);
      if (c < 14)
        continue;

      tmp.tm_sec = sec;
      tmp.tm_min = min;
      tmp.tm_hour = hh;
      tmp.tm_mday = dd;
      tmp.tm_mon = mm - 1;
      tmp.tm_year = year - 1900;

      data->array[n].t = mktime(&tmp);
      data->array[n].v_zonal = vipe1;
      data->array[n].v_vert = vipn1;

      if (++n >= data->ntot)
        {
          fprintf(stderr, "julia_read_data: error: ntot too small\n");
          exit(1);
        }
    }

  fclose(fp);

  data->n = n;

  return s;
} /* julia_read_data() */

/*
julia_read_idx()
  Read all available Jicamarca data

Inputs: idx_filename - index file
        data         - (output) append data here
                       or set to NULL for new
                       allocation
*/

julia_data *
julia_read_idx(const char *idx_filename, julia_data *data)
{
  int s;
  FILE *fp;
  char buffer[JULIA_MAX_BUFFER];

  fp = fopen(idx_filename, "r");
  if (!fp)
    {
      fprintf(stderr, "julia_read_idx: unable to open index file %s: %s\n",
              idx_filename, strerror(errno));
      return 0;
    }

  if (data == NULL)
    data = julia_alloc(JULIA_MAX_DATA);

  while (fgets(buffer, sizeof(buffer), fp) != 0)
    {
      buffer[strlen(buffer) - 1] = '\0';
      s = julia_read_data(buffer, data);
      if (s)
        break;
    }

  fclose(fp);

  return data;
} /* julia_read_idx() */

/*
julia_read_avg_data()
  Read a JULIA avg data file and add to data structure
*/

int
julia_read_avg_data(const char *filename, julia_data *data)
{
  int s = 0;
  char buffer[JULIA_MAX_BUFFER];
  FILE *fp;
  int year, month, day, hour, min, sec;
  double glat, glon, inttms, alt, dalt, vipn2, dvipn2, vipe1, dvipe1;
  size_t n = data->n;
  struct tm tm_p;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "julia_read_avg_data: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  putenv("TZ=GMT");

  while (fgets(buffer, sizeof(buffer), fp) != 0)
    {
      int c;

      c = sscanf(buffer, "%d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                 &year,
                 &month,
                 &day,
                 &hour,
                 &min,
                 &sec,
                 &glat,
                 &glon,
                 &inttms,
                 &alt,
                 &dalt,
                 &vipn2,
                 &dvipn2,
                 &vipe1,
                 &dvipe1);
      if (c < 15)
        continue;

      tm_p.tm_sec = sec;
      tm_p.tm_min = min;
      tm_p.tm_hour = hour;
      tm_p.tm_mday = day;
      tm_p.tm_mon = month - 1;
      tm_p.tm_year = year - 1900;
      tm_p.tm_isdst = 0;

      data->array[n].t = mktime(&tm_p);
      data->array[n].v_zonal = vipe1;
      data->array[n].v_vert = vipn2;

      if (++n >= data->ntot)
        {
          fprintf(stderr, "julia_read_avg_data: ntot too small (%zu)\n", data->ntot);
          break;
        }
    }

  fclose(fp);

  data->n = n;

  return s;
} /* julia_read_avg_data() */

/*
julia_read_avg_idx()
  Read all Jicamarca avg data from index file

Inputs: idx_filename - index file
*/

julia_data *
julia_read_avg_idx(const char *idx_filename, julia_data *data)
{
  int s;
  FILE *fp;
  char buffer[JULIA_MAX_BUFFER];

  fp = fopen(idx_filename, "r");
  if (!fp)
    {
      fprintf(stderr, "julia_read_avg_idx: unable to open index file %s: %s\n",
              idx_filename, strerror(errno));
      return 0;
    }

  if (data == NULL)
    data = julia_alloc(JULIA_MAX_DATA);

  while (fgets(buffer, sizeof(buffer), fp) != 0)
    {
      buffer[strlen(buffer) - 1] = '\0';
      s = julia_read_avg_data(buffer, data);
      if (s)
        break;
    }

  fclose(fp);

  return data;
} /* julia_read_avg_idx() */

julia_data *
julia_alloc(const size_t n)
{
  julia_data *data;

  data = calloc(1, sizeof(julia_data));
  if (!data)
    return 0;

  data->array = malloc(n * sizeof(julia_datum));

  data->ntot = n;
  data->n = 0;

  return data;
} /* julia_alloc() */

void
julia_free(julia_data *data)
{
  if (data->array)
    free(data->array);

  free(data);
} /* julia_free() */

int
julia_sort(julia_data *data)
{
  qsort(data->array, data->n, sizeof(julia_datum), julia_compare_sort);
  return 0;
} /* julia_sort() */

/*
julia_search()
  Search JULIA data for a measurement within a specified window of a
given time

Inputs: t      - timestamp to search for
        window - allowed time window
        data   - JULIA data
        index  - (output) index into array of data if found

Return: 0 if found, -1 if not found
*/

int
julia_search(const time_t t, const double window, const julia_data *data,
             size_t *index)
{
  int s = 0;
  julia_datum *ptr;

  /* XXX Global variable */
  JULIA_TIME_WINDOW = window;

  ptr = bsearch(&t, data->array, data->n, sizeof(julia_datum),
                julia_compare_search);
  if (ptr)
    {
      size_t idx = ptr - data->array;
      assert(ptr->t == data->array[idx].t);
      *index = idx;
    }
  else
    s = -1; /* not found */

  return s;
}

static int
julia_compare_sort(const void *a, const void *b)
{
  julia_datum *da = (julia_datum *) a;
  julia_datum *db = (julia_datum *) b;
  const time_t t1 = da->t;
  const time_t t2 = db->t;

  return (t1 - t2);
} /* julia_compare_sort() */

static int
julia_compare_search(const void *a, const void *b)
{
  const time_t t1 = *(const time_t *) a;
  julia_datum *db = (julia_datum *) b;
  const time_t t2 = db->t;
  double diff = (t2 - t1) / 60.0;

  if (fabs(diff) < JULIA_TIME_WINDOW)
    return 0; /* found */
  else if (t1 < t2)
    return -1;
  else
    return 1;
} /* julia_compare_search() */
