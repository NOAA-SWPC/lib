/*
 * jicmag.c
 * Patrick Alken
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

#include "jicmag.h"

static int jicmag_compare(const void *a, const void *b);

/*
jicmag_read_data()
  Read data from a Jicamarca data file

Inputs: filename - data file
        w        - kp workspace
*/

int
jicmag_read_data(const char *filename, jicmag_data *data)
{
  int s = 0;
  FILE *fp;
  char buf[JICMAG_MAX_BUFFER];
  size_t n;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "jicmag_read_data: fopen: cannot open %s: %s\n",
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
      int dd, mm, year, hh, min;
      double D, H;

      c = sscanf(buf, "%d %d %d %d %d %lf %lf",
                 &dd,
                 &mm,
                 &year,
                 &hh,
                 &min,
                 &D,
                 &H);
      if (c < 7)
        continue;

      tmp.tm_sec = 0;
      tmp.tm_min = min;
      tmp.tm_hour = hh;
      tmp.tm_mday = dd;
      tmp.tm_mon = mm - 1;
      tmp.tm_year = year - 1900;

      data->array[n].t = mktime(&tmp);
      data->array[n].H = H;

      if (++n >= data->ntot)
        {
          fprintf(stderr, "jicmag_read_data: error: ntot too small\n");
          exit(1);
        }
    }

  fclose(fp);

  data->n = n;

  return s;
} /* jicmag_read_data() */

/*
jicmag_read_idx()
  Read all available Jicamarca data

Inputs: idx_filename - index file
*/

jicmag_data *
jicmag_read_idx(const char *idx_filename)
{
  int s;
  FILE *fp;
  char buffer[JICMAG_MAX_BUFFER];
  jicmag_data *data;

  fp = fopen(idx_filename, "r");
  if (!fp)
    {
      fprintf(stderr, "jicmag_read_idx: unable to open index file %s: %s\n",
              idx_filename, strerror(errno));
      return 0;
    }

  data = jicmag_alloc(JICMAG_MAX_DATA);

  while (fgets(buffer, sizeof(buffer), fp) != 0)
    {
      buffer[strlen(buffer) - 1] = '\0';
      s = jicmag_read_data(buffer, data);
      if (s)
        break;
    }

  fclose(fp);

  return data;
} /* jicmag_read_idx() */

jicmag_data *
jicmag_alloc(const size_t n)
{
  jicmag_data *data;

  data = calloc(1, sizeof(jicmag_data));
  if (!data)
    return 0;

  data->array = malloc(n * sizeof(jicmag_datum));

  data->ntot = n;
  data->n = 0;

  return data;
} /* jicmag_alloc() */

void
jicmag_free(jicmag_data *data)
{
  if (data->array)
    free(data->array);

  free(data);
} /* jicmag_free() */

int
jicmag_sort(jicmag_data *data)
{
  qsort(data->array, data->n, sizeof(jicmag_datum), jicmag_compare);
  return 0;
} /* jicmag_sort() */

static int
jicmag_compare(const void *a, const void *b)
{
  jicmag_datum *da = (jicmag_datum *) a;
  jicmag_datum *db = (jicmag_datum *) b;
  const time_t t1 = da->t;
  const time_t t2 = db->t;

  return (t1 - t2);
} /* jicmag_compare() */
