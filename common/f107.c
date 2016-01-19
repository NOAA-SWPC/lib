/*
 * f107.c
 *
 * This module reads data from a F10.7 data file containing
 * F10.7 index values and then returns the correct value for
 * a given day.
 *
 * The data files are taken from the SPIDR database
 *
 * The file format is described in the header of the file
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "common.h"
#include "f107.h"

/*
f107_read_data()
  Read F107 data from a Penticton data file

Inputs: filename - F10.7 data file
        w        - f107 workspace

Notes: there are a few days (notably in 2008) which have missing data;
these are stored in the data structure as -1.0 values
*/

static int
f107_read_data(const char *filename, f107_workspace *w)
{
  int s = 0;
  FILE *fp;
  int ret;
  char buf[F107_MAX_BUFFER];
  int day;
  f107_data *fptr;
  int year, month;
  double f107; /* F10.7 value */
  char str[F107_MAX_BUFFER];

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "fopen: cannot open %s: %s\n",
              filename, strerror(errno));
      return 1;
    }

  while (fgets(buf, F107_MAX_BUFFER, fp) != NULL)
    {
      if (*buf == '#')
        continue;

      /*
       * file format is:
       *
       * yyyy-mm-dd hh:mm f107
       */

      /* read day number */
      ret = sscanf(buf, "%d-%d-%d %05s %lf", &year, &month, &day, str, &f107);
      if (ret < 5)
        continue;

      assert(day >= 1 && day <= 31);
      assert(month >= 1 && month <= 12);

      if (year < F107_FIRST_YEAR)
        {
          fprintf(stderr, "f107_read_data: invalid year: %d\n", year);
          continue;
        }

      fptr = &(w->f107[year - F107_FIRST_YEAR]);

      /* store value */
      fptr->data[month - 1][day - 1] = f107;
    }

  fclose(fp);

  return s;
} /* f107_read_data() */

/*
f107_alloc()
  Allocate a f107 workspace and read in data from filename

Inputs: filename - F10.7 data file

Return: pointer to workspace
*/

f107_workspace *
f107_alloc(const char *filename)
{
  f107_workspace *w;
  size_t i, j, k;

  w = (f107_workspace *) calloc(1, sizeof(f107_workspace));
  if (!w)
    return 0;

  for (i = 0; i < F107_MAX_YEAR; ++i)
    {
      f107_data *fptr = &(w->f107[i]);
      for (j = 0; j < 12; ++j)
        for (k = 0; k < 31; ++k)
          fptr->data[j][k] = -1.0;
    }

  /* read in available NOAA-formatted data */
  f107_read_data(filename, w);

  return (w);
} /* f107_alloc() */

void
f107_free(f107_workspace *w)
{
  if (!w)
    return;

  free(w);
} /* f107_free() */

/*
f107_get()
  Obtain the F10.7 value for a given timestamp

Inputs: t      - UT timestamp
        result - where to store F10.7 value
        w      - f107 workspace

Return: 0 on success, 1 on failure
*/

int
f107_get(time_t t, double *result, f107_workspace *w)
{
  int year;
  struct tm *tm_p;
  double val;
  f107_data *fptr;

  tm_p = gmtime(&t);
  year = tm_p->tm_year + 1900;
  fptr = &(w->f107[year - F107_FIRST_YEAR]);

  val = fptr->data[tm_p->tm_mon][tm_p->tm_mday - 1];
  if (val < 0.0)
    {
      fprintf(stderr, "f107_get: data not found for timestamp %ld: %s",
              t, asctime(tm_p));
      return 1;
    }

  *result = val;

  return 0;
} /* f107_get() */

/*
f107a_get()
  Obtain the F10.7A value for a given day (the 81-day centered average of
F10.7)

Inputs: t      - timestamp
        result - where to store F10.7A value
        w      - f107 workspace

Return: 0 on success, 1 on failure
*/

int
f107a_get(time_t t, double *result, f107_workspace *w)
{
  double sum, tmp;
  time_t tmin, tmax, tt;
  int wlen = 81; /* number of days to average */
  int s;

  tmax = t + (wlen/2) * 86400;
  tmin = tmax - (wlen - 1) * 86400;

  sum = 0.0;
  for (tt = tmin; tt <= tmax; tt += 86400)
    {
      s = f107_get(tt, &tmp, w);
      if (s)
        continue;

      sum += tmp;
    }

  *result = sum / (double) wlen;

  return 0;
} /* f107a_get() */

/*
f107a_get2()
  Obtain the F10.7A value for a given day (the 27-day backward average of
F10.7)

Inputs: t      - timestamp
        result - where to store F10.7A value
        w      - f107 workspace

Return: 0 on success, 1 on failure
*/

int
f107a_get2(time_t t, double *result, f107_workspace *w)
{
  double sum, tmp;
  time_t tmin, tmax, tt;
  int wlen = 27; /* number of days to average */
  int s;

  /* backward average */
  tmax = t;
  tmin = tmax - (wlen - 1) * 86400;

  sum = 0.0;
  for (tt = tmin; tt <= tmax; tt += 86400)
    {
      s = f107_get(tt, &tmp, w);
      if (s)
        continue;

      sum += tmp;
    }

  *result = sum / (double) wlen;

  return 0;
} /* f107a_get2() */
