/*
 * ut1.c
 *
 * This module contains routine for converting from UTC
 * to UT1 time, using UT1-UTC observations read from a data
 * file. Data file is produced from the EOP C04 dataset from:
 *
 * http://hpiers.obspm.fr/eop-pc/index.php?index=C04&lang=en
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <time.h>

#include "ut1.h"

static int ut1_get_val(const time_t t, double *result, ut1_workspace *w);

static int
ut1_read_data(const char *filename, ut1_workspace *w)
{
  int s = 0;
  FILE *fp;
  char buf[UT1_MAX_BUFFER];

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "ut1_read_data: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  while (fgets(buf, UT1_MAX_BUFFER, fp) != NULL)
    {
      int c;
      int year, month, day;
      double dut1;
      ut1_data *ptr;

      if (*buf == '#')
        continue;

      c = sscanf(buf, "%d %d %d %lf", &year, &month, &day, &dut1);
      if (c < 4)
        continue;

      assert(day >= 1 && day <= 31);
      assert(month >= 1 && month <= 12);

      if (year < UT1_FIRST_YEAR)
        {
          fprintf(stderr, "ut1_read_data: invalid year: %d\n", year);
          continue;
        }

      ptr = &(w->ut1[year - UT1_FIRST_YEAR]);

      /* store value and convert to s */
      ptr->data[month - 1][day - 1] = dut1 / 1000.0;
    }

  fclose(fp);

  return s;
}

ut1_workspace *
ut1_alloc(const char *filename)
{
  ut1_workspace *w;
  size_t i, j, k;
  int s;

  w = calloc(1, sizeof(ut1_workspace));
  if (!w)
    return 0;

  /* initialize data */
  for (i = 0; i < UT1_MAX_YEAR; ++i)
    {
      ut1_data *ptr = &(w->ut1[i]);
      for (j = 0; j < 12; ++j)
        {
          for (k = 0; k < 31; ++k)
            ptr->data[j][k] = UT1_BAD_VALUE;
        }
    }

  s = ut1_read_data(filename, w);
  if (s)
    {
      ut1_free(w);
      return 0;
    }

  return w;
} /* ut1_alloc() */

void
ut1_free(ut1_workspace *w)
{
  free(w);
} /* ut1_free() */

/*
ut1_get()
  Look up dUT1 = UT1 - UTC value for a given time

Inputs: t      - timestamp
        result - UT1 - UTC (seconds)
        w      - workspace

Return: success or error
*/

int
ut1_get(const time_t t, double *result, ut1_workspace *w)
{
  int s = 0;
  double val1, val2;
  time_t t2 = t + 86400;

  s = ut1_get_val(t, &val1, w);
  if (s)
    return s;

  s = ut1_get_val(t2, &val2, w);
  if (s)
    return s;

  /**result = interp1d((double)t, (double)t2, val1, val2, );*/

  *result = val1;

  return s;
} /* ut1_get() */

static int
ut1_get_val(const time_t t, double *result, ut1_workspace *w)
{
  int s = 0;
  int year;
  struct tm *tm_p;
  ut1_data *ptr;
  double val;

  tm_p = gmtime(&t);
  year = tm_p->tm_year + 1900;

  if (year < UT1_FIRST_YEAR || year > UT1_FIRST_YEAR + UT1_MAX_YEAR)
    {
      fprintf(stderr, "ut1_get_val: invalid year: %d\n", year);
      return 1;
    }

  ptr = &(w->ut1[year - UT1_FIRST_YEAR]);
  val = ptr->data[tm_p->tm_mon][tm_p->tm_mday - 1];

  if (val == UT1_BAD_VALUE)
    {
      fprintf(stderr, "ut1_get_val: data not available for timestamp %ld\n", t);
      return 1;
    }

  *result = val;

  return s;
} /* ut1_get_val() */
