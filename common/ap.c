/*
 * ap.c
 * Patrick Alken
 *
 * This module reads data from a SPIDR data file containing
 * AP index values and then returns the correct value for
 * a given day.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <string.h>

#include "ap.h"

/*
read_data()
  Read AP data from a SPIDR data file
*/

static size_t
read_data(const char *filename, ap_workspace *w)
{
  FILE *fp;
  size_t n;
  int ret;
  char buf[AP_MAX_BUFFER];
  char date[11]; /* date in yyyy-mm-dd format */
  char hour[6];  /* hour in hh:mm format */
  double ap;     /* ap value */

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "fopen: cannot open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  n = 0;
  while (fgets(buf, AP_MAX_BUFFER, fp) != NULL)
    {
      /* ignore comments */
      if (*buf == '#')
        continue;

      /*
       * file format is:
       *
       * yyyy-mm-dd hh:mm value qualifier description
       *
       * with eight data values per day (taken at 00:00, 03:00,
       * 06:00, 09:00, 12:00, 15:00, 18:00, 21:00 UTC)
       */
      ret = sscanf(buf,
                   "%10s %5s %lf%*[^\n]",
                   date,
                   hour,
                   &ap);
      if (ret <= 0)
        continue;

      if (w->t0 == 0)
        {
          struct tm tm0;
          int year, month, day;
          char tz_str[] = "TZ=GMT";

          ret = sscanf(date, "%d-%d-%d", &year, &month, &day);
          if (ret <= 0)
            {
              fprintf(stderr, "read_data: AP data file is corrupted\n");
              return 0;
            }

          /*
           * this is the first data line - record the first
           * date timestamp
           */
          tm0.tm_sec = 0;
          tm0.tm_min = 0;
          tm0.tm_hour = 0;
          tm0.tm_mday = day;
          tm0.tm_mon = month - 1;
          tm0.tm_year = year - 1900;
          tm0.tm_isdst = 0;

          putenv(tz_str);
          w->t0 = mktime(&tm0);
        }

      w->ap_data[n] = ap;

      if (++n >= AP_MAX_DATA)
        {
          fprintf(stderr, "read_data: AP_MAX_DATA not large enough\n");
          return n;
        }
    }

  fclose(fp);

  return n;
} /* read_data() */

/*
ap_alloc()
  Allocate a ap workspace and read in data from filename
*/

ap_workspace *
ap_alloc(const char *filename)
{
  ap_workspace *w;
  size_t i;

  w = (ap_workspace *) calloc(1, sizeof(ap_workspace));
  if (!w)
    return 0;

  w->ap_data = (double *) malloc(sizeof(double) * AP_MAX_DATA);
  for (i = 0; i < AP_MAX_DATA; ++i)
    w->ap_data[i] = 0.0;

  w->t0 = 0;

  w->n = read_data(filename, w);

  return (w);
} /* ap_alloc() */

void
ap_free(ap_workspace *w)
{
  if (!w)
    return;

  if (w->ap_data)
    free(w->ap_data);

  free(w);
} /* ap_free() */

/*
ap_get()
  Obtain the AP value for a given timestamp

Inputs: t      - timestamp in UTC seconds since epoch
        result - where to store AP value
        w      - ap workspace

Return: 0 on success, 1 on failure
*/

int
ap_get(time_t t, double *result, ap_workspace *w)
{
  size_t days;
  size_t idx;
  size_t offset;
  size_t dt;
  double hours;

  dt = (size_t) (t - w->t0);

  days = dt / 86400;
  idx = days * 8;

  if ((t < w->t0) || (idx >= w->n))
    {
      fprintf(stderr, "ap_get: data not available for given t (%ld)\n", t);
      return 1;
    }

  /*
   * there are 8 AP readings each day, so determine which
   * eighth of the day this time falls into
   */
  hours = (double)dt / 86400.0 - (double) days;
  hours *= 24.0;
  offset = (size_t) round(8.0 / 24.0 * hours);

  *result = w->ap_data[idx + offset];

  return 0;
} /* ap_get() */
