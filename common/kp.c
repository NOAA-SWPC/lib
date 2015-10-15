/*
 * kp.c
 * Patrick Alken
 *
 * This module reads data from a SPIDR data file containing
 * KP index values and then returns the correct value for
 * a given day.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "kp.h"

/*
read_data()
  Read KP data from a SPIDR data file
*/

static size_t
read_data(const char *filename, kp_workspace *w)
{
  FILE *fp;
  size_t n;
  int ret;
  char buf[KP_MAX_BUFFER];
  char date[11]; /* date in yyyy-mm-dd format */
  char hour[6];  /* hour in hh:mm format */
  double kp;     /* kp value */

  fp = fopen(filename, "r");
  if (!fp)
    {
      perror("fopen");
      return 0;
    }

  n = 0;
  while (fgets(buf, KP_MAX_BUFFER, fp) != NULL)
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
                   &kp);
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
              fprintf(stderr, "read_data: KP data file is corrupted\n");
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

      w->kp_data[n] = kp;

      if (++n >= KP_MAX_DATA)
        {
          fprintf(stderr, "read_data: KP_MAX_DATA not large enough\n");
          return n;
        }
    }

  fclose(fp);

  return n;
} /* read_data() */

/*
kp_alloc()
  Allocate a kp workspace and read in data from filename
*/

kp_workspace *
kp_alloc(const char *filename)
{
  kp_workspace *w;
  size_t i;

  w = (kp_workspace *) calloc(1, sizeof(kp_workspace));
  if (!w)
    return 0;

  w->kp_data = (double *) malloc(sizeof(double) * KP_MAX_DATA);
  for (i = 0; i < KP_MAX_DATA; ++i)
    w->kp_data[i] = 0.0;

  w->t0 = 0;

  w->n = read_data(filename, w);

  return (w);
} /* kp_alloc() */

void
kp_free(kp_workspace *w)
{
  if (!w)
    return;

  if (w->kp_data)
    free(w->kp_data);

  free(w);
} /* kp_free() */

/*
kp_get()
  Obtain the KP value for a given timestamp

Inputs: t      - timestamp in UTC seconds since epoch
        result - where to store KP value
        w      - kp workspace

Return: 0 on success, 1 on failure
*/

int
kp_get(time_t t, double *result, kp_workspace *w)
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
      fprintf(stderr, "kp_get: data not available for given day (t = %d)\n", t);
      return 1;
    }

  /*
   * there are 8 KP readings each day, so determine which
   * eighth of the day this time falls into
   */
  hours = (double)dt / 86400.0 - (double) days;
  hours *= 24.0;
  offset = (size_t) round(8.0 / 24.0 * hours);

  *result = w->kp_data[idx + offset];

  return 0;
} /* kp_get() */
