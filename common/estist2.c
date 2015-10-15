/*
 * estist2.c
 * Patrick Alken
 *
 * This module reads Est/Ist data files and then returns the correct value
 * for a given timestamp.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "estist2.h"

#include "common.h"

static long estist2_get_val(const char *buf, int bcol, int ecol);
static int estist2_read_data(const char *filename, estist2_workspace *w);

/*
estist2_read_data()
  Read estist2 data from data file

Inputs: filename - Est/Ist data file
        w        - estist2 workspace
*/

static int
estist2_read_data(const char *filename, estist2_workspace *w)
{
  int s = 0;
  FILE *fp;
  char buf[ESTIST2_MAX_BUFFER];
  estist2_data *dptr;
  int year, month, day;
  int hh;
  int c;
  time_t t;
  struct tm *tm_p;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "estist2_read_data: fopen: cannot open %s: %s\n",
              filename, strerror(errno));
      return 1;
    }

  while (fgets(buf, ESTIST2_MAX_BUFFER, fp) != NULL)
    {
      double I_st, E_st, fday, dst;

      if (*buf == '#')
        continue;

      c = sscanf(buf, "%lf %lf %lf %lf",
                 &fday,
                 &dst,
                 &E_st,
                 &I_st);
      if (c < 4)
        continue;

      t = fday2timet(fday);

      tm_p = gmtime(&t);
      hh = tm_p->tm_hour;
      day = tm_p->tm_mday;
      month = tm_p->tm_mon + 1;
      year = tm_p->tm_year + 1900;

      /* ACE data starts at year ESTIST2_FIRST_YEAR */
      dptr = &(w->estist[year - ESTIST2_FIRST_YEAR]);

      assert(month >= 1 && month <= 12);
      assert(day >= 1 && day <= 31);

      dptr->E_st[month - 1][day - 1][hh] = E_st;
      dptr->I_st[month - 1][day - 1][hh] = I_st;
    }

  fclose(fp);

  return s;
} /* estist2_read_data() */

/*
estist2_alloc()
  Allocate a ace workspace and read in data from filename

Inputs: idx_filename - index of ESTIST2 files

Return: pointer to workspace
*/

estist2_workspace *
estist2_alloc(const char *idx_filename)
{
  estist2_workspace *w;
  size_t i, j, k, l;
  estist2_data *dptr;
  int s;

  w = (estist2_workspace *) calloc(1, sizeof(estist2_workspace));
  if (!w)
    return 0;

  for (i = 0; i < ESTIST2_MAX_YEAR; ++i)
    {
      dptr = &(w->estist[i]);
      for (j = 0; j < 12; ++j)
        for (k = 0; k < 31; ++k)
          for (l = 0; l < 24; ++l)
            {
              dptr->E_st[j][k][l] = ESTIST2_BAD_VALUE;
              dptr->I_st[j][k][l] = ESTIST2_BAD_VALUE;
            }
    }

  /* read in available data */
  s = estist2_read_data(idx_filename, w);
  if (s)
    {
      estist2_free(w);
      return 0;
    }

  return (w);
} /* estist2_alloc() */

void
estist2_free(estist2_workspace *w)
{
  if (!w)
    return;

  free(w);
} /* estist2_free() */

/*
estist2_get()
  Obtain the Est/Ist value for a given timestamp

Inputs: t      - UT timestamp
        E_st   - (output) E_st in nT
        I_st   - (output) I_st in nT
        w      - ace workspace

Return: 0 on success, 1 on failure
*/

int
estist2_get(time_t t, double *E_st, double *I_st, estist2_workspace *w)
{
  int s = 0;
  int year;
  estist2_data *dptr;
  struct tm *tm_p;
  double ut;
  size_t k;

  tm_p = gmtime(&t);

  year = tm_p->tm_year + 1900;

  if (year < ESTIST2_FIRST_YEAR || year > ESTIST2_FIRST_YEAR + ESTIST2_MAX_YEAR)
    {
      fprintf(stderr, "estist2_get: invalid year: %d\n", year);
      return 1;
    }

  ut = (double) tm_p->tm_hour +
       (double) tm_p->tm_min / 60.0 +
       (double) tm_p->tm_sec / 3600.0;

  dptr = &(w->estist[year - ESTIST2_FIRST_YEAR]);

  k = (size_t) ut;
  assert(k < 24);

  *E_st = dptr->E_st[tm_p->tm_mon][tm_p->tm_mday - 1][k];
  *I_st = dptr->I_st[tm_p->tm_mon][tm_p->tm_mday - 1][k];

  if (*E_st == ESTIST2_BAD_VALUE)
    {
      fprintf(stderr, "estist2_get: unknown E_st value for timestamp %ld\n", t);
    }

  return s;
} /* estist2_get() */
