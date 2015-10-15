/*
 * omni.c
 * Patrick Alken
 *
 * This module reads data from a OMNI2 data file containing solar wind
 * parameters
 *
 * The data files are taken from:
 *
 * ftp://nssdcftp.gsfc.nasa.gov/spacecraft_data/omni
 *
 * The file format is described in the omni2.text file in that dir
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>

#include "common.h"
#include "omni.h"

/*
omni_read_data()
  Read OMNI data from a Penticton data file

Inputs: filename - OMNI data file
        w        - omni workspace
*/

static int
omni_read_data(const char *filename, omni_workspace *w)
{
  int s = 0;
  FILE *fp;
  int ret;
  char buf[OMNI_MAX_BUFFER];
  omni_data *optr;
  int year;
  double Bx, By, Bz, V;
  int hour, doy;
  int dummy;
  double dummyf;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "fopen: cannot open %s: %s\n",
              filename, strerror(errno));
      return 1;
    }

  while (fgets(buf, OMNI_MAX_BUFFER, fp) != NULL)
    {
      if (*buf == '#')
        continue;

      ret = sscanf(buf, "%d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   &year,
                   &doy,
                   &hour,
                   &dummy,
                   &dummy,
                   &dummy,
                   &dummy,
                   &dummy,
                   &dummyf,
                   &dummyf,
                   &dummyf,
                   &dummyf,
                   &Bx,
                   &By,
                   &Bz,
                   &dummyf,
                   &dummyf,
                   &dummyf,
                   &dummyf,
                   &dummyf,
                   &dummyf,
                   &dummyf,
                   &dummyf,
                   &dummyf,
                   &V);
      if (ret < 25)
        continue;

      assert(doy >= 1 && doy <= 366);
      assert(hour >= 0 && hour <= 23);

      if (year < OMNI_FIRST_YEAR)
        {
          fprintf(stderr, "omni_read_data: invalid year: %d\n", year);
          continue;
        }

      optr = &(w->data[year - OMNI_FIRST_YEAR][doy - 1][hour]);

      /* store value */
      optr->Bx = Bx;
      optr->By = By;
      optr->Bz = Bz;
      optr->V = V;
    }

  fclose(fp);

  return s;
} /* omni_read_data() */

/*
omni_alloc()
  Allocate a omni workspace and read in data from filename

Inputs: datadir - OMNI data directory

Return: pointer to workspace
*/

omni_workspace *
omni_alloc(const char *datadir)
{
  omni_workspace *w;
  size_t i, j, k;
  int year;
  char filename[OMNI_MAX_BUFFER];

  w = (omni_workspace *) calloc(1, sizeof(omni_workspace));
  if (!w)
    return 0;

  for (i = 0; i < OMNI_MAX_YEAR; ++i)
    {
      for (j = 0; j < 366; ++j)
        {
          for (k = 0; k < 24; ++k)
            {
              omni_data *optr = &(w->data[i][j][k]);
              optr->Bx = OMNI_MISSING;
              optr->By = OMNI_MISSING;
              optr->Bz = OMNI_MISSING;
              optr->V = OMNI_MISSING;
            }
        }
    }

  for (year = 1980; year < 2012; ++year)
    {
      sprintf(filename, "%s/omni2_%04d.dat", datadir, year);
      omni_read_data(filename, w);
    }

  return (w);
} /* omni_alloc() */

void
omni_free(omni_workspace *w)
{
  free(w);
} /* omni_free() */

/*
omni_get()
  Obtain the F10.7 value for a given timestamp

Inputs: t      - UT timestamp
        B      - (output) where to store IMF field in GSE coordinates (nT)
        V      - (output) where to store wind velocity (km/s)
        w      - omni workspace

Return: 0 on success, 1 on failure
*/

int
omni_get(time_t t, double B[4], double *V, omni_workspace *w)
{
  int year, doy, hour;
  struct tm *tm_p;
  omni_data *optr;

  tm_p = gmtime(&t);
  year = tm_p->tm_year + 1900;
  doy = tm_p->tm_yday;
  hour = tm_p->tm_hour;

  optr = &(w->data[year - OMNI_FIRST_YEAR][doy][hour]);

  if (optr->Bx == OMNI_MISSING)
    {
      fprintf(stderr, "omni_get: data not found for timestamp %ld: %s",
              t, asctime(tm_p));
      return 1;
    }

  B[0] = optr->Bx;
  B[1] = optr->By;
  B[2] = optr->Bz;
  B[3] = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

  *V = optr->V;

  return 0;
} /* omni_get() */
