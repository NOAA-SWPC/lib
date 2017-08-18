/*
 * iaga.c
 *
 * Routines for reading IAGA formatted data files
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <strings.h>
#include <ctype.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "grobs.h"

static int iaga_read_HDZF(FILE *fp, grobs_data *data);
static int iaga_read_XYZF(FILE *fp, grobs_data *data);

grobs_data *
grobs_iaga_read(const char *filename, grobs_data *data)
{
  FILE *fp;
  char buf[GROBS_MAX_BUFFER];
  char s1[GROBS_MAX_BUFFER];
  char s2[GROBS_MAX_BUFFER];
  int header_cnt = 0;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "fopen: cannot open %s: %s\n",
              filename, strerror(errno));
      return NULL;
    }

  /* allocate space for 1-minute data */
  if (data == NULL)
    data = grobs_alloc(GROBS_MAX_YEAR * 527040);

  while (fgets(buf, GROBS_MAX_BUFFER, fp) != NULL)
    {
      int c;
      double val;
      char *bufptr = buf;

      while (isspace(*bufptr))
        bufptr++;

      if (*bufptr == '#')
        continue;

      if (!strncasecmp(bufptr, "geodetic latitude", 17))
        {
          c = sscanf(bufptr, "%s %s %lf", s1, s2, &val);
          if (c == 3)
            {
              ++header_cnt;
              data->glat = val;
            }
        }
      else if (!strncasecmp(bufptr, "geodetic longitude", 18))
        {
          c = sscanf(bufptr, "%s %s %lf", s1, s2, &val);
          if (c == 3)
            {
              ++header_cnt;
              data->glon = val;
            }
        }
      else if (!strncasecmp(bufptr, "reported", 8))
        {
          c = sscanf(bufptr, "%s %s", s1, s2);
          if (c == 2)
            {
              ++header_cnt;
              if (!strncmp(s2, "HDZF", 4))
                {
                  iaga_read_HDZF(fp, data);
                }
              else if (!strncmp(s2, "XYZF", 4) || !strncmp(s2, "XYZG", 4))
                {
                  iaga_read_XYZF(fp, data);
                }
              else
                {
                  fprintf(stderr, "grobs_iaga_read: unknown format type: %s\n", s2);
                  header_cnt = 0;
                }
            }

          /* this is the last header line read, so break after
           * reading the file */
          break;
        }
    }

  fclose(fp);

  if (header_cnt != 3)
    {
      fprintf(stderr, "grobs_iaga_read: failed to read header\n");
    }

  return data;
}

static int
iaga_read_HDZF(FILE *fp, grobs_data *data)
{
  int ret;
  char buf[GROBS_MAX_BUFFER];
  size_t n = data->n;

  while (fgets(buf, GROBS_MAX_BUFFER, fp) != NULL)
    {
      int year, month, day;
      int hour, min, doy;
      double fsec;
      double Z, D, F, H, Drad;

      /* read day number */
      ret = sscanf(buf, "%d-%d-%d %d:%d:%lf %d %lf %lf %lf %lf",
                   &year, &month, &day, &hour, &min, &fsec,
                   &doy, &H, &D, &Z, &F);
      if (ret < 11)
        continue;

      if (day < 1 || day > 31)
        continue;

      if (month < 1 || month > 12)
        continue;

      /* check missing data */
      if (fabs(H) > 90000.0)
        continue;
      if (fabs(D) > 90000.0)
        continue;
      if (fabs(Z) > 90000.0)
        continue;
      if (fabs(F) > 90000.0)
        continue;

      /* D is in arcminutes */
      Drad = (D / 60.0) * M_PI / 180.0;

      /* store value */
      data->t[n] = date2timet((int) fsec, min, hour, day, month, year);
      data->X[n] = H * cos(Drad);
      data->Y[n] = H * sin(Drad);
      data->Z[n] = Z;
      data->H[n] = H;
      data->D[n] = Drad * 180.0 / M_PI;
      data->I[n] = atan2(Z, H);

      if (++n >= data->ntot)
        {
          fprintf(stderr, "iaga_read: not enough space allocated: %zu\n", data->ntot);
          break;
        }
    }

  data->n = n;

  return 0;
}

static int
iaga_read_XYZF(FILE *fp, grobs_data *data)
{
  int ret;
  char buf[GROBS_MAX_BUFFER];
  size_t n = data->n;

  while (fgets(buf, GROBS_MAX_BUFFER, fp) != NULL)
    {
      int year, month, day;
      int hour, min, doy;
      double fsec;
      double X, Y, Z, F, H;

      /* read day number */
      ret = sscanf(buf, "%d-%d-%d %d:%d:%lf %d %lf %lf %lf %lf",
                   &year, &month, &day, &hour, &min, &fsec,
                   &doy, &X, &Y, &Z, &F);
      if (ret < 11)
        continue;

      if (day < 1 || day > 31)
        continue;

      if (month < 1 || month > 12)
        continue;

      /* check missing data */
      if (fabs(X) > 90000.0)
        continue;
      if (fabs(Y) > 90000.0)
        continue;
      if (fabs(Z) > 90000.0)
        continue;
      if (fabs(F) > 90000.0)
        continue;

      H = gsl_hypot(X, Y);

      /* store value */
      data->t[n] = date2timet((int) fsec, min, hour, day, month, year);
      data->X[n] = X;
      data->Y[n] = Y;
      data->Z[n] = Z;
      data->H[n] = H;
      data->D[n] = atan2(Y, X) * 180.0 / M_PI;
      data->I[n] = atan2(Z, H);

      if (++n >= data->ntot)
        {
          fprintf(stderr, "iaga_read: not enough space allocated: %zu\n", data->ntot);
          break;
        }
    }

  data->n = n;

  return 0;
}
