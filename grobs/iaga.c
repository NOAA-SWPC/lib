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

#include <gsl/gsl_math.h>

#include "common.h"
#include "grobs.h"

grobs_data *
iaga_read_HDZF(const char *filename, grobs_data *data)
{
  FILE *fp;
  int ret;
  char buf[GROBS_MAX_BUFFER];
  size_t n;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "fopen: cannot open %s: %s\n",
              filename, strerror(errno));
      return NULL;
    }

  if (data == NULL)
    {
      data = calloc(1, sizeof(grobs_data));
      n = 0;
    }
  else
    {
      n = data->n;
    }

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
      data->D[n] = Drad * 180.0 / M_PI;
      data->I[n] = atan2(Z, H);

      ++n;
    }

  fclose(fp);

  data->n = n;

  return data;
}

grobs_data *
iaga_read_XYZF(const char *filename, grobs_data *data)
{
  FILE *fp;
  int ret;
  char buf[GROBS_MAX_BUFFER];
  size_t n;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "fopen: cannot open %s: %s\n",
              filename, strerror(errno));
      return NULL;
    }

  if (data == NULL)
    {
      data = calloc(1, sizeof(grobs_data));
      n = 0;
    }
  else
    {
      n = data->n;
    }

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
      data->D[n] = atan2(Y, X) * 180.0 / M_PI;
      data->I[n] = atan2(Z, H);

      ++n;
    }

  fclose(fp);

  data->n = n;

  return data;
}
