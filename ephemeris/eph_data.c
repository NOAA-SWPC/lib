/*
 * eph_data.c
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <cdf.h>
#include <zlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "common.h"
#include "eph_data.h"

void
eph_data_free(eph_data *data)
{
  free(data);
}

/* read ephemeris (gzip) data file as provided by Bruce Bowman */
eph_data *
eph_data_read_bowman(const char *filename)
{
  eph_data *data;
  gzFile fp;
  char buffer[2048];
  size_t n = 0;

  fp = gzopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "eph_data_read_bowman: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  data = malloc(sizeof(eph_data));
  data->flags = EPH_DATA_FLG_ECI;

  while (gzgets(fp, buffer, 2048) != 0)
    {
      int c;
      int year, doy;
      int month, day, hour, min, sec, msec;
      double r_ECI[3], latc, lonc, fsec, alt, lt, rho;
      double VX, VY, VZ;

      if (*buffer == '#')
        continue;

      c = sscanf(buffer, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                 &year,
                 &doy,
                 &hour,
                 &min,
                 &fsec,
                 &alt,
                 &lt,
                 &latc,
                 &lonc,
                 &rho,
                 &r_ECI[0],
                 &r_ECI[1],
                 &r_ECI[2],
                 &VX,
                 &VY,
                 &VZ);
      if (c < 16)
        continue;

      sec = (int) fsec;
      msec = (int) ((fsec - sec) * 1000.0);

      doy2md(year, doy, &month, &day);

      data->t[n] = computeEPOCH(year, month, day, hour, min, sec, msec);
      data->X[n] = r_ECI[0];
      data->Y[n] = r_ECI[1];
      data->Z[n] = r_ECI[2];
      data->VX[n] = VX;
      data->VY[n] = VY;
      data->VZ[n] = VZ;
      data->latitude[n] = latc;
      data->longitude[n] = lonc;

      if (++n >= EPH_MAX_DATA)
        {
          fprintf(stderr, "EPH_MAX_DATA too small\n");
          break;
        }
    }

  gzclose(fp);

  data->n = n;

  return data;
} /* eph_data_read_bowman() */

/* read ephemeris data file as provided by TENA */
eph_data *
eph_data_read_tena(const char *filename)
{
  eph_data *data;
  FILE *fp;
  char buffer[2048];
  size_t n = 0;
  int gotheader = 0;
  int year = 0, doy = 0;
  char *ptr;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "eph_data_read_tena: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  data = malloc(sizeof(eph_data));
  data->flags = EPH_DATA_FLG_ECEF;

  while (fgets(buffer, 2048, fp) != 0)
    {
      int c;
      int dummy;
      int month, day, hour, min, sec, msec;
      double fsec;
      double X, Y, Z;
      double VX, VY, VZ;

      /* headers start as the first character on the line; there
       * can be multiple headers per file */
      if (*buffer == '2')
        {
          c = sscanf(buffer, "%d %d %d %d\n",
                     &year,
                     &doy,
                     &dummy,
                     &dummy);

          if (c == 4)
            {
              doy2md(year, doy, &month, &day);
              gotheader = 1;
            }

          continue;
        }

      if (!gotheader)
        continue;

      c = sscanf(buffer, "%lf %lf %lf %lf\n",
                 &fsec,
                 &X,
                 &Y,
                 &Z);
      if (c < 4)
        {
          fprintf(stderr, "eph_data_read_tena: error reading position\n");
          continue;
        }

      /* read the next line (with velocities) */

      ptr = fgets(buffer, 2048, fp);
      if (ptr == NULL)
        {
          fprintf(stderr, "eph_data_read_tena: error reading velocity line\n");
          continue;
        }

      c = sscanf(buffer, "%lf %lf %lf %lf\n",
                 &fsec,
                 &VX,
                 &VY,
                 &VZ);
      if (c < 4)
        {
          fprintf(stderr, "eph_data_read_tena: error reading velocity\n");
          continue;
        }

      hour = (int) (fsec / 3600.0);
      fsec -= hour * 3600;
      min = (int) (fsec / 60.0);
      fsec -= min * 60;
      sec = (int) fsec;
      msec = 0;

      data->t[n] = computeEPOCH(year, month, day, hour, min, sec, msec);
      data->X[n] = X;
      data->Y[n] = Y;
      data->Z[n] = Z;
      data->VX[n] = VX;
      data->VY[n] = VY;
      data->VZ[n] = VZ;

      if (++n >= EPH_MAX_DATA)
        {
          fprintf(stderr, "EPH_MAX_DATA too small\n");
          break;
        }
    }

  fclose(fp);

  data->n = n;

  return data;
}
