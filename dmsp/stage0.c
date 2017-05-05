/*
 * stage0.c
 *
 * Convert ASCII DMSP data files into CDF format with possible fixing
 * of ephemeris values
 *
 * Usage: ./stage0 <-i input_ascii_gz_file> [-o output_cdf_file]
 *                 [-n nasa_ephemeris_file] [-s sgp4_ephemeris_file]
 *                 [-b bowman_ephemeris_file]
 *
 * The measurements are read from the ascii file, ephemeris values are
 * possibly modified (lat,lon,alt) if needed, and output is written to CDF
 *
 * Altitudes are stored in CDF file using reference radius:
 *
 * R_EARTH_KM = 6371.2 km
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <sys/time.h>
#include <zlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>

#include <indices/indices.h>
#include <satdata/satdata.h>

#include "common.h"
#include "eci.h"
#include "eph.h"
#include "eph_data.h"
#include "hermite.h"

/* read DMSP SSM Fred Rich ASCII data file */
size_t
dmsp_read_MFR(const char *filename, satdata_mag *data)
{
  const double R = R_EARTH_KM;
  char buf[SATDATA_MAX_BUF];
  size_t n = 0;
  gzFile fp;

  fp = gzopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "unable to open %s: %s\n", filename, strerror(errno));
      return 0;
    }

  while (gzgets(fp, buf, SATDATA_MAX_BUF) != 0)
    {
      int c;
      double date, sec, lat, lon, alt, X, Y, Z;
      char as[10], istr[10], dstr[10];
      int mid;
      int year, month, day, doy;
      int hh, mm, ss, msec;

      if (*buf == '#')
        continue;

      c = sscanf(buf, "%lf %lf %d %s %lf %lf %lf %s %s %lf %lf %lf",
                 &date,
                 &sec,
                 &mid,
                 as,
                 &lat,
                 &lon,
                 &alt,
                 istr,
                 dstr,
                 &X,
                 &Y,
                 &Z);
      if (c < 12)
        continue;

      year = (int) (date / 1.0e9);
      doy = (int) ((date - year * 1.0e9) / 1.0e6);

      hh = (int) (sec / 3600.0);
      mm = (int) ((sec - hh * 3600.0) / 60.0);
      ss = (int) (sec - hh * 3600.0 - mm * 60.0);
      msec = (int) ((sec - (int) sec) * 1000.0);

      doy2md(year, doy, &month, &day);

      data->t[n] = computeEPOCH(year, month, day, hh, mm, ss, msec);
      data->latitude[n] = lat;
      data->longitude[n] = lon;
      data->r[n] = R + alt;

      /*
       * the (X,Y,Z) in the DMSP ASCII files are not NEC so they are
       * stored differently below
       */
      SATDATA_VEC_X(data->B_VFM, n) = Y; /* velocity direction */
      SATDATA_VEC_Y(data->B_VFM, n) = Z; /* orbit normal */
      SATDATA_VEC_Z(data->B_VFM, n) = X; /* down */

      data->F[n] = gsl_hypot3(X, Y, Z);

      if (++n >= data->ntot)
        {
          fprintf(stderr, "dmsp_read_MFR: file %s contains too many data records\n",
                  filename);
          return n;
        }
    }

  gzclose(fp);

  data->n = n;
  data->R = R;

  return n;
}

int
interp_eph(satdata_mag *data, eph_data *eph)
{
  int s = 0;
  size_t i;
  eph_workspace *w = eph_alloc(eph);

  for (i = 0; i < data->n; ++i)
    {
      double r_sph[3]; /* position in spherical coordinates */
      double r, lat, lon;

      /* interpolate ephemeris data to time ti */
      s = eph_interp_sph(data->t[i], r_sph, w);
      if (s)
        {
          data->flags[i] |= SATDATA_FLG_NOEPH;
          continue;
        }

      /* convert to degrees */
      r = r_sph[0];
      lat = 90.0 - r_sph[1] * 180.0 / M_PI;
      lon = r_sph[2] * 180.0 / M_PI;

      data->r[i] = r;
      data->latitude[i] = lat;
      data->longitude[i] = wrap180(lon);
    }

  eph_free(w);

  return 0;
}

int
main(int argc, char *argv[])
{
  char *infile = NULL;
  char *outfile = NULL;
  satdata_mag *data = NULL, *data_final;
  eph_data *eph = NULL;
  int c;
  struct timeval tv0, tv1;

  while ((c = getopt(argc, argv, "i:o:b:t:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'b':
            fprintf(stderr, "main: reading Bowman ephemerides from %s...", optarg);
            gettimeofday(&tv0, NULL);
            eph = eph_data_read_bowman(optarg);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu read, %g seconds)\n", eph->n, time_diff(tv0, tv1));
            break;

          case 't':
            fprintf(stderr, "main: reading TENA ephemerides from %s...", optarg);
            gettimeofday(&tv0, NULL);
            eph = eph_data_read_tena(optarg);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu read, %g seconds)\n", eph->n, time_diff(tv0, tv1));
            break;

          case 'o':
            outfile = optarg;
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i DMSP_ascii_gz_file> [-o output_cdf_file] [-b bowman_ephemeris_file] [-t tena_ephemeris_file]\n",
              argv[0]);
      exit(1);
    }

  data = satdata_mag_alloc(86400);
  data_final = satdata_mag_alloc(86400);

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);
  dmsp_read_MFR(infile, data);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
                   time_diff(tv0, tv1));

  if (data->n > 0)
    {
      if (eph)
        {
          fprintf(stderr, "main: interpolating ephemeris of %s...", infile);
          interp_eph(data, eph);
          fprintf(stderr, "done\n");
        }

      fprintf(stderr, "main: copying data...");
      satdata_select_filter(data, data_final);
      fprintf(stderr, "done (%zu data thrown out due to no ephemeris)\n",
              data->n - data_final->n);
 
      if (outfile && data_final->n > 0)
        {
          fprintf(stderr, "main: writing %s...", outfile);
          gettimeofday(&tv0, NULL);
          satdata_dmsp_write(0, outfile, data_final);
          gettimeofday(&tv1, NULL);
          fprintf(stderr, "done (%zu records written, %g seconds)\n", data_final->n,
                  time_diff(tv0, tv1));
        }
    }

  satdata_mag_free(data);
  satdata_mag_free(data_final);
  if (eph)
    eph_data_free(eph);

  return 0;
}
