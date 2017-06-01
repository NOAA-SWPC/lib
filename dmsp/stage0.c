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
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <sys/time.h>
#include <zlib.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_test.h>

#include <indices/indices.h>
#include <satdata/satdata.h>

#include "common.h"
#include "eci.h"
#include "eph.h"
#include "eph_data.h"
#include "ellipsoid.h"
#include "hermite.h"
#include "quat.h"

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

/*
calc_spacecraft_basis()
  Calculate spacecraft basis to rotate a vector from the spacecraft (S/C) frame
to NEC. We define spacecraft-fixed basis vectors as

s1 = s2 x s3 (velocity direction)
s2 = (s3 x v) / | s3 x v |
s3 = -rhat (geocentric downward)

Note: While the DMSP attitude control is supposed to keep the satellite fixed
wrt the geodetic normal, the MFR data files provide the geocentric vertical
component, so it appears we do not need to compute the geodetic direction ourselves,
and can just use rhat

Inputs: r_ECEF - position (X,Y,Z) in Cartesian ECEF (km)
        v_ECEF - velocity (VX,VY,VZ) in Cartesian ECEF (km/s)
        s1     - (output) s1 unit vector (ECEF)
        s2     - (output) s2 unit vector (ECEF)
        s3     - (output) s3 unit vector (ECEF)
*/

static int
calc_spacecraft_basis(const double r_ECEF[3], const double v_ECEF[3],
                      double s1[3], double s2[3], double s3[3])
{
  int s = 0;
  size_t i;
  double norm;

#if 1 /* use s3 = -e_mu */

  /* compute ECEF components of ellipsoid basis vectors, storing e_mu in s3 */
  /*ellipsoid_basis_mu(r_ECEF, WGS84_MU, s3, s1, s2);*/
  ellipsoid_basis(r_ECEF, s3, s1, s2);

#else /* use s3 = -rhat */

  /* store rhat in s3 */
  ecef2sph_basis(r_ECEF, s3, s1, s2);

#endif

  /* reverse s3 to point downward */
  for (i = 0; i < 3; ++i)
    s3[i] *= -1.0;

  /* s2 = (s3 x v) / | s3 x v | */
  vec_cross(s3, v_ECEF, s2);

  norm = vec_norm(s2);
  for (i = 0; i < 3; ++i)
    s2[i] /= norm;

  /* s1 = s2 x s3 */
  vec_cross(s2, s3, s1);

  return s;
}

/*
calc_quaternions_ECEF()
  Calculate quaternions to rotate a vector from the spacecraft (S/C) frame
to NEC.

Inputs: r_ECEF - position (X,Y,Z) in Cartesian ECEF (km)
        v_ECEF - velocity (VX,VY,VZ) in Cartesian ECEF (km/s)
        q      - (output) quaternions for rotation
*/

int
calc_quaternions_ECEF(const double r_ECEF[3], const double v_ECEF[3], double q[4])
{
  int s = 0;
  size_t i;
  double s1[3], s2[3], s3[3];        /* spacecraft-fixed unit basis vectors (ECEF) */
  double rhat[3], that[3], phat[3];  /* spherical unit basis vectors (ECEF) */
  double nhat[3], ehat[3], chat[3];  /* NEC unit basis vectors (ECEF) */
  double R_data[9];
  gsl_matrix_view R = gsl_matrix_view_array(R_data, 3, 3);

  /* compute spacecraft-fixed basis vectors in ECEF */
  calc_spacecraft_basis(r_ECEF, v_ECEF, s1, s2, s3);

  /* compute spherical basis vectors in ECEF */
  ecef2sph_basis(r_ECEF, rhat, that, phat);

  /* convert to NEC vectors */
  for (i = 0; i < 3; ++i)
    {
      nhat[i] = -that[i];
      ehat[i] = phat[i];
      chat[i] = -rhat[i];
    }

  /* build rotation matrix from S/C to NEC */

  gsl_matrix_set(&R.matrix, 0, 0, vec_dot(s1, nhat));
  gsl_matrix_set(&R.matrix, 0, 1, vec_dot(s2, nhat));
  gsl_matrix_set(&R.matrix, 0, 2, vec_dot(s3, nhat));

  gsl_matrix_set(&R.matrix, 1, 0, vec_dot(s1, ehat));
  gsl_matrix_set(&R.matrix, 1, 1, vec_dot(s2, ehat));
  gsl_matrix_set(&R.matrix, 1, 2, vec_dot(s3, ehat));

  gsl_matrix_set(&R.matrix, 2, 0, vec_dot(s1, chat));
  gsl_matrix_set(&R.matrix, 2, 1, vec_dot(s2, chat));
  gsl_matrix_set(&R.matrix, 2, 2, vec_dot(s3, chat));

  /* convert to quaternions */
  quat_R2q(&R.matrix, q);

  /* XXX sanity check */
  {
    double Rq_data[9];
    gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
    gsl_vector_view v1 = gsl_vector_view_array(R_data, 9);
    gsl_vector_view v2 = gsl_vector_view_array(Rq_data, 9);
    double norm;

    quat_q2R(q, &Rq.matrix);

    gsl_vector_sub(&v1.vector, &v2.vector);
    norm = gsl_blas_dnrm2(&v1.vector);

    if (norm > 10.0 * GSL_DBL_EPSILON)
      fprintf(stderr, "error: || R - Rq || = %.12e\n", gsl_blas_dnrm2(&v1.vector));
  }

  return s;
}

int
interp_eph(satdata_mag *data, eph_data *eph)
{
  int s = 0;
  size_t i, j;
  eph_workspace *w = eph_alloc(eph);

  for (i = 0; i < data->n; ++i)
    {
      double pos[3], vel[3]; /* position and velocity (ECI or ECEF) */
      double r, theta, phi;
      double q[4];           /* quaternions for rotation S/C to NEC */

      /* interpolate ephemeris data to time ti */
      s = eph_interp(data->t[i], pos, vel, w);
      if (s)
        {
          data->flags[i] |= SATDATA_FLG_NOEPH;
          continue;
        }

      if (w->data->flags & EPH_DATA_FLG_ECEF)
        {
          /* compute (r,theta,phi) for this point from ECEF position */
          r = gsl_hypot3(pos[0], pos[1], pos[2]);
          theta = acos(pos[2] / r);
          phi = atan2(pos[1], pos[0]);

          calc_quaternions_ECEF(pos, vel, q);
        }
      else
        {
          fprintf(stderr, "interp_eph: ECI not yet supported\n");
          return -1;
        }

      data->r[i] = r;
      data->latitude[i] = 90.0 - theta * 180.0 / M_PI;
      data->longitude[i] = wrap180(phi * 180.0 / M_PI);

      for (j = 0; j < 4; ++j)
        data->q[4*i + j] = q[j];
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
