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

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_blas.h>

#include <indices/indices.h>
#include <satdata/satdata.h>

#include "common.h"
#include "eci.h"
#include "eph.h"
#include "eph_data.h"
#include "ellipsoid.h"
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

/*
calc_spacecraft_basis()
  Calculate spacecraft basis to rotate a vector from the spacecraft (S/C) frame
to NEC. We define spacecraft-fixed basis vectors as

s1 = v_t / |v_t|
s2 = s3 x s1
s3 = -rhat (geocentric downward)

where

v_t = v - (v . s3) s3 (component of velocity normal to s3)

Note: While the DMSP attitude control is supposed to keep the satellite fixed
wrt the geodetic normal, the MFR data files provide the geocentric vertical
component, so it appears we do not need to compute the geodetic direction ourselves,
and can just use rhat

Inputs: theta  - colatitude (radians)
        phi    - longitude (radians)
        r_ECEF - position (X,Y,Z) in Cartesian ECEF (km)
        v_ECEF - velocity (VX,VY,VZ) in Cartesian ECEF (km/s)
        s1     - (output) s1 unit vector (ECEF)
        s2     - (output) s2 unit vector (ECEF)
        s3     - (output) s3 unit vector (ECEF)
*/

static int
calc_spacecraft_basis(const double theta, const double phi,
                      const double r_ECEF[3], const double v_ECEF[3],
                      double s1[3], double s2[3], double s3[3])
{
  int s = 0;
  size_t i;
  double vt[3]; /* component of velocity normal to s3 */
  double v3;    /* component of v_ECEF along s3 */

#if 0 /* use s3 = -e_mu */

  /* compute ECEF components of ellipsoid basis vectors, storing e_mu in s3 */
  ellipsoid_basis_mu(r_ECEF, WGS84_MU, s3, s1, s2);

#else /* use s3 = -rhat */

  /* store rhat in s3 */
  sph_basis(theta, phi, s3, s1, s2);

#endif

  /* reverse s3 to point downward */
  for (i = 0; i < 3; ++i)
    s3[i] *= -1.0;

  /* compute component of velocity along s3 */
  v3 = vec_dot(v_ECEF, s3);

  /* compute v_t = v - v3 * s3 */
  for (i = 0; i < 3; ++i)
    vt[i] = v_ECEF[i] - v3 * s3[i];

  /* s1 = v_t / |v_t| */
  for (i = 0; i < 3; ++i)
    s1[i] = vt[i] / vec_norm(vt);

  /* s2 = s3 x s1 */
  sphcross(s3, s1, s2);

  return s;
}

static int
euler_Rq(const double *q, gsl_matrix *Rq)
{
  const double q1 = q[0];
  const double q2 = q[1];
  const double q3 = q[2];
  const double q4 = q[3];

  gsl_matrix_set(Rq, 0, 0, 1.0 - 2.0*q2*q2 - 2.0*q3*q3);
  gsl_matrix_set(Rq, 0, 1, 2.0*(q1*q2 + q3*q4));
  gsl_matrix_set(Rq, 0, 2, 2.0*(q1*q3 - q2*q4));

  gsl_matrix_set(Rq, 1, 0, 2.0*(q1*q2 - q3*q4));
  gsl_matrix_set(Rq, 1, 1, 1.0 - 2.0*q1*q1 - 2.0*q3*q3);
  gsl_matrix_set(Rq, 1, 2, 2.0*(q2*q3 + q1*q4));

  gsl_matrix_set(Rq, 2, 0, 2.0*(q1*q3 + q2*q4));
  gsl_matrix_set(Rq, 2, 1, 2.0*(q2*q3 - q1*q4));
  gsl_matrix_set(Rq, 2, 2, 1.0 - 2.0*q1*q1 - 2.0*q2*q2);

  return GSL_SUCCESS;
}

/*
calc_quaternions_ECEF()
  Calculate quaternions to rotate a vector from the spacecraft (S/C) frame
to NEC.

Inputs: theta  - colatitude (radians)
        phi    - longitude (radians)
        r_ECEF - position (X,Y,Z) in Cartesian ECEF (km)
        v_ECEF - velocity (VX,VY,VZ) in Cartesian ECEF (km/s)
        q      - (output) quaternions for rotation
*/

int
calc_quaternions_ECEF(const double theta, const double phi, const double r_ECEF[3], const double v_ECEF[3], double q[4])
{
  int s = 0;
  size_t i;
  double s1[3], s2[3], s3[3];        /* spacecraft-fixed unit basis vectors (ECEF) */
  double rhat[3], that[3], phat[3];  /* spherical unit basis vectors (ECEF) */
  double nhat[3], ehat[3], chat[3];  /* NEC unit basis vectors (ECEF) */
  double R11, R12, R13;
  double R21, R22, R23;
  double R31, R32, R33;
  double Trace, A;

  /* compute spacecraft-fixed basis vectors in ECEF */
  calc_spacecraft_basis(theta, phi, r_ECEF, v_ECEF, s1, s2, s3);

  /* compute spherical basis vectors in ECEF */
  sph_basis(theta, phi, rhat, that, phat);

  /* convert to NEC vectors */
  for (i = 0; i < 3; ++i)
    {
      nhat[i] = -that[i];
      ehat[i] = phat[i];
      chat[i] = -rhat[i];
    }

  /* build rotation matrix from S/C to NEC */

  R11 = vec_dot(s1, nhat);
  R12 = vec_dot(s2, nhat);
  R13 = vec_dot(s3, nhat);

  R21 = vec_dot(s1, ehat);
  R22 = vec_dot(s2, ehat);
  R23 = vec_dot(s3, ehat);

  R31 = vec_dot(s1, chat);
  R32 = vec_dot(s2, chat);
  R33 = vec_dot(s3, chat);

  Trace = R11 + R22 + R33;
  A = 0.25 * (1.0 - Trace);

  q[0] = 0.5 * R11 + A;
  q[1] = 0.5 * R22 + A;
  q[2] = 0.5 * R33 + A;
  q[3] = 0.25 * (1.0 + Trace);

  if ((q[0] >= q[1]) && (q[0] >= q[2]) && (q[0] >= q[3]))
    {
      q[0] = sqrt(q[0]);
      q[1] = (R21 + R12) / (4.0 * q[0]);
      q[2] = (R31 + R13) / (4.0 * q[0]);
      q[3] = (R23 - R32) / (4.0 * q[0]);
    }

  if ((q[1] >= q[0]) && (q[1] >= q[2]) && (q[1] >= q[3]))
    {
      q[1] = sqrt(q[1]);
      q[2] = (R32 + R23) / (4.0 * q[1]);
      q[3] = (R31 - R13) / (4.0 * q[1]);
      q[0] = (R12 + R21) / (4.0 * q[1]);
    }

  if ((q[2] >= q[0]) && (q[2] >= q[1]) && (q[2] >= q[3]))
    {
      q[2] = sqrt(q[2]);
      q[3] = (R12 - R21) / (4.0 * q[2]);
      q[0] = (R13 + R31) / (4.0 * q[2]);
      q[1] = (R23 + R32) / (4.0 * q[2]);
    }

  if ((q[3] >= q[0]) && (q[3] >= q[1]) && (q[3] >= q[2]))
    {
      q[3] = sqrt(q[3]);
      q[0] = (R23 - R32) / (4.0 * q[3]);
      q[1] = (R31 - R13) / (4.0 * q[3]);
      q[2] = (R12 - R21) / (4.0 * q[3]);
    }

  {
    double R_data[9], Rq_data[9];
    gsl_matrix_view R = gsl_matrix_view_array(R_data, 3, 3);
    gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
    gsl_vector_view v1 = gsl_vector_view_array(R_data, 9);
    gsl_vector_view v2 = gsl_vector_view_array(Rq_data, 9);
    double norm;

    gsl_matrix_set(&R.matrix, 0, 0, R11);
    gsl_matrix_set(&R.matrix, 0, 1, R12);
    gsl_matrix_set(&R.matrix, 0, 2, R13);

    gsl_matrix_set(&R.matrix, 1, 0, R21);
    gsl_matrix_set(&R.matrix, 1, 1, R22);
    gsl_matrix_set(&R.matrix, 1, 2, R23);

    gsl_matrix_set(&R.matrix, 2, 0, R31);
    gsl_matrix_set(&R.matrix, 2, 1, R32);
    gsl_matrix_set(&R.matrix, 2, 2, R33);

    euler_Rq(q, &Rq.matrix);

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

          calc_quaternions_ECEF(theta, phi, pos, vel, q);
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
