/*
 * compare.c
 * Patrick Alken
 *
 * Compare Swarm ASMV with VFM data
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <assert.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include <satdata/satdata.h>
#include <indices/indices.h>

#include "bsearch.h"
#include "common.h"
#include "euler.h"

/* define to compare vectors in NEC frame, instead of VFM */
#define COMPARE_NEC      1

#define MAX_DATA         100000

typedef struct
{
  double t[MAX_DATA];  /* timestamp of ascending node (CDF_EPOCH) */
  double lt[MAX_DATA]; /* LT of ascending node in hours */
  size_t n;
} lt_data;

int
lt_init(const satdata_mag *data, lt_data *ltdata)
{
  int s = 0;
  size_t i;
  size_t n = 0;

  for (i = 0; i < data->n - 1; ++i)
    {
      /* search for ascending equator crossing */
      if (data->latitude[i] <= 0.0 &&
          data->latitude[i + 1] > 0.0)
        {
          double phi = data->longitude[i] * M_PI / 180.0;
          double lt = satdata_localtime(data->t[i], phi);

          ltdata->t[n] = data->t[i];
          ltdata->lt[n] = lt;
          ++n;
        }
    }

  ltdata->n = n;

  return s;
} /* lt_init() */

/*
lt_angles()
  Determine sun angles for a given time

Inputs: sat_idx - index into data
        data    - satellite data
        alpha   - (output) alpha angle in radians in [0,2pi]
        beta    - (output) beta angle in radians in [0,2pi]
        ltdata  - lt data

Return: success/error
*/

int
lt_angles(const size_t sat_idx, const satdata_mag *data, double *alpha, double *beta,
          const lt_data *ltdata)
{
  int s = 0;
  size_t idx;
  int satdir = satdata_mag_satdir(sat_idx, data);
  double t = data->t[sat_idx];
  double theta = M_PI / 2.0 - data->latitude[sat_idx] * M_PI / 180.0;

  if (t < ltdata->t[0] || t > ltdata->t[ltdata->n - 1])
    {
      /*fprintf(stderr, "lt_angles: error: t out of range\n");*/
      return -1;
    }

  idx = bsearch_double(ltdata->t, t, 0, ltdata->n - 1);

  /* LT of ascending node */
  *beta = 2.0 * M_PI * (ltdata->lt[idx] / 24.0);

  /* angle from equator of ascending node */
  if (satdir == 1)
    *alpha = wrap2pi(M_PI / 2.0 - theta);
  else
    *alpha = wrap2pi(M_PI / 2.0 + theta);

  return s;
} /* lt_angles() */

int
compare_asmv(const char *filename, const satdata_mag *data, const satdata_mag *data_asmv,
             const lt_data *ltdata)
{
  int s = 0;
  size_t i, j, k;
  const double tol = 1.0e-3;
  double rms[3] = { 0.0, 0.0, 0.0 };
  size_t n = 0;
  FILE *fp;

  if (data->n == 0 || data_asmv->n == 0)
    {
      fprintf(stderr, "compare_asmv: error: empty dataset\n");
      return -1;
    }

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "compare_asmv: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: time (UT)\n", i++);
  fprintf(fp, "# Field %zu: time (decimal year)\n", i++);
  fprintf(fp, "# Field %zu: UT (hours)\n", i++);
  fprintf(fp, "# Field %zu: local time (hours)\n", i++);
  fprintf(fp, "# Field %zu: season (day of year)\n", i++);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: satellite direction\n", i++);
  fprintf(fp, "# Field %zu: alpha (degrees)\n", i++);
  fprintf(fp, "# Field %zu: beta (degrees)\n", i++);
  fprintf(fp, "# Field %zu: scalar field (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual X field, VFM frame (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual Y field, VFM frame (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual Z field, VFM frame (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual F field (nT)\n", i++);
  fprintf(fp, "# Field %zu: scalar residual |B_asm| - F (nT)\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      double B_vfm[4], B_asmv[4], B_diff[3];
      double alpha, beta;

      j = bsearch_double(data_asmv->t, data->t[i], 0, data_asmv->n - 1);
      if (data->t[i] != data_asmv->t[j])
        continue;

      if (data->flags[i] != 0 || data_asmv->flags[j] != 0)
        continue;

      /* make sure both datasets match exactly */
      assert(data->t[i] == data_asmv->t[j]);
      gsl_test_rel(data_asmv->latitude[j], data->latitude[i], tol, "i=%zu latitude", i);
      gsl_test_rel(data_asmv->altitude[j], data->altitude[i], tol, "i=%zu altitude", i);
      gsl_test_rel(data_asmv->longitude[j], data->longitude[i], tol, "i=%zu longitude", i);
      gsl_test_rel(data_asmv->F[j], data->F[i], tol, "i=%zu F", i);

      s = lt_angles(i, data, &alpha, &beta, ltdata);
      if (s)
        continue;

#if COMPARE_NEC

      /* VFM vector in VFM frame */
      B_vfm[0] = SATDATA_VEC_X(data->B, i);
      B_vfm[1] = SATDATA_VEC_Y(data->B, i);
      B_vfm[2] = SATDATA_VEC_Z(data->B, i);
      B_vfm[3] = vec_norm(B_vfm);

#else

      /* VFM vector in VFM frame */
      B_vfm[0] = SATDATA_VEC_X(data->B_VFM, i);
      B_vfm[1] = SATDATA_VEC_Y(data->B_VFM, i);
      B_vfm[2] = SATDATA_VEC_Z(data->B_VFM, i);
      B_vfm[3] = vec_norm(B_vfm);

#endif

      /*
       * The ASMV 'B' array will either be VFM or NEC, depending on whether
       * user supplied Euler angles
       */

      /* ASMV vector in VFM frame */
      B_asmv[0] = SATDATA_VEC_X(data_asmv->B, j);
      B_asmv[1] = SATDATA_VEC_Y(data_asmv->B, j);
      B_asmv[2] = SATDATA_VEC_Z(data_asmv->B, j);
      B_asmv[3] = vec_norm(B_asmv);

      if (B_vfm[3] < 50.0 || B_asmv[3] < 50.0)
        continue;

      /* normalize so that |B_asm| = |B_vfm| */
      for (k = 0; k < 3; ++k)
        B_asmv[k] *= B_vfm[3] / B_asmv[3];
      B_asmv[3] = vec_norm(B_asmv);

      /* compute residuals in VFM frame */
      for (k = 0; k < 4; ++k)
        B_diff[k] = B_vfm[k] - B_asmv[k];

      if (fabs(B_diff[0]) > 50.0 ||
          fabs(B_diff[1]) > 50.0 ||
          fabs(B_diff[2]) > 50.0)
        continue;

      for (k = 0; k < 3; ++k)
        rms[k] += pow(B_diff[k], 2.0);

      ++n;

      if (i % 20 == 0)
        {
          time_t unix_time = satdata_epoch2timet(data->t[i]);
          double phi = data->longitude[i] * M_PI / 180.0;
          double lt = get_localtime(unix_time, phi);

          fprintf(fp, "%ld %.8f %.2f %.2f %6.2f %10.4f %8.4f %10.4f %10.4f %d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
                  unix_time,
                  satdata_epoch2year(data->t[i]),
                  get_ut(unix_time),
                  lt,
                  get_season(unix_time),
                  data->longitude[i],
                  data->latitude[i],
                  data->altitude[i],
                  data->qdlat[i],
                  satdata_mag_satdir(i, data),
                  alpha * 180.0 / M_PI,
                  beta * 180.0 / M_PI,
                  data->F[i],
                  B_diff[0],
                  B_diff[1],
                  B_diff[2],
                  B_diff[3],
                  B_asmv[3] - data->F[i]);
        }
    }

  if (n > 0)
    {
      for (k = 0; k < 3; ++k)
        rms[k] = sqrt(rms[k] / (double) n);
    }

  fprintf(stderr, "number of data = %zu\n", n);
  fprintf(stderr, "X rms = %.2f [nT]\n", rms[0]);
  fprintf(stderr, "Y rms = %.2f [nT]\n", rms[1]);
  fprintf(stderr, "Z rms = %.2f [nT]\n", rms[2]);

  fclose(fp);

  return s;
}

int
compare_data(const char *filename, const satdata_mag *data1, const satdata_mag *data2,
             const lt_data *ltdata)
{
  int s = 0;
  size_t i, j, k;
  const double tol = 1.0e-3;
  double rms[3] = { 0.0, 0.0, 0.0 };
  size_t n = 0;
  FILE *fp;

  if (data1->n == 0 || data2->n == 0)
    {
      fprintf(stderr, "compare_data: error: empty dataset\n");
      return -1;
    }

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "compare_data: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: time (UT)\n", i++);
  fprintf(fp, "# Field %zu: time (decimal year)\n", i++);
  fprintf(fp, "# Field %zu: UT (hours)\n", i++);
  fprintf(fp, "# Field %zu: local time (hours)\n", i++);
  fprintf(fp, "# Field %zu: season (day of year)\n", i++);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: satellite direction\n", i++);
  fprintf(fp, "# Field %zu: alpha (degrees)\n", i++);
  fprintf(fp, "# Field %zu: beta (degrees)\n", i++);
  fprintf(fp, "# Field %zu: scalar field (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual X field, VFM frame (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual Y field, VFM frame (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual Z field, VFM frame (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual F field (nT)\n", i++);
  fprintf(fp, "# Field %zu: scalar residual |B_asm| - F (nT)\n", i++);

  for (i = 0; i < data1->n; ++i)
    {
      double B_vfm1[4], B_vfm2[4], B_diff[3];
      double alpha, beta;

      j = bsearch_double(data2->t, data1->t[i], 0, data2->n - 1);
      if (data1->t[i] != data2->t[j])
        continue;

      if (data1->flags[i] != 0 || data2->flags[j] != 0)
        continue;

      /* make sure both datasets match exactly */
      assert(data1->t[i] == data2->t[j]);
      gsl_test_rel(data2->latitude[j], data1->latitude[i], tol, "i=%zu latitude", i);
      gsl_test_rel(data2->altitude[j], data1->altitude[i], tol, "i=%zu altitude", i);
      gsl_test_rel(data2->longitude[j], data1->longitude[i], tol, "i=%zu longitude", i);
      gsl_test_rel(data2->F[j], data1->F[i], tol, "i=%zu F", i);

      s = lt_angles(i, data1, &alpha, &beta, ltdata);
      if (s)
        continue;

#if COMPARE_NEC

      /* VFM vector in NEC frame */
      B_vfm1[0] = SATDATA_VEC_X(data1->B, i);
      B_vfm1[1] = SATDATA_VEC_Y(data1->B, i);
      B_vfm1[2] = SATDATA_VEC_Z(data1->B, i);

      B_vfm2[0] = SATDATA_VEC_X(data2->B, j);
      B_vfm2[1] = SATDATA_VEC_Y(data2->B, j);
      B_vfm2[2] = SATDATA_VEC_Z(data2->B, j);

#else

      /* VFM vector in VFM frame */
      B_vfm1[0] = SATDATA_VEC_X(data1->B_VFM, i);
      B_vfm1[1] = SATDATA_VEC_Y(data1->B_VFM, i);
      B_vfm1[2] = SATDATA_VEC_Z(data1->B_VFM, i);

      /* VFM vector in VFM frame */
      B_vfm2[0] = SATDATA_VEC_X(data2->B_VFM, j);
      B_vfm2[1] = SATDATA_VEC_Y(data2->B_VFM, j);
      B_vfm2[2] = SATDATA_VEC_Z(data2->B_VFM, j);

#endif

      B_vfm1[3] = vec_norm(B_vfm1);
      B_vfm2[3] = vec_norm(B_vfm2);

      if (B_vfm1[3] < 50.0 || B_vfm2[3] < 50.0)
        continue;

      /* compute residuals in VFM frame */
      for (k = 0; k < 4; ++k)
        B_diff[k] = B_vfm1[k] - B_vfm2[k];

      if (fabs(B_diff[0]) > 50.0 ||
          fabs(B_diff[1]) > 50.0 ||
          fabs(B_diff[2]) > 50.0)
        continue;

      for (k = 0; k < 3; ++k)
        rms[k] += pow(B_diff[k], 2.0);

      ++n;

      if (i % 20 == 0)
        {
          time_t unix_time = satdata_epoch2timet(data1->t[i]);
          double phi = data1->longitude[i] * M_PI / 180.0;
          double lt = get_localtime(unix_time, phi);

          fprintf(fp, "%ld %.8f %.2f %.2f %6.2f %10.4f %8.4f %10.4f %10.4f %d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
                  unix_time,
                  satdata_epoch2year(data1->t[i]),
                  get_ut(unix_time),
                  lt,
                  get_season(unix_time),
                  data1->longitude[i],
                  data1->latitude[i],
                  data1->altitude[i],
                  data1->qdlat[i],
                  satdata_mag_satdir(i, data1),
                  alpha * 180.0 / M_PI,
                  beta * 180.0 / M_PI,
                  data1->F[i],
                  B_diff[0],
                  B_diff[1],
                  B_diff[2],
                  B_diff[3],
                  B_vfm2[3] - data1->F[i]);
        }
    }

  for (k = 0; k < 3; ++k)
    rms[k] = sqrt(rms[k] / (double) n);

  fprintf(stderr, "number of data = %zu\n", n);
  fprintf(stderr, "X rms = %.2f [nT]\n", rms[0]);
  fprintf(stderr, "Y rms = %.2f [nT]\n", rms[1]);
  fprintf(stderr, "Z rms = %.2f [nT]\n", rms[2]);

  fclose(fp);

  return s;
}

int
main(int argc, char *argv[])
{
  satdata_mag *data1 = NULL;
  satdata_mag *data2 = NULL;
  satdata_mag *data = NULL;
  satdata_mag *data_asmv = NULL;
  euler_workspace *euler_asmv_p = NULL;
  euler_workspace *euler_p1 = NULL;
  euler_workspace *euler_p2 = NULL;
  char *outfile = "data.dat";
  struct timeval tv0, tv1;
  lt_data ltdata;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "asmv_file", required_argument, NULL, 'a' },
          { "input_file", required_argument, NULL, 'i' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:e:f:i:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            if (!data1)
              {
                data1 = satdata_swarm_read_idx(optarg, 0);
                data = data1;
              }
            else if (!data2)
              {
                data2 = satdata_swarm_read_idx(optarg, 0);
                data = data2;
              }
            if (!data)
              exit(1);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
                    time_diff(tv0, tv1));
            break;

          case 'a':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data_asmv = satdata_swarm_read_idx(optarg, 0);
            if (!data_asmv)
              exit(1);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data_asmv->n,
                    time_diff(tv0, tv1));
            break;

          case 'e':
            fprintf(stderr, "main: reading Euler angles from %s...", optarg);
            euler_asmv_p = euler_read(optarg);
            if (!euler_asmv_p)
              exit(1);
            fprintf(stderr, "done (%zu sets of angles read)\n", euler_asmv_p->n);
            break;

          case 'f':
            fprintf(stderr, "main: reading Euler angles from %s...", optarg);
            if (!euler_p1)
              {
                euler_p1 = euler_read(optarg);
                if (!euler_p1)
                  exit(1);
                fprintf(stderr, "done (%zu sets of angles read)\n", euler_p1->n);
              }
            else
              {
                euler_p2 = euler_read(optarg);
                if (!euler_p2)
                  exit(1);
                fprintf(stderr, "done (%zu sets of angles read)\n", euler_p2->n);
              }
            break;

          default:
            break;
        }
    }

  if (!data1)
    {
      fprintf(stderr, "Usage: %s <-a swarm_asmv_index_file> <-i swarm_index_file> <-e euler_angle_asmv_file> <-f swarm_euler_angle_file>\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: initialzing LT structure...");
  lt_init(data1, &ltdata);
  fprintf(stderr, "done\n");

  if (euler_asmv_p && data_asmv)
    {
      fprintf(stderr, "main: rotating ASMV measurements with new Euler angles...");
      euler_apply(data_asmv, euler_asmv_p);
      fprintf(stderr, "done\n");
    }

  if (euler_p1 && data1)
    {
      fprintf(stderr, "main: rotating Swarm measurements with new Euler angles...");
      euler_apply(data1, euler_p1);
      fprintf(stderr, "done\n");
    }

  if (euler_p2 && data2)
    {
      fprintf(stderr, "main: rotating Swarm measurements with new Euler angles...");
      euler_apply(data2, euler_p2);
      fprintf(stderr, "done\n");
    }

  if (data_asmv)
    {
      fprintf(stderr, "main: writing data to %s...", outfile);
      compare_asmv(outfile, data1, data_asmv, &ltdata);
      fprintf(stderr, "done\n");
    }
  else if (data2)
    {
      fprintf(stderr, "main: writing data to %s...", outfile);
      compare_data(outfile, data1, data2, &ltdata);
      fprintf(stderr, "done\n");
    }

#if 0
  fprintf(stderr, "main: writing data to %s...", outfile);
  print_data(outfile, data, data_asmv);
  fprintf(stderr, "done\n");
#endif

  if (data1)
    satdata_mag_free(data1);
  if (data2)
    satdata_mag_free(data2);
  if (data_asmv)
    satdata_mag_free(data_asmv);

  return 0;
} /* main() */
