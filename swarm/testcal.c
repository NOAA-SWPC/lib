/*
 * Test Swarm scalar calibration
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>

#include <satdata/satdata.h>

#include "common.h"
#include "magcal.h"
#include "msynth.h"

#include "filter.c"

gsl_vector *
docal(const satdata_mag *data)
{
  int s = 0;
  const size_t p = MAGCAL_P + 1; /* add 1 for timing shift */
  magcal_workspace *magcal_p;
  gsl_vector *c;                 /* calibration coefficients */

  magcal_p = magcal_alloc(data->n);
  c = gsl_vector_alloc(p);

  /* initialize coefficient vector */
  magcal_initcond(c);

  /* compute calibration parameters */
  s = magcal_proc(c, data, magcal_p);

  magcal_free(magcal_p);

  if (s)
    {
      gsl_vector_free(c);
      c = NULL;
    }

  return c;
}

int
print_data(const gsl_vector *c, satdata_mag *data)
{
  size_t i;
  msynth_workspace *msynth_p = msynth_read("swarmC.cof");

  i = 1;
  printf("# Field %zu: time (years)\n", i++);
  printf("# Field %zu: latitude (degrees)\n", i++);
  printf("# Field %zu: longitude (degrees)\n", i++);
  printf("# Field %zu: altitude (km)\n", i++);
  printf("# Field %zu: original X residual (nT)\n", i++);
  printf("# Field %zu: original Y residual (nT)\n", i++);
  printf("# Field %zu: original Z residual (nT)\n", i++);
  printf("# Field %zu: new X residual (nT)\n", i++);
  printf("# Field %zu: new Y residual (nT)\n", i++);
  printf("# Field %zu: new Z residual (nT)\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      double tyr = satdata_epoch2year(data->t[i]);
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double alt = data->altitude[i];
      double B_int[4];
      double B_orig[3], B_cal[4];

      if (data->flags[i])
        continue;

      msynth_eval(tyr, alt + data->R,
                  theta, phi, B_int, msynth_p);

      B_orig[0] = SATDATA_VEC_X(data->B, i);
      B_orig[1] = SATDATA_VEC_Y(data->B, i);
      B_orig[2] = SATDATA_VEC_Z(data->B, i);

      magcal_apply_cal(c, B_orig, B_cal);

      printf("%f %f %f %f %.12e %.12e %.12e %.12e %.12e %.12e\n",
             tyr,
             data->latitude[i],
             data->longitude[i],
             data->altitude[i],
             B_orig[0] - B_int[0],
             B_orig[1] - B_int[1],
             B_orig[2] - B_int[2],
             B_cal[0] - B_int[0],
             B_cal[1] - B_int[1],
             B_cal[2] - B_int[2]);
    }

  msynth_free(msynth_p);

  return 0;
}

int
print_parameters(const gsl_vector *m, const satdata_mag *data)
{
  int s = 0;
  double t0 = satdata_epoch2year(data->t[0]),
         t1 = satdata_epoch2year(data->t[data->n - 1]),
         tavg = (t1 + t0) / 2.0;
  double SX = gsl_vector_get(m, MAGCAL_IDX_SX);
  double SY = gsl_vector_get(m, MAGCAL_IDX_SY);
  double SZ = gsl_vector_get(m, MAGCAL_IDX_SZ);
  double OX = gsl_vector_get(m, MAGCAL_IDX_OX);
  double OY = gsl_vector_get(m, MAGCAL_IDX_OY);
  double OZ = gsl_vector_get(m, MAGCAL_IDX_OZ);
  double AXY = gsl_vector_get(m, MAGCAL_IDX_AXY);
  double AXZ = gsl_vector_get(m, MAGCAL_IDX_AXZ);
  double AYZ = gsl_vector_get(m, MAGCAL_IDX_AYZ);
  FILE *fp = stdout;

  /* print parameter values */
  fprintf(fp, "%f %.12e %.12e %.12e %.12e %.12e %.12e %f %f %f\n",
         tavg,
         SX,
         SY,
         SZ,
         OX,
         OY,
         OZ,
         AXY * 180.0 / M_PI,
         AXZ * 180.0 / M_PI,
         AYZ * 180.0 / M_PI);
  fflush(fp);

  return s;
}

int
main(int argc, char *argv[])
{
  satdata_mag *data;
  struct timeval tv0, tv1;
  int c;
  char *infile = NULL;
  char *outfile = NULL;
  gsl_vector *coeffs;

  while ((c = getopt(argc, argv, "i:o:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'o':
            outfile = optarg;
            break;

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i swarm_idx_file> [-o output_cdf_file]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = satdata_swarm_read_idx(infile, 0);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
          time_diff(tv0, tv1));

  swarm_filter(5.0, 22.0, data);

  {
    size_t nflagged = satdata_nflagged(data);
    size_t nleft = data->n - nflagged;

    if (nleft < 2000)
      {
        fprintf(stderr, "main: insufficient data points for calibration\n");
        exit(1);
      }
  }

  coeffs = docal(data);
  if (coeffs == NULL)
    {
      fprintf(stderr, "main: error occurred in calibration\n");
      exit(1);
    }

  if (outfile)
    {
      magcal_apply(coeffs, data);

      fprintf(stderr, "main: writing calibrated data to %s...", outfile);
      satdata_swarm_write(outfile, data);
      fprintf(stderr, "done\n");

      print_parameters(coeffs, data);
    }
  /*else
    print_data(coeffs, data);*/

  return 0;
}
