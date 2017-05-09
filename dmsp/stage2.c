/*
 * stage2.c
 *
 * 1. Read DMSP cdf file
 * 2. Filter out day-time data and disturbed times
 * 3. Compute timing shift from filtered scalar data
 * 4. Compute 9 calibration parameters (scale factors, offsets,
 *    non-orthogonality angles) from vector data
 * 5. Iterate steps 3-4 until time shift and scalar calibration parameters
 *    converge
 * 6. Search calibrated scalar residuals for large outliers and flag
 * 7. Re-compute time shift and scalar calibration ignoring outliers from step 6.
 * 8. Apply time shift and calibration to dataset
 * 9. Compute Euler angles and apply to dataset
 *
 * Usage: ./stage2 <-i residual_index_file> <-o residual_output_file>
 *                 [-f] [-p parameter_file]
 *
 *
 * Files output:
 * If WRITE_DEBUG_FILES is defined, the following files are written:
 *
 * step0.dat - initial dataset before data selection (ASCII)
 * step1.dat - dataset after dF selection (ASCII)
 * step2.dat - dataset after QD latitude selection (ASCII)
 * step3.dat - dataset after local time selection selection (ASCII)
 * step4.dat - dataset after WMM data selection (ASCII)
 * step5.dat - final calibrated dataset with outliers removed (ASCII)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_spline.h>

#include <satdata/satdata.h>
#include <indices/indices.h>

#include "common.h"
#include "eph.h"
#include "magcal.h"
#include "msynth.h"
#include "track.h"

#include "stage2_filter.c"

/*
correct_track()
  Correct a track for possible (multiple) data jumps

Inputs: tol    - tolerance to detect a jump (nT)
        offset - (input/output)
                 on output,
                   offset(:,0) = offsets to add to X VFM data
                   offset(:,1) = offsets to add to Y VFM data
                   offset(:,2) = offsets to add to Z VFM data
                   offset(:,3) = offsets to add to F data
        data   - satellite data
        tptr   - track data
*/

int
correct_track(const double tol, gsl_matrix * offset, satdata_mag * data, track_data * tptr)
{
  int s = 0;
  double alpha[4] = { 0.0, 0.0, 0.0, 0.0 };
  size_t j;

  for (j = 0; j < tptr->n - 1; ++j)
    {
      size_t didx = tptr->start_idx + j;
      double dt = (data->t[didx + 1] - data->t[didx]) * 1.0e-3;
      double *offset_F = gsl_matrix_ptr(offset, j, 3);
      double *next_offset_F = gsl_matrix_ptr(offset, j + 1, 3);
      double delta[4];
      double B_VFM[4], B_VFM_next[4], B_model[4], B_model_next[4];

      *offset_F += alpha[3];

#if 0
      if (fabs(dt) > 1.0)
        continue; /* data gap */
#endif

      satdata_mag_model(didx, B_model, data);
      satdata_mag_model(didx + 1, B_model_next, data);

      B_VFM[3] = data->F[didx] + *offset_F - alpha[3];

      B_VFM_next[3] = data->F[didx + 1] + *next_offset_F;

      /* compute difference between adjacent samples and check for jump */
      delta[3] = (B_VFM_next[3] - B_model_next[3]) - (B_VFM[3] - B_model[3]);

      if (fabs(delta[3]) > tol)
        alpha[3] -= delta[3];
    }

  return s;
}

int
fix_track_jumps(const char *filename, satdata_mag *data, track_workspace *w)
{
  int s = 0;
  const double tol = 4.0; /* maximum allowed jump in nT for any vector component between adjacent samples */
  FILE *fp;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "fix_track_jumps: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: timestamp\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: field model F (nT)\n", i++);
  fprintf(fp, "# Field %zu: original F (nT)\n", i++);
  fprintf(fp, "# Field %zu: corrected F (nT)\n", i++);

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      gsl_matrix *offset = gsl_matrix_calloc(tptr->n, 4);
      size_t j;

      s = correct_track(tol, offset, data, tptr);
      if (s)
        return s;

      for (j = 0; j < tptr->n; ++j)
        {
          size_t didx = tptr->start_idx + j;
          time_t unix_time = satdata_epoch2timet(data->t[didx]);
          double offset_F = gsl_matrix_get(offset, j, 3);
          double B_model[4];

          satdata_mag_model(didx, B_model, data);

          fprintf(fp, "%ld %.4f %.4f %.4f %.4f\n",
                  unix_time,
                  data->qdlat[didx],
                  B_model[3],
                  data->F[didx],
                  data->F[didx] + offset_F);
        }

      fprintf(fp, "\n\n");

      gsl_matrix_free(offset);
    }

  fclose(fp);

  return s;
}

int
main(int argc, char *argv[])
{
  const char *track_file = "track_data.dat";
  char *outfile = NULL;
  satdata_mag *data = NULL;
  eph_data *eph = NULL;
  track_workspace *track_p;
  magcal_workspace *magcal_p;
  gsl_vector *coef = gsl_vector_alloc(MAGCAL_P);
  struct timeval tv0, tv1;
  double lt_min = 6.0;  /* local time interval for data selection at low/mid latitudes */
  double lt_max = 18.0;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "lt_min", required_argument, NULL, 'c' },
          { "lt_max", required_argument, NULL, 'd' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:o:b:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = satdata_dmsp_read_idx(optarg, 1);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
                    time_diff(tv0, tv1));
            break;

          case 'o':
            outfile = optarg;
            break;

          case 'b':
            fprintf(stderr, "main: reading Bowman ephemerides from %s...", optarg);
            gettimeofday(&tv0, NULL);
            eph = eph_data_read_bowman(optarg);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu read, %g seconds)\n", eph->n, time_diff(tv0, tv1));

          case 'c':
            lt_min = atof(optarg);
            break;

          case 'd':
            lt_max = atof(optarg);
            break;

          default:
            break;
        }
    }

  if (!data)
    {
      fprintf(stderr, "Usage: %s <-i dmsp_index_file> <-b bowman_ephemeris_file> [-o output_file] [--lt_min lt_min] [--lt_max lt_max]\n",
              argv[0]);
      exit(1);
    }

  if (outfile)
    fprintf(stderr, "output file = %s\n", outfile);

  track_p = track_alloc();

  fprintf(stderr, "main: initializing tracks...");
  gettimeofday(&tv0, NULL);
  track_init(data, NULL, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: filtering tracks for quiet periods...");
  gettimeofday(&tv0, NULL);
  stage2_filter(track_p, data);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* add unflagged data to magcal workspace */
  fprintf(stderr, "main: adding data for scalar calibration...");
  gettimeofday(&tv0, NULL);

  {
    size_t nflagged = satdata_nflagged(data);
    size_t n = data->n - nflagged;
    size_t i, j;

    magcal_p = magcal_alloc(n);

    for (i = 0; i < track_p->n; ++i)
      {
        track_data *tptr = &(track_p->tracks[i]);
        size_t start_idx = tptr->start_idx;
        size_t end_idx = tptr->end_idx;

        if (tptr->flags)
          continue;

        for (j = start_idx; j <= end_idx; ++j)
          {
            double B_VFM[3], B_model[4];

            if (data->flags[j])
              continue;

            if (fabs(data->qdlat[j]) > 55.0)
              continue;

            B_VFM[0] = SATDATA_VEC_X(data->B_VFM, j);
            B_VFM[1] = SATDATA_VEC_Y(data->B_VFM, j);
            B_VFM[2] = SATDATA_VEC_Z(data->B_VFM, j);

            satdata_mag_model(j, B_model, data);

            magcal_add_datum(data->t[j], B_VFM, B_model[3], magcal_p);
          }
      }
  }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* set initial values of calibration parameters */
  gsl_vector_set(coef, MAGCAL_IDX_SX, 1.0);
  gsl_vector_set(coef, MAGCAL_IDX_SY, 1.0);
  gsl_vector_set(coef, MAGCAL_IDX_SZ, 1.0);
  gsl_vector_set(coef, MAGCAL_IDX_OX, 0.0);
  gsl_vector_set(coef, MAGCAL_IDX_OY, 0.0);
  gsl_vector_set(coef, MAGCAL_IDX_OZ, 0.0);
  gsl_vector_set(coef, MAGCAL_IDX_AXY, M_PI / 2.0);
  gsl_vector_set(coef, MAGCAL_IDX_AXZ, M_PI / 2.0);
  gsl_vector_set(coef, MAGCAL_IDX_AYZ, M_PI / 2.0);

  fprintf(stderr, "main: performing scalar calibration...");
  gettimeofday(&tv0, NULL);
  magcal_proc(coef, magcal_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

#if 0
  fprintf(stderr, "main: fixing track jumps...");
  fix_track_jumps(track_file, data, track_p);
  fprintf(stderr, "done (data written to %s)\n", track_file);
#endif

  gsl_vector_free(coef);
  magcal_free(magcal_p);
  track_free(track_p);
  satdata_mag_free(data);
  eph_data_free(eph);

  return 0;
}
