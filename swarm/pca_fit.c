/*
 * pca_fit.c
 *
 * Fit PCA modes to satellite data track by track
 *
 * Pre-processing steps are:
 * 1. Instrument flags
 * 2. Track rms test
 * 3. Downsample by factor 15
 *
 * Usage: ./pca_fit -s swarmA.idx -r swarmC.idx -x swarmB.idx
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <satdata/satdata.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>

#include "bin3d.h"
#include "common.h"
#include "msynth.h"
#include "oct.h"
#include "pca.h"
#include "track.h"

/* define this to allow only 1 camera for data selection */
#define POLTOR_ONE_CAMERA          0

/* define to fit pca modes to C/B data */
#define FIT_PCA_C                  1
#define FIT_PCA_B                  1

/* relative weightings of different components */
#define PCA_WEIGHT_X               (5.0)
#define PCA_WEIGHT_Y               (1.0)
#define PCA_WEIGHT_Z               (5.0)

/* assign higher weight to low-latitude data for better EEJ fit */
#define PCA_WEIGHT_EEJ             (3.0)

typedef struct
{
  int all;            /* print all tracks */
  double lt_min;      /* local time minimum (hours) */
  double lt_max;      /* local time maximum (hours) */
  double alt_min;     /* altitude minimum (km) */
  double alt_max;     /* altitude maximum (km) */
  double qd_min;      /* minimum required QD latitude (deg) */
  double qd_max;      /* maximum required QD latitude (deg) */
  double lon_min;     /* minimum longitude (deg) */
  double lon_max;     /* maximum longitude (deg) */
  double kp_min;      /* minimum kp */
  double kp_max;      /* maximum kp */
  size_t downsample;  /* downsampling factor */
  double alpha;       /* smoothing factor for high latitudes */
  double thresh[4];   /* rms thresholds */
} preprocess_parameters;

satdata_mag *
read_swarm(const char *filename)
{
  size_t nflag;
  satdata_mag *data;
  struct timeval tv0, tv1;

  fprintf(stderr, "read_swarm: reading %s...", filename);
  gettimeofday(&tv0, NULL);
  data = satdata_swarm_read_idx(filename, 0);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu data read, %g seconds)\n",
          data->n, time_diff(tv0, tv1));

  /* check for instrument flags since we use Stage1 data */

  fprintf(stderr, "read_swarm: filtering for instrument flags...");
  nflag = satdata_swarm_filter_instrument(1, data);
  fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
          nflag, data->n, (double)nflag / (double)data->n * 100.0);

  return data;
} /* read_swarm() */

int
print_track_stats(const satdata_mag *data, const track_workspace *track_p)
{
  size_t nflagged = satdata_nflagged(data);
  size_t nleft = data->n - nflagged;
  size_t nflagged_track = track_nflagged(track_p);
  size_t nleft_track = track_p->n - nflagged_track;

  fprintf(stderr, "preprocess_data: total flagged data: %zu/%zu (%.1f%%)\n",
          nflagged, data->n, (double)nflagged / (double)data->n * 100.0);
  fprintf(stderr, "preprocess_data: total remaining data: %zu/%zu (%.1f%%)\n",
          nleft, data->n, (double)nleft / (double)data->n * 100.0);

  fprintf(stderr, "preprocess_data: total flagged tracks: %zu/%zu (%.1f%%)\n",
          nflagged_track, track_p->n, (double)nflagged_track / (double)track_p->n * 100.0);
  fprintf(stderr, "preprocess_data: total remaining tracks: %zu/%zu (%.1f%%)\n",
          nleft_track, track_p->n, (double)nleft_track / (double)track_p->n * 100.0);

  return 0;
} /* print_track_stats() */

/*
preprocess_data()

Inputs: params - preprocess parameters
          lt_min     - minimum local time
          lt_max     - maximum local time
          downsample - downsampling factor
        data   - satellite data

Return: pointer to sorted track workspace (should be freed by caller)
*/

track_workspace *
preprocess_data(const preprocess_parameters *params, satdata_mag *data)
{
  struct timeval tv0, tv1;
  track_workspace *track_p = track_alloc();

  fprintf(stderr, "preprocess_data: initializing tracks...");
  gettimeofday(&tv0, NULL);
  track_init(data, NULL, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  if (params->all)
    {
      print_track_stats(data, track_p);
      return track_p; /* no further filtering */
    }

  {
    const char *rmsfile = "satrms.dat";
    size_t nrms;

    nrms = track_flag_rms(rmsfile, params->thresh, NULL, data, track_p);
    fprintf(stderr, "preprocess_data: flagged (%zu/%zu) (%.1f%%) tracks due to high rms\n",
            nrms, track_p->n, (double) nrms / (double) track_p->n * 100.0);
  }

  /* flag local time */
  if (params->lt_min >= 0.0 && params->lt_max >= 0.0)
    {
      size_t nlt = track_flag_lt(params->lt_min, params->lt_max, NULL, data, track_p);

      fprintf(stderr, "preprocess_data: flagged data outside LT window [%g,%g]: %zu/%zu (%.1f%%) tracks flagged)\n",
              params->lt_min, params->lt_max,
              nlt, track_p->n, (double)nlt / (double)track_p->n * 100.0);
    }

  /* flag altitude */
  {
    size_t nalt = track_flag_meanalt(params->alt_min, params->alt_max, data, track_p);

    fprintf(stderr, "preprocess_data: flagged data due to altitude: %zu/%zu (%.1f%%) data flagged)\n",
            nalt, data->n, (double)nalt / (double)data->n * 100.0);
  }

  /* flag high kp data */
  {
    size_t nkp = track_flag_kp(params->kp_min, params->kp_max, data, track_p);

    fprintf(stderr, "preprocess_data: flagged data outside kp window [%g,%g]: %zu/%zu (%.1f%%) data flagged)\n",
            params->kp_min, params->kp_max,
            nkp, data->n, (double)nkp / (double)data->n * 100.0);
  }

  /* flag longitude */
  {
    size_t nlon = track_flag_lon(params->lon_min, params->lon_max, NULL, data, track_p);

    fprintf(stderr, "preprocess_data: flagged data due to longitude: %zu/%zu (%.1f%%) tracks flagged)\n",
            nlon, track_p->n, (double)nlon / (double)track_p->n * 100.0);
  }

  /* print track statistics */
  {
    char *jstat_file = "track_jump_stats.dat";

    fprintf(stderr, "preprocess_data: printing jump track statistics to %s...", jstat_file);
    track_print_stats_flag(jstat_file, TRACK_FLG_JUMP, track_p);
    fprintf(stderr, "done\n");
  }

  /* downsample data */
  {
    size_t i;

    fprintf(stderr, "preprocess_data: downsampling data by factor %zu...", params->downsample);

    for (i = 0; i < data->n; ++i)
      {
        if (i % params->downsample != 0)
          data->flags[i] |= SATDATA_FLG_DOWNSAMPLE;
      }

    fprintf(stderr, "done\n");
  }

#if 0
  if (params->alpha > 0.0)
    {
      fprintf(stderr, "preprocess_data: smoothing high-latitude residuals...");
      track_smooth(params->alpha, data, track_p);
      fprintf(stderr, "done\n");
    }

  fprintf(stderr, "preprocess_data: removing magnetospheric offsets...");
  gettimeofday(&tv0, NULL);
  track_fix_offsets(data, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
#endif

  print_track_stats(data, track_p);

  return track_p;
} /* preprocess_data() */

int
build_matrix_row(const double r, const double theta, const double phi,
                 gsl_vector *X, gsl_vector *Y, gsl_vector *Z,
                 pca_workspace *w)
{
  const size_t p = X->size;
  size_t i;

  for (i = 0; i < p; ++i)
    {
      double B[3];

      pca_pc_B(i, r, theta, phi, B, w);

      gsl_vector_set(X, i, B[0]);
      gsl_vector_set(Y, i, B[1]);
      gsl_vector_set(Z, i, B[2]);
    }

  return 0;
}

int
write_pca_track(const int header, FILE *fp, const satdata_mag *data,
                const gsl_vector *c, const size_t track_idx,
                track_workspace *track_p, pca_workspace *pca_p)
{
  track_data *tptr;
  size_t j;

  if (header)
    {
      j = 1;
      fprintf(fp, "# Field %zu: QD latitude (degrees)\n", j++);
      fprintf(fp, "# Field %zu: B_x observed (nT)\n", j++);
      fprintf(fp, "# Field %zu: B_y observed (nT)\n", j++);
      fprintf(fp, "# Field %zu: B_z observed (nT)\n", j++);
      fprintf(fp, "# Field %zu: B_x modeled (nT)\n", j++);
      fprintf(fp, "# Field %zu: B_y modeled (nT)\n", j++);
      fprintf(fp, "# Field %zu: B_z modeled (nT)\n", j++);

      return 0;
    }

  tptr = &(track_p->tracks[track_idx]);

  for (j = 0; j < tptr->n; ++j)
    {
      size_t didx = j + tptr->start_idx;
      double r = data->altitude[didx] + data->R;
      double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
      double phi = data->longitude[didx] * M_PI / 180.0;
      double B_res[3], B_model[3];

      if (SATDATA_BadData(data->flags[didx]))
        continue;

      /* compute magnetic residual */
      track_residual(j, didx, B_res, data, track_p);

      /* compute PCA model at this point */
      pca_B(c, r, theta, phi, B_model, pca_p);

      fprintf(fp, "%f %f %f %f %f %f %f\n",
              data->qdlat[didx],
              B_res[0],
              B_res[1],
              B_res[2],
              B_model[0],
              B_model[1],
              B_model[2]);
    }

  fprintf(fp, "\n\n");

  return 0;
}

/* add a track Bx,By,Bz into LS matrix and rhs vector */
int
add_LS(const track_data *tptr, const satdata_mag *data, track_workspace *track_p,
       gsl_matrix *X, gsl_vector *rhs, gsl_vector *wts, size_t *ridx, pca_workspace *pca_p)
{
  size_t rowidx = *ridx;
  size_t j;

  for (j = 0; j < tptr->n; ++j)
    {
      size_t didx = j + tptr->start_idx;
      double r = data->altitude[didx] + data->R;
      double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
      double phi = data->longitude[didx] * M_PI / 180.0;
      double B_res[3];
      gsl_vector_view vx, vy, vz;
      double wj = 1.0;

      if (!SATDATA_AvailableData(data->flags[didx]))
        continue;

      /* fit only low-latitude data */
      if (fabs(data->qdlat[didx]) > 40.0)
        continue;

      /* upweight equatorial data to get a better EEJ fit */
      if (fabs(data->qdlat[didx]) <= 10.0)
        wj = PCA_WEIGHT_EEJ;

      /* compute magnetic residual */
      track_residual(j, didx, B_res, data, track_p);

      /* set rhs vector */
      gsl_vector_set(rhs, rowidx, B_res[0]);
      gsl_vector_set(rhs, rowidx + 1, B_res[1]);
      gsl_vector_set(rhs, rowidx + 2, B_res[2]);

      /* set weight vector */
      gsl_vector_set(wts, rowidx, PCA_WEIGHT_X * wj);
      gsl_vector_set(wts, rowidx + 1, PCA_WEIGHT_Y * wj);
      gsl_vector_set(wts, rowidx + 2, PCA_WEIGHT_Z * wj);

      /* build 3 rows of the LS matrix */
      vx = gsl_matrix_row(X, rowidx);
      vy = gsl_matrix_row(X, rowidx + 1);
      vz = gsl_matrix_row(X, rowidx + 2);

      build_matrix_row(r, theta, phi, &vx.vector, &vy.vector, &vz.vector, pca_p);

      rowidx += 3;
    }

  *ridx = rowidx;

  return 0;
}

int
pcafit(const size_t track_idx, const satdata_mag *data, track_workspace *track_p,
       const size_t track2_idx, const satdata_mag *data2, track_workspace *track2_p,
       const size_t track3_idx, const satdata_mag *data3, track_workspace *track3_p,
       gsl_vector *c, pca_workspace *pca_p)
{
  track_data *tptr = &(track_p->tracks[track_idx]);
  track_data *tptr2 = NULL;
  track_data *tptr3 = NULL;
  const size_t p = c->size;                  /* number of PCs / model parameters */
  size_t n = 3 * tptr->n;                    /* number of data */
  size_t rowidx = 0;

  gsl_multifit_linear_workspace *multifit_p;
  gsl_matrix *X;
  gsl_vector *rhs, *wts;

  if (track2_p)
    {
      tptr2 = &(track2_p->tracks[track2_idx]);
      n += 3 * tptr2->n;
    }

  if (track3_p)
    {
      tptr3 = &(track3_p->tracks[track3_idx]);
      n += 3 * tptr3->n;
    }

  multifit_p = gsl_multifit_linear_alloc(n, p);
  X = gsl_matrix_alloc(n, p);
  rhs = gsl_vector_alloc(n);
  wts = gsl_vector_alloc(n);

  /* add Swarm A data to LS matrix */
  add_LS(tptr, data, track_p, X, rhs, wts, &rowidx, pca_p);

  /* add Swarm C data to LS matrix */
  if (tptr2)
    add_LS(tptr2, data2, track2_p, X, rhs, wts, &rowidx, pca_p);

  /* add Swarm B data to LS matrix */
  if (tptr3)
    add_LS(tptr3, data3, track3_p, X, rhs, wts, &rowidx, pca_p);

  /* solve regularized LS system */
  {
    double lambda;
    gsl_matrix_view A = gsl_matrix_submatrix(X, 0, 0, rowidx, p);
    gsl_vector_view b = gsl_vector_subvector(rhs, 0, rowidx);
    gsl_vector_view w = gsl_vector_subvector(wts, 0, rowidx);
    double rcond, rnorm, snorm;

    const size_t nL = 200;
    gsl_vector *reg_param = gsl_vector_alloc(nL);
    gsl_vector *rho = gsl_vector_alloc(nL);
    gsl_vector *eta = gsl_vector_alloc(nL);
    size_t reg_idx;

    /* apply weight matrix */
    gsl_multifit_linear_applyW(&A.matrix, &w.vector, &b.vector,
                               &A.matrix, &b.vector);

#if 0
    print_octave(&A.matrix, "A2");
    printv_octave(&b.vector, "b2");
    exit(1);
#endif

    /* compute SVD */
    gsl_multifit_linear_svd(&A.matrix, multifit_p);

    rcond = gsl_multifit_linear_rcond(multifit_p);
    fprintf(stderr, "cond(A) = %g...", 1.0 / rcond);

    gsl_multifit_linear_lcurve(&b.vector, reg_param, rho, eta, multifit_p);
    gsl_multifit_linear_lcorner(rho, eta, &reg_idx);
    lambda = gsl_vector_get(reg_param, reg_idx);
#if 0
    lambda = GSL_MAX(5.0, lambda);
#endif

    /* solve LS system */
    gsl_multifit_linear_solve(lambda, &A.matrix, &b.vector, c, &rnorm, &snorm, multifit_p);

    fprintf(stderr, "lambda = %g...", lambda);

    gsl_vector_free(reg_param);
    gsl_vector_free(rho);
    gsl_vector_free(eta);
  }

  gsl_multifit_linear_free(multifit_p);
  gsl_matrix_free(X);
  gsl_vector_free(rhs);
  gsl_vector_free(wts);

  return 0;
}

int
main_proc(const satdata_mag *data, const satdata_mag *data2, const satdata_mag *data3,
          track_workspace *track1, track_workspace *track2, track_workspace *track3)
{
  int s = 0;
#if 0 /* XXX */
  const size_t p = 43;   /* number of PCs / model parameters */
#else
  const size_t p = 14;   /* number of PCs / model parameters */
#endif
  size_t i, j, k;
  size_t nflagged, nunflagged;
  pca_workspace *pca_p = pca_alloc();
  gsl_vector *c = gsl_vector_alloc(p);
  const char *file1 = "data1.txt";
  const char *file2 = "data2.txt";
  const char *file3 = "data3.txt";
  const char *file_chi = "chi.txt";
  FILE *fp1 = fopen(file1, "w");
  FILE *fp2 = fopen(file2, "w");
  FILE *fp3 = fopen(file3, "w");
  FILE *fp_chi = fopen(file_chi, "w");
  size_t idx = 0;
  char buf[2048];

  nflagged = track_nflagged(track1);
  nunflagged = track1->n - nflagged;
  fprintf(stderr, "Total tracks:    %zu\n", track1->n);
  fprintf(stderr, "Total flagged:   %zu\n", nflagged);
  fprintf(stderr, "Total unflagged: %zu\n", nunflagged);

  /* write headers */
  write_pca_track(1, fp1, NULL, NULL, 0, NULL, pca_p);
  write_pca_track(1, fp2, NULL, NULL, 0, NULL, pca_p);
  write_pca_track(1, fp3, NULL, NULL, 0, NULL, pca_p);

  for (i = 0; i < track1->n; ++i)
    {
      track_data *tptr = &(track1->tracks[i]);
      track_data *tptr2, *tptr3;
      double dphi;
      time_t unix_time;

      if (tptr->flags != 0)
        continue;

      /* find Swarm C crossing within 1 min and 1.7 deg of A */
      s = track_find(tptr->t_eq, tptr->lon_eq, 1.0, 1.7, &j, track2);
      if (s)
        continue;

      tptr2 = &(track2->tracks[j]);
      if (tptr2->flags != 0)
        continue;

      /* find Swarm B crossing */
      s = track_find(tptr->t_eq, tptr->lon_eq, 5.0, 20.0, &k, track3);
      if (s)
        continue;

      tptr3 = &(track3->tracks[k]);
      if (tptr3->flags != 0)
        continue;

      dphi = wrap180(tptr->lon_eq - tptr3->lon_eq);
      if (fabs(dphi) < 10.0)
        continue;

      unix_time = satdata_epoch2timet(tptr->t_eq);
      sprintf(buf, "%s", ctime(&unix_time));
      buf[strlen(buf) - 1] = '\0';

      fprintf(stderr, "main_proc: found track triple, %s, dt = %f [min], dphi = %f [deg]\n",
              buf,
              (tptr->t_eq - tptr3->t_eq) / 60000.0,
              dphi);
      fprintf(stderr, "\t LT A = %f\n", tptr->lt_eq);
      fprintf(stderr, "\t LT C = %f\n", tptr2->lt_eq);
      fprintf(stderr, "\t LT B = %f\n", tptr3->lt_eq);
      fprintf(stderr, "\t longitude A = %f [deg]\n", tptr->lon_eq);
      fprintf(stderr, "\t longitude C = %f [deg]\n", tptr2->lon_eq);
      fprintf(stderr, "\t longitude B = %f [deg]\n", tptr3->lon_eq);

      fprintf(stderr, "main_proc: fitting PCs to track %zu/%zu (index %zu)...", i + 1, track1->n, idx);

      pcafit(i, data, track1,
#if FIT_PCA_C
             j, data2, track2,
#else
             0, NULL, NULL,
#endif
#if FIT_PCA_B
             k, data3, track3,
#else
             0, NULL, NULL,
#endif
             c, pca_p);

      fprintf(stderr, "done\n");

      /* write PCA along-track prediction for all satellites */
      write_pca_track(0, fp1, data, c, i, track1, pca_p);
      write_pca_track(0, fp2, data2, c, j, track2, pca_p);
      write_pca_track(0, fp3, data3, c, k, track3, pca_p);

      /* write PCA current map */
      fprintf(stderr, "main_proc: writing current map for index %zu to %s...", idx, file_chi);
      pca_print_map(fp_chi, R_EARTH_KM + 450.0, c, pca_p);
      fprintf(stderr, "done\n");

      ++idx;
    }

  pca_free(pca_p);
  gsl_vector_free(c);

  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp_chi);

  return s;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --swarm_file | -s swarm_index_file          - Swarm index file\n");
  fprintf(stderr, "\t --swarm_file2 | -r swarm_index_file2        - Swarm index file 2\n");
  fprintf(stderr, "\t --swarm_file3 | -x swarm_index_file3        - Swarm index file 3\n");
  fprintf(stderr, "\t --champ_file | -c champ_index_file          - CHAMP index file\n");
  fprintf(stderr, "\t --champ_plp_file | -b champ_plp_index_file  - CHAMP PLP index file\n");
  fprintf(stderr, "\t --all | -a                                  - print all tracks (no filtering)\n");
  fprintf(stderr, "\t --downsample | -d downsample                - downsampling factor\n");
  fprintf(stderr, "\t --output_file | -o output_file              - output file\n");
  fprintf(stderr, "\t --lt_min | -j lt_min                        - local time minimum\n");
  fprintf(stderr, "\t --lt_max | -k lt_max                        - local time maximum\n");
  fprintf(stderr, "\t --alt_min | -l alt_min                      - altitude minimum\n");
  fprintf(stderr, "\t --alt_max | -m alt_max                      - altitude maximum\n");
  fprintf(stderr, "\t --lon_min | -t lon_min                      - longitude minimum\n");
  fprintf(stderr, "\t --lon_max | -u lon_max                      - longitude maximum\n");
  fprintf(stderr, "\t --kp_min | -v kp_min                        - kp minimum\n");
  fprintf(stderr, "\t --kp_max | -w kp_max                        - kp maximum\n");
  fprintf(stderr, "\t --alpha | -q alpha                          - smoothing factor for high latitudes\n");
}

int
main(int argc, char *argv[])
{
  char *data_file = "track_data.dat";
  char *stats_file = "track_stats.dat";
  char *output_file = NULL;
  satdata_mag *data = NULL;
  satdata_mag *data2 = NULL;
  satdata_mag *data3 = NULL;
  satdata_lp *lp_data = NULL;
  struct timeval tv0, tv1;
  track_workspace *track_p, *track2_p, *track3_p;
  preprocess_parameters params;

  /* defaults */
  params.all = 0;
  params.lt_min = -1.0;
  params.lt_max = -1.0;
  params.alt_min = 0.0;
  params.alt_max = 1000.0;
  params.qd_min = -30.0;
  params.qd_max = 30.0;
  params.lon_min = -200.0;
  params.lon_max = 200.0;
  params.kp_min = 0.0;
  params.kp_max = 2.0;
  params.downsample = 1;
  params.alpha = -1.0;
  params.thresh[0] = 210.0;
  params.thresh[1] = 170.0;
  params.thresh[2] = 150.0;
  params.thresh[3] = 160.0;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "all", no_argument, NULL, 'a' },
          { "champ_file", required_argument, NULL, 'c' },
          { "swarm_file", required_argument, NULL, 's' },
          { "champ_plp_file", required_argument, NULL, 'b' },
          { "downsample", required_argument, NULL, 'd' },
          { "lt_min", required_argument, NULL, 'j' },
          { "lt_max", required_argument, NULL, 'k' },
          { "alt_min", required_argument, NULL, 'l' },
          { "alt_max", required_argument, NULL, 'm' },
          { "lon_min", required_argument, NULL, 't' },
          { "lon_max", required_argument, NULL, 'u' },
          { "kp_min", required_argument, NULL, 'v' },
          { "kp_max", required_argument, NULL, 'w' },
          { "output_file", required_argument, NULL, 'o' },
          { "alpha", required_argument, NULL, 'q' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "ab:c:d:j:k:l:m:o:q:r:s:t:u:v:x:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            params.all = 1;
            break;

          case 's':
            data = read_swarm(optarg);
            break;

          case 'r':
            data2 = read_swarm(optarg);
            break;

          case 'x':
            data3 = read_swarm(optarg);
            break;

          case 'c':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = satdata_champ_read_idx(optarg, 0);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu data read, %g seconds)\n",
                    data->n, time_diff(tv0, tv1));

            /* check for instrument flags since we use Stage1 data */
            {
              size_t nflag;
              size_t champ_flags = 0;

#if POLTOR_ONE_CAMERA
              /* allow only 1 camera in data selection */
              champ_flags = SATDATA_FLG_ONESC;
#endif

              fprintf(stderr, "main: filtering for instrument flags...");
              nflag = satdata_champ_filter_instrument(1, champ_flags, data);
              fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
                      nflag, data->n, (double)nflag / (double)data->n * 100.0);
            }

            break;

          case 'b':
            fprintf(stderr, "main: reading CHAMP PLP data from %s...", optarg);
            lp_data = satdata_champ_plp_read_idx(optarg, 0);
            fprintf(stderr, "done\n");
            break;

          case 'd':
            params.downsample = (size_t) atoi(optarg);
            break;

          case 'o':
            output_file = optarg;
            break;

          case 'j':
            params.lt_min = atof(optarg);
            break;

          case 'k':
            params.lt_max = atof(optarg);
            break;

          case 'l':
            params.alt_min = atof(optarg);
            break;

          case 'm':
            params.alt_max = atof(optarg);
            break;

          case 't':
            params.lon_min = atof(optarg);
            break;

          case 'u':
            params.lon_max = atof(optarg);
            break;

          case 'v':
            params.kp_min = atof(optarg);
            break;

          case 'w':
            params.kp_max = atof(optarg);
            break;

          case 'q':
            params.alpha = atof(optarg);
            break;

          default:
            break;
        }
    }

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  fprintf(stderr, "main: LT minimum       = %.1f\n", params.lt_min);
  fprintf(stderr, "main: LT maximum       = %.1f\n", params.lt_max);
  fprintf(stderr, "main: altitude minimum = %.1f\n", params.alt_min);
  fprintf(stderr, "main: altitude maximum = %.1f\n", params.alt_max);
  fprintf(stderr, "main: QD minimum       = %.1f\n", params.qd_min);
  fprintf(stderr, "main: QD maximum       = %.1f\n", params.qd_max);
  fprintf(stderr, "main: lon minimum      = %.1f\n", params.lon_min);
  fprintf(stderr, "main: lon maximum      = %.1f\n", params.lon_max);
  fprintf(stderr, "main: kp minimum       = %.1f\n", params.kp_min);
  fprintf(stderr, "main: kp maximum       = %.1f\n", params.kp_max);
  fprintf(stderr, "main: smoothing alpha  = %f\n", params.alpha);

  if (lp_data)
    {
      fprintf(stderr, "main: adding electron densities to magnetic data...");
      satdata_mag_fill_ne(data, lp_data);
      fprintf(stderr, "done\n");
    }

  track_p = preprocess_data(&params, data);

  if (data2)
    track2_p = preprocess_data(&params, data2);

  if (data3)
    track3_p = preprocess_data(&params, data3);

  main_proc(data, data2, data3, track_p, track2_p, track3_p);

  satdata_mag_free(data);
  track_free(track_p);

  if (data2)
    {
      satdata_mag_free(data2);
      track_free(track2_p);
    }

  return 0;
} /* main() */
