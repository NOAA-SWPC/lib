/*
 * necorrJ.c
 *
 * Correlate gravity line currents with electron density peaks
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
#include <flow/flow.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>

#include "apex.h"
#include "bin3d.h"
#include "bsearch.h"
#include "common.h"
#include "interp.h"
#include "peak.h"
#include "track.h"

/* maximum number of peaks in J data */
#define MAX_PEAKS                  10

#define MIN_KP                     (0.0)
#define MAX_KP                     (2.0)

/* number of line currents in each profile */
#define NCURR                      240

/* maximum number of profiles per satellite */
#define MAX_PROFILES               10000

#define MAX_BUFFER                 4096

typedef struct
{
  time_t t;        /* timestamp of equator crossing */
  double phi;      /* longitude of equator crossing (deg) */
  double J[NCURR]; /* line current values in mA/m */

  double qdlat_peak;
  double J_peak;

  size_t peak_idx[MAX_PEAKS]; /* index into J of peak locations */
  size_t npeak;               /* number of peaks found */
} curr_profile;

typedef struct
{
  curr_profile profiles[MAX_PROFILES];
  size_t n;
} current_data;

int
read_lc(const char *filename, current_data *data)
{
  FILE *fp;
  char buffer[MAX_BUFFER];

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "read_lc: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  data->n = 0;

  while (fgets(buffer, MAX_BUFFER, fp) != NULL)
    {
      char *str = buffer;
      int offset, c;
      size_t i;
      curr_profile *prof = &data->profiles[data->n];

      if (*str == '#')
        continue;

      c = sscanf(str, "%ld %lf %n",
                 &(prof->t), &(prof->phi), &offset);
      if (c < 2)
        continue;

      str += offset;
      for (i = 0; i < NCURR; ++i)
        {
          c = sscanf(str, "%lf %n", &(prof->J[i]), &offset);
          if (c < 1)
            {
              fprintf(stderr, "read_lc: error reading %s\n", filename);
              fclose(fp);
              return -1;
            }
          str += offset;
        }

      ++(data->n);
    }

  fclose(fp);

  return 0;
}

/* interpolate Ne data near time t to latitude qdlat */
double
interp_ne(const time_t t, const double qdlat, const satdata_lp *data)
{
  double t_cdf = satdata_timet2epoch(t);
  double Ne;

  if (t_cdf >= data->t[0] && t_cdf <= data->t[data->n - 1])
    {
      int idx = (int) bsearch_double(data->t, t_cdf, 0, data->n - 1);
      int half_window = 15; /* number of minutes in half time window */
      int idx1 = GSL_MAX(0, idx - half_window*60*2);
      int idx2 = GSL_MIN(data->n - 1, idx + half_window*60*2);
      int qidx;
      
      if (data->qdlat[idx1] < data->qdlat[idx2])
        qidx = (int) bsearch_double(data->qdlat, qdlat, idx1, idx2);
      else
        qidx = (int) bsearch_desc_double(data->qdlat, qdlat, idx1, idx2);

      Ne = interp1d(data->qdlat[qidx], data->qdlat[qidx + 1],
                    data->ne[qidx], data->ne[qidx + 1], qdlat);
    }
  else
    Ne = 0.0;

  return Ne;
}

/*
find_J_peaks()
  Go through each profile in 'data' and search for
1 peak in the north hemisphere, and 1 peak in the south hemisphere.

Inputs: peak_qd_min - search for peaks in [peak_qd_min,peak_qd_max]
        peak_qd_max - see above
        minheight   - minimum height needed for peak (mA/m)
        peak_file   - output file for peaks
        prof_file   - output file for profiles
        lp_data     - Langmuir probe data
        data        - profile data

On output:
1) If exactly 1 peak is found, 'npeak' is set to 1 and
the qdlat_peak and J_peak fields are set appropriately
*/

int
find_J_peaks(const double peak_qd_min, const double peak_qd_max,
             const double minheight, const char *peak_file,
             const char *prof_file, const satdata_lp *lp_data,
             current_data *data)
{
  int s;
  FILE *fp_peak, *fp_prof;
  const size_t width_points = 20;                                 /* number of points in half-width of peak */
  const size_t smooth_window = width_points / 2;                  /* smoothing window size */
  const double minslope = 0.5 * pow((double) width_points, -2.0); /* minimum slope of first derivative */
  peak_workspace *peak_p = peak_alloc(NCURR);

  const double qdlat_max = 30.0;
  const double qdlat_step = 2.0 * qdlat_max / (NCURR - 1.0);
  double *qdlat = malloc(NCURR * sizeof(double));
  gsl_vector_view qdvec = gsl_vector_view_array(qdlat, NCURR);
  size_t i, j;
  size_t cnt = 0;

  fp_peak = fopen(peak_file, "w");
  if (!fp_peak)
    {
      fprintf(stderr, "find_J_peaks: unable to open %s: %s\n",
              peak_file, strerror(errno));
      return -1;
    }

  fp_prof = fopen(prof_file, "w");
  if (!fp_prof)
    {
      fprintf(stderr, "find_J_peaks: unable to open %s: %s\n",
              prof_file, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp_peak, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp_peak, "# Field %zu: current strength (mA/m)\n", i++);
  fprintf(fp_peak, "# Field %zu: derivative of current (mA/m/deg)\n", i++);
  fprintf(fp_peak, "# Field %zu: smoothed derivative of current (mA/m/deg)\n", i++);
  fprintf(fp_peak, "# Field %zu: electron density (cm^{-3})\n", i++);
  fprintf(fp_peak, "# Field %zu: current peak density (mA/m)\n", i++);

  i = 1;
  fprintf(fp_prof, "# Field %zu: timestamp of equator crossing (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
  fprintf(fp_prof, "# Field %zu: longitude of equator crossing (degrees)\n", i++);
  fprintf(fp_prof, "# Field %zu: QD latitude of J peak (degrees)\n", i++);
  fprintf(fp_prof, "# Field %zu: current density of J peak (mA/m)\n", i++);

  for (i = 0; i < NCURR; ++i)
    qdlat[i] = -qdlat_max + i * qdlat_step;

  /* loop over current profiles */
  for (i = 0; i < data->n; ++i)
    {
      curr_profile *prof = &(data->profiles[i]);
      gsl_vector_view J = gsl_vector_view_array(prof->J, NCURR);

      prof->npeak = 0;

      /* initialize peak search */
      peak_init(smooth_window, &qdvec.vector, &J.vector, peak_p);

      do
        {
          s = peak_find(1, minslope, minheight, &qdvec.vector, &J.vector, peak_p);

          if (s == GSL_CONTINUE)
            {
              size_t pidx = peak_p->pidx; /* index into arrays of peak location */
              double center = qdlat[pidx];

              /* test if peak center is in [peak_qd_min,peak_qd_max] */
              if (center >= peak_qd_min && center <= peak_qd_max)
                {
                  /* found a peak, store location */
                  prof->peak_idx[prof->npeak] = pidx;
                  ++(prof->npeak);

                  if (prof->npeak >= MAX_PEAKS)
                    break;
                }
            }
        }
      while (s == GSL_CONTINUE);

      if (prof->npeak == 1)
        {
          size_t idx = prof->peak_idx[0];

          prof->qdlat_peak = qdlat[idx];
          prof->J_peak = prof->J[idx];

          fprintf(fp_prof, "%ld %8.3f %8.3f %6.2f\n",
                  prof->t,
                  prof->phi,
                  prof->qdlat_peak,
                  prof->J_peak);

          ++cnt;
        }
      else
        {
          /* more/less than 1 peak found, ignore */
          prof->npeak = 0;
        }

      for (j = 0; j < NCURR; ++j)
        {
          double J = 1.0 / 0.0;
          double Ne = interp_ne(prof->t, qdlat[j], lp_data);

          if (prof->npeak == 1)
            {
              if (fabs(qdlat[j] - prof->qdlat_peak) < 10.0)
                J = prof->J[j];
            }

          fprintf(fp_peak, "%f %f %f %f %e %f\n",
                  qdlat[j],
                  prof->J[j],
                  peak_deriv(j, peak_p),
                  peak_sderiv(j, peak_p),
                  Ne,
                  J);
        }

      fprintf(fp_peak, "\n\n");

    }

  fprintf(stderr, "find_J_peaks: found peak in %zu/%zu profiles\n",
          cnt, data->n);

  peak_free(peak_p);
  free(qdlat);

  fclose(fp_peak);
  fclose(fp_prof);

  return 0;
}

size_t
find_ne_peaks(const double qd_min, const double qd_max, const size_t width_points,
              const gsl_vector *qdlat, const gsl_vector *ne, size_t idx[],
              peak_workspace *peak_p)
{
  int s;
  const size_t smooth_window = width_points / 2;                  /* smoothing window size */
  const double minslope = 0.4 * pow((double) width_points, -2.0); /* minimum slope of first derivative */
  const double minheight = 8.0e5;                                 /* minimum peak height */
  size_t npeak = 0;

  /* initialize peak search */
  peak_init(smooth_window, qdlat, ne, peak_p);

  do
    {
      s = peak_find(1, minslope, minheight, qdlat, ne, peak_p);

      if (s == GSL_CONTINUE)
        {
          size_t pidx = peak_p->pidx; /* index into arrays of peak location */
          double center = gsl_vector_get(qdlat, pidx);

          /* test if peak center is in [qd_min,qd_max] */
          if (fabs(center) >= qd_min && fabs(center) <= qd_max)
            {
              /* found a peak, store Gaussian parameters */
              idx[npeak++] = pidx;

              if (npeak >= MAX_PEAKS)
                break;
            }
        }
    }
  while (s == GSL_CONTINUE);

  return npeak;
}

size_t
find_F_peaks(const double qd_min, const double qd_max, const size_t width_points,
             const gsl_vector *qdlat, const gsl_vector *F, size_t idx[],
             peak_workspace *peak_p)
{
  int s;
  const size_t smooth_window = width_points / 2;                  /* smoothing window size */
  const double minslope = 1.0 * pow((double) width_points, -2.0); /* minimum slope of first derivative */
  const double minheight = -1.0e6;                                /* minimum peak height */
  size_t npeak = 0;

  /* initialize peak search */
  peak_init(smooth_window, qdlat, F, peak_p);

  do
    {
      /*s = peak_find(-1, minslope, minheight, qdlat, F, peak_p);*/
      s = peak_find(0, minslope, minheight, qdlat, F, peak_p);

      if (s == GSL_CONTINUE)
        {
          size_t pidx = peak_p->pidx; /* index into arrays of peak location */
          double center = gsl_vector_get(qdlat, pidx);

          /* test if peak center is in [qd_min,qd_max] */
          if (center >= qd_min && center <= qd_max)
            {
              /* found a peak, store Gaussian parameters */
              idx[npeak++] = pidx;

              if (npeak >= MAX_PEAKS)
                break;
            }
        }
    }
  while (s == GSL_CONTINUE);

  return npeak;
}

int
correlate(const char *peak_file, const char *corr_file, const satdata_mag *data, const track_workspace *w)
{
  int s = 0;
  size_t i, j;
  const size_t width_points = 100;                                /* number of points in half-width of peak */
  const size_t fit_width = width_points * 1.5;                    /* gaussian fit window size */
  const size_t half_width = fit_width / 2;
  FILE *fp_peak, *fp_corr;

  size_t ne_idx[MAX_PEAKS];                                       /* indices for centers of Ne peaks */
  size_t F_idx[MAX_PEAKS];                                        /* indices for centers of F peaks */
  gsl_vector *ne_coef[MAX_PEAKS];                                 /* save gauss parameters for each peak */
  const size_t p = 3;

  const double qd_min = 4.0;                                      /* search for peaks only in [qd_min,qd_max] */
  const double qd_max = 30.0;

  fp_peak = fopen(peak_file, "w");
  fp_corr = fopen(corr_file, "w");

  i = 1;
  fprintf(fp_peak, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp_peak, "# Field %zu: electron density (cm^{-3})\n", i++);
  fprintf(fp_peak, "# Field %zu: scalar field residual (nT)\n", i++);
  fprintf(fp_peak, "# Field %zu: N_e first derivative (cm^{-3})\n", i++);
  fprintf(fp_peak, "# Field %zu: smoothed N_e first derivative (cm^{-3})\n", i++);
  fprintf(fp_peak, "# Field %zu: F first derivative (nT)\n", i++);
  fprintf(fp_peak, "# Field %zu: smoothed F first derivative (nT)\n", i++);
  fprintf(fp_peak, "# Field %zu: N_e data used for 1st peak fit (cm^{-3})\n", i++);
  fprintf(fp_peak, "# Field %zu: N_e data used for 2nd peak fit (cm^{-3})\n", i++);
  fprintf(fp_peak, "# Field %zu: scalar field peak 1 used for correlation (nT)\n", i++);
  fprintf(fp_peak, "# Field %zu: scalar field peak 2 used for correlation (nT)\n", i++);

  i = 1;
  fprintf(fp_corr, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
  fprintf(fp_corr, "# Field %zu: time (decimal year)\n", i++);
  fprintf(fp_corr, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp_corr, "# Field %zu: local time (hours)\n", i++);
  fprintf(fp_corr, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp_corr, "# Field %zu: season (day of year)\n", i++);
  fprintf(fp_corr, "# Field %zu: N_e peak center in northern hemisphere (degrees QD latitude)\n", i++);
  fprintf(fp_corr, "# Field %zu: N_e peak height in northern hemisphere (cm^{-3})\n", i++);
  fprintf(fp_corr, "# Field %zu: N_e peak center in southern hemisphere (degrees QD latitude)\n", i++);
  fprintf(fp_corr, "# Field %zu: N_e peak height in southern hemisphere (cm^{-3})\n", i++);
  fprintf(fp_corr, "# Field %zu: F peak center in northern hemisphere (degrees QD latitude)\n", i++);
  fprintf(fp_corr, "# Field %zu: F peak height in northern hemisphere (nT)\n", i++);
  fprintf(fp_corr, "# Field %zu: F peak center in southern hemisphere (degrees QD latitude)\n", i++);
  fprintf(fp_corr, "# Field %zu: F peak height in southern hemisphere (nT)\n", i++);
  fprintf(fp_corr, "# Field %zu: corr(F,N_e) in northern hemisphere\n", i++);
  fprintf(fp_corr, "# Field %zu: corr(F,N_e) in southern hemisphere\n", i++);

  for (i = 0; i < MAX_PEAKS; ++i)
    ne_coef[i] = gsl_vector_alloc(p);

  for (i = 0; i < w->n; ++i)
    {
      time_t unix_time;
      track_data *tptr = &(w->tracks[i]);
      size_t sidx = tptr->start_idx;
      size_t n = tptr->n;
      peak_workspace *peak_ne = peak_alloc(n);
      peak_workspace *peak_F = peak_alloc(n);
      gsl_vector_view qd = gsl_vector_view_array(&(data->qdlat[sidx]), n);
      gsl_vector_view ne = gsl_vector_view_array(&(data->ne[sidx]), n);
      gsl_vector_view F = gsl_vector_view_array(tptr->Bf, n);
      size_t npeak_ne, npeak_F;
      double corr1, corr2;
      size_t ne_idx0[2], ne_idx1[2];
      size_t F_idx0[2], F_idx1[2];
      double ne_center1, ne_center2;

      /* discard flagged tracks */
      if (tptr->flags)
        continue;

      npeak_ne = find_ne_peaks(qd_min, qd_max, width_points, &qd.vector, &ne.vector, ne_idx, peak_ne);

      if (npeak_ne != 2)
        continue;

      /* store QD latitudes of peak centers */
      ne_center1 = gsl_vector_get(&qd.vector, ne_idx[0]);
      ne_center2 = gsl_vector_get(&qd.vector, ne_idx[1]);

      if (ne_center1 * ne_center2 > 0.0)
        {
          fprintf(stderr, "correlate: detected two peaks in the same hemisphere\n");
          continue;
        }

      /* search for corresponding peaks in F data */
      {
        const double qd_width = 5.0;

        /* search for an F peak within ne_center +/- qd_width */
        npeak_F = find_F_peaks(ne_center1 - qd_width, ne_center1 + qd_width, width_points, &qd.vector, &F.vector, F_idx, peak_F);
        if (npeak_F != 1)
          continue; /* no corresponding F peak found */

        npeak_F = find_F_peaks(ne_center2 - qd_width, ne_center2 + qd_width, width_points, &qd.vector, &F.vector, &(F_idx[1]), peak_F);
        if (npeak_F != 1)
          continue; /* no corresponding F peak found */
      }

      ne_idx0[0] = ne_idx[0] - half_width;
      ne_idx1[0] = ne_idx[0] + half_width;
      ne_idx0[1] = ne_idx[1] - half_width;
      ne_idx1[1] = ne_idx[1] + half_width;

      F_idx0[0] = F_idx[0] - half_width;
      F_idx1[0] = F_idx[0] + half_width;
      F_idx0[1] = F_idx[1] - half_width;
      F_idx1[1] = F_idx[1] + half_width;

      fprintf(stderr, "correlate: track %zu: found %zu N_e peaks\n", i + 1, npeak_ne);

      {
        size_t n1 = ne_idx1[0] - ne_idx0[0] + 1; /* number of points for correlation of peak 1 */
        size_t n2 = ne_idx1[1] - ne_idx0[1] + 1; /* number of points for correlation of peak 2 */

        corr1 = gsl_stats_correlation(&ne.vector.data[ne_idx0[0]], 1, &tptr->Bf[F_idx0[0]], 1, n1);
        corr2 = gsl_stats_correlation(&ne.vector.data[ne_idx0[1]], 1, &tptr->Bf[F_idx0[1]], 1, n2);

        /* make sure corr1 is always the peak in the northern hemisphere */
        if (ne_center1 < 0.0)
          {
            SWAP(corr1, corr2);
          }
      }

      for (j = 0; j < n; ++j)
        {
          double qdj = gsl_vector_get(&qd.vector, j);
          double ne1, ne2, F1, F2;

          if (j >= ne_idx0[0] && j <= ne_idx1[0])
            ne1 = gsl_vector_get(&ne.vector, j);
          else
            ne1 = 1.0/0.0;

          if (j >= ne_idx0[1] && j <= ne_idx1[1])
            ne2 = gsl_vector_get(&ne.vector, j);
          else
            ne2 = 1.0/0.0;

          if (j >= F_idx0[0] && j <= F_idx1[0])
            F1 = tptr->Bf[j];
          else
            F1 = 1.0/0.0;

          if (j >= F_idx0[1] && j <= F_idx1[1])
            F2 = tptr->Bf[j];
          else
            F2 = 1.0/0.0;

          fprintf(fp_peak, "%f %.12e %f %.12e %.12e %f %f %.12e %.12e %f %f\n",
                  qdj,
                  gsl_vector_get(&ne.vector, j),
                  tptr->Bf[j],
                  peak_deriv(j, peak_ne),
                  peak_sderiv(j, peak_ne),
                  peak_deriv(j, peak_F),
                  peak_sderiv(j, peak_F),
                  ne1,
                  ne2,
                  F1,
                  F2);
        }

      fprintf(fp_peak, "\n\n");

      unix_time = satdata_epoch2timet(tptr->t_eq);
      fprintf(fp_corr, "%ld %f %f %f %f %f %f %e %f %e %f %f %f %f %f %f\n",
              unix_time,
              satdata_epoch2year(tptr->t_eq),
              tptr->lon_eq,
              tptr->lt_eq,
              tptr->meanalt,
              get_season(unix_time),
              gsl_vector_get(&qd.vector, ne_idx[0]),
              gsl_vector_get(&ne.vector, ne_idx[0]),
              gsl_vector_get(&qd.vector, ne_idx[1]),
              gsl_vector_get(&ne.vector, ne_idx[1]),
              gsl_vector_get(&qd.vector, F_idx[0]),
              gsl_vector_get(&F.vector, F_idx[0]),
              gsl_vector_get(&qd.vector, F_idx[1]),
              gsl_vector_get(&F.vector, F_idx[1]),
              corr1,
              corr2);

      /*exit(1);*/

      peak_free(peak_ne);
      peak_free(peak_F);
    }

  for (i = 0; i < MAX_PEAKS; ++i)
    gsl_vector_free(ne_coef[i]);

  fclose(fp_peak);
  fclose(fp_corr);

  return s;
}

int
fill_qd(satdata_lp *lp_data)
{
  size_t i;
  time_t t = satdata_epoch2timet(lp_data->t[0]);
  int year = get_year(t);
  apex_workspace *apex_p = apex_alloc(year);

  for (i = 0; i < lp_data->n; ++i)
    {
      double alon, alat, qdlat;
      double phi = lp_data->longitude[i] * M_PI / 180.0;
      double theta = M_PI / 2.0 - lp_data->latitude[i] * M_PI / 180.0;
      double r = lp_data->r[i] * 1.0e3;

      apex_transform(theta, phi, r, &alon, &alat, &qdlat,
                     NULL, NULL, NULL, apex_p);
      lp_data->qdlat[i] = qdlat;
    }

  apex_free(apex_p);

  return 0;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --curr_file | -j line_current_file          - Line current file\n");
  fprintf(stderr, "\t --champ_plp_file | -b champ_plp_index_file  - CHAMP PLP index file\n");
  fprintf(stderr, "\t --swarm_lp_file | -f swarm_lp_index_file    - Swarm LP index file\n");
  fprintf(stderr, "\t --all | -a                                  - print all tracks (no filtering)\n");
}

int
main(int argc, char *argv[])
{
  const char *peak_file1 = "peak1.dat";
  const char *prof_file1 = "prof1.dat";
  const char *corr_file = "corr.dat";
  satdata_lp *lp_data = NULL;
  current_data sat1;
  peak_workspace *peak_workspace_p;

  sat1.n = 0;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "curr_file", required_argument, NULL, 'j' },
          { "champ_plp_file", required_argument, NULL, 'b' },
          { "swarm_plp_file", required_argument, NULL, 'f' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "b:j:f:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'j':
            fprintf(stderr, "main: reading %s...", optarg);
            read_lc(optarg, &sat1);
            fprintf(stderr, "done (%zu profiles read)\n", sat1.n);
            break;

          case 'b':
            fprintf(stderr, "main: reading CHAMP PLP data from %s...", optarg);
            lp_data = satdata_champ_plp_read_idx(optarg, 0);
            fprintf(stderr, "done\n");
            break;

          case 'f':
            fprintf(stderr, "main: reading Swarm LP data from %s...", optarg);
            lp_data = satdata_swarm_lp_read_idx(optarg, 0);
            fprintf(stderr, "done\n");
            break;

          default:
            print_help(argv);
            exit(1);
            break;
        }
    }

  if (sat1.n == 0)
    {
      print_help(argv);
      exit(1);
    }

  if (lp_data)
    {
      fprintf(stderr, "main: filling QD latitudes for LP data...");
      fill_qd(lp_data);
      fprintf(stderr, "done\n");
    }

  peak_workspace_p = peak_alloc(NCURR);

  fprintf(stderr, "main: searching for northern J peaks in satellite 1...");
  find_J_peaks(5.0, 23.0, 1.0, peak_file1, prof_file1, lp_data, &sat1);
  fprintf(stderr, "done\n");
  
  fprintf(stderr, "main: peak data written to %s\n", peak_file1);
  fprintf(stderr, "main: profile data written to %s\n", prof_file1);

#if 0
  fprintf(stderr, "main: computing correlations...");
  correlate(peak_file, corr_file, data, track_p);
  fprintf(stderr, "done (outputs: %s, %s)\n", peak_file, corr_file);
#endif

  if (lp_data)
    satdata_lp_free(lp_data);

  peak_free(peak_workspace_p);

  return 0;
} /* main() */
