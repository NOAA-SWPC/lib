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
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

#include "apex.h"
#include <common/bin.h>
#include <common/bin2d.h>
#include <common/bsearch.h>
#include <common/common.h>
#include <common/interp.h>

#include "curvefit.h"
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

size_t find_ne_peaks(const double qd_min, const double qd_max, const size_t width_points,
                     const gsl_vector *qdlat, const gsl_vector *ne, size_t idx[],
                     peak_workspace *peak_p);

typedef struct
{
  time_t t;        /* timestamp of equator crossing */
  double phi;      /* longitude of equator crossing (deg) */
  double J[NCURR]; /* line current values in mA/m */

  double qdlat_J;
  double peak_J;

  double qdlat_Ne; /* QD latitude of Ne peak */
  double peak_Ne;  /* magnitude of Ne peak cm^{-3} */
  double ctr_Ne;   /* crest-to-trough of Ne peak cm^{-3} */

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
interp_ne(const time_t t, const double qdlat, const satdata_mag *data)
{
  double t_cdf = satdata_timet2epoch(t);
  double Ne;

  if (t_cdf >= data->t[0] && t_cdf <= data->t[data->n - 1])
    {
      int idx = (int) bsearch_double(data->t, t_cdf, 0, data->n - 1);
      int half_window = 15; /* number of minutes in half time window */
      int idx1 = GSL_MAX(0, idx - half_window*60*2);
      int idx2 = GSL_MIN(data->n - 1, idx + half_window*60*2);
      int qidx = (int) bsearch_double(data->qdlat, qdlat, idx1, idx2);

      Ne = interp1d(data->qdlat[qidx], data->qdlat[qidx + 1],
                    data->ne[qidx], data->ne[qidx + 1], qdlat);
    }
  else
    Ne = 0.0;

  return Ne;
}

/*
search_Ne_peak()
  Search for Ne density peak for track corresponding to time
t within [peak_qd_min,peak_qd_max] QD latitude
*/

int
search_Ne_peak(const time_t t, const double peak_qd_min,
               const double peak_qd_max, const track_workspace *track_p,
               const satdata_mag *data, curr_profile *prof)
{
  int s = 0;
  size_t i;
  int track_idx = -1;

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      time_t t_eq = satdata_epoch2timet(tptr->t_eq);

      if (abs(t - t_eq) < 10)
        {
          track_idx = (int) i;
          break;
        }
    }

  if (track_idx < 0)
    {
      fprintf(stderr, "search_Ne_peak: track not found for time %ld\n", t);
      return -1;
    }

  {
    const size_t width_points = 100;        /* number of points in half-width of peak */
    track_data *tptr = &(track_p->tracks[track_idx]);
    size_t sidx = tptr->start_idx;
    size_t n = tptr->n;
    gsl_vector_view qd = gsl_vector_view_array(&(data->qdlat[sidx]), n);
    gsl_vector_view ne = gsl_vector_view_array(&(data->ne[sidx]), n);
    size_t npeak;
    peak_workspace *peak_p = peak_alloc(n);
    size_t ne_idx[MAX_PEAKS];

    npeak = find_ne_peaks(peak_qd_min, peak_qd_max, width_points,
                          &qd.vector, &ne.vector, ne_idx, peak_p);

    /* if more than 1 peak found, choose the largest */
    if (npeak > 1)
      {
        double max_Ne = -1.0e6;
        size_t max_idx = 0;

        for (i = 0; i < npeak; ++i)
          {
            size_t idx = ne_idx[i];
            double Ne = gsl_vector_get(&ne.vector, idx);

            if (Ne > max_Ne)
              {
                max_Ne = Ne;
                max_idx = idx;
              }
          }

        ne_idx[0] = max_idx;
        npeak = 1;
      }

    if (npeak == 1)
      {
        gsl_vector *qdvec = gsl_vector_alloc(n);
        gsl_vector *nevec = gsl_vector_alloc(n);
        double eq_Ne; /* value of Ne at equator */

        gsl_vector_memcpy(qdvec, &qd.vector);
        gsl_vector_memcpy(nevec, &ne.vector);
        gsl_sort_vector2(qdvec, nevec);

        eq_Ne = interp_xy(qdvec->data, nevec->data, n, 0.0);

        prof->qdlat_Ne = gsl_vector_get(&qd.vector, ne_idx[0]);
        prof->peak_Ne = gsl_vector_get(&ne.vector, ne_idx[0]);
        prof->ctr_Ne = prof->peak_Ne - eq_Ne;

        gsl_vector_free(qdvec);
        gsl_vector_free(nevec);

        s = 0;
      }
    else
      s = -1; /* peak not found */

    peak_free(peak_p);
  }

  return s;
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
        mag_data    - magnetic data
        track_p     - track workspace
        data        - profile data

On output:
1) If exactly 1 peak is found, 'npeak' is set to 1 and
the qdlat_J and peak_J fields are set appropriately
*/

int
find_J_peaks(const double peak_qd_min, const double peak_qd_max,
             const double minheight, const char *peak_file,
             const char *prof_file, const satdata_mag *mag_data,
             const track_workspace *track_p, current_data *data)
{
  int s;
  FILE *fp_peak, *fp_prof;
  const size_t width_points = 15;                                 /* number of points in half-width of peak */
  const size_t smooth_window = width_points / 2;                  /* smoothing window size */
  const double minslope = 0.4 * pow((double) width_points, -2.0); /* minimum slope of first derivative */
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
  fprintf(fp_peak, "# Field %zu: Ne peak density (cm^{-3})\n", i++);

  i = 1;
  fprintf(fp_prof, "# Field %zu: timestamp of equator crossing (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
  fprintf(fp_prof, "# Field %zu: longitude of equator crossing (degrees)\n", i++);
  fprintf(fp_prof, "# Field %zu: QD latitude of J peak (degrees)\n", i++);
  fprintf(fp_prof, "# Field %zu: current density of J peak (mA/m)\n", i++);
  fprintf(fp_prof, "# Field %zu: QD latitude of Ne peak (degrees)\n", i++);
  fprintf(fp_prof, "# Field %zu: magnitude of Ne peak (cm^{-3})\n", i++);
  fprintf(fp_prof, "# Field %zu: crest-trough of Ne peak (cm^{-3})\n", i++);

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
              double J = prof->J[pidx];

              /* test if peak center is in [peak_qd_min,peak_qd_max] */
              if (fabs(J) <= 30.0 && center >= peak_qd_min && center <= peak_qd_max)
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

      /* if more than 1 peak found, choose the highest peak */
      if (prof->npeak > 1)
        {
          double max_J = -1.0e6;
          size_t max_idx = 0;

          for (j = 0; j < prof->npeak; ++j)
            {
              size_t idx = prof->peak_idx[j];

              if (prof->J[idx] > max_J)
                {
                  max_J = prof->J[idx];
                  max_idx = idx;
                }
            }

          prof->peak_idx[0] = max_idx;
          prof->npeak = 1;
        }

      if (prof->npeak == 1)
        {
          size_t idx = prof->peak_idx[0];

          prof->qdlat_J = qdlat[idx];
          prof->peak_J = prof->J[idx];
          prof->qdlat_Ne = 1.0 / 0.0;
          prof->peak_Ne = 1.0 / 0.0;
          prof->ctr_Ne = 1.0 / 0.0;

          if (mag_data)
            {
              /*s = search_Ne_peak(prof->t, peak_qd_min, peak_qd_max, track_p, mag_data, prof);*/
              s = search_Ne_peak(prof->t, prof->qdlat_J - 5.0, prof->qdlat_J + 5.0, track_p, mag_data, prof);
            }

          fprintf(fp_prof, "%ld %8.3f %8.3f %6.2f %8.3f %e %e\n",
                  prof->t,
                  prof->phi,
                  prof->qdlat_J,
                  prof->peak_J,
                  prof->qdlat_Ne,
                  prof->peak_Ne,
                  prof->ctr_Ne);

          ++cnt;
        }
      else
        {
          /* more/less than 1 peak found, ignore */
          prof->npeak = 0;
        }

      for (j = 0; j < NCURR; ++j)
        {
          double J_peak = 1.0 / 0.0;
          double Ne_peak = 1.0 / 0.0;
          double Ne = mag_data ? interp_ne(prof->t, qdlat[j], mag_data) : 0.0;

          if (prof->npeak == 1)
            {
              if (fabs(qdlat[j] - prof->qdlat_J) < 10.0)
                J_peak = prof->J[j];

              if (isfinite(prof->qdlat_Ne) &&
                  fabs(qdlat[j] - prof->qdlat_Ne) < 10.0)
                Ne_peak = Ne;
            }

          fprintf(fp_peak, "%f %f %f %f %e %f %e\n",
                  qdlat[j],
                  prof->J[j],
                  peak_deriv(j, peak_p),
                  peak_sderiv(j, peak_p),
                  Ne,
                  J_peak,
                  Ne_peak);
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

/* correlate centers of J/Ne peaks and magnitudes of J/Ne peaks */
int
docorr(current_data *data, double *corr_center, double *corr_peak)
{
  double qd_J[10000], qd_Ne[10000];
  double peak_J[10000], peak_Ne[10000];
  size_t n = 0;
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      curr_profile *prof = &(data->profiles[i]);

      if (prof->npeak == 0)
        continue;

      if (!isfinite(prof->qdlat_Ne))
        continue;

      qd_J[n] = prof->qdlat_J;
      qd_Ne[n] = prof->qdlat_Ne;

      peak_J[n] = prof->peak_J;
      peak_Ne[n] = prof->peak_Ne;

      ++n;
    }

  *corr_center = gsl_stats_correlation(qd_J, 1, qd_Ne, 1, n);
  *corr_peak = gsl_stats_correlation(peak_J, 1, peak_Ne, 1, n);

  fprintf(stderr, "docorr: n = %zu\n", n);

  return 0;
}

int
analyze_hemisphere(const double peak_qd_min, const double peak_qd_max,
                   const char *peak_file,
                   const char *prof_file, const satdata_mag *mag_data,
                   const track_workspace *track_p, current_data *data)
{
  fprintf(stderr, "main: searching for northern J peaks in satellite...");
  find_J_peaks(peak_qd_min, peak_qd_max, -100.0, peak_file, prof_file, mag_data, track_p, data);
  fprintf(stderr, "done\n");

  {
    double corr_center, corr_peak;

    docorr(data, &corr_center, &corr_peak);

    fprintf(stderr, "main: north hemisphere corr(QD(J), QD(Ne))     = %f\n",
            corr_center);
    fprintf(stderr, "main: north hemisphere corr(peak(J), peak(Ne)) = %f\n",
            corr_peak);
  }
  
  fprintf(stderr, "main: peak data written to %s\n", peak_file);
  fprintf(stderr, "main: profile data written to %s\n", prof_file);

  return 0;
}

/* find profile which is closest in time to t */
int
find_profile_t(const time_t t, current_data *data, size_t *prof_idx)
{
  int s = -1;
  size_t i;
  time_t dt_min = 1e6;
  int idx = -1;

  for (i = 0; i < data->n; ++i)
    {
      curr_profile *prof = &(data->profiles[i]);
      time_t dt = prof->t - t;

      if (abs(dt) < dt_min)
        {
          dt_min = abs(dt);
          idx = (int) i;
          s = 0;
        }
    }

  *prof_idx = (size_t) idx;

  return s;
}

int
correlateJ(const char *data_file, const char *corr_file, current_data *data1, current_data *data2)
{
  int s = 0;
  size_t i, j;
  bin2d_workspace *bin2d_p;
  bin_workspace *bin_p;
  const double phi_min = -5.0;
  const double phi_max = 25.0;
  const size_t nphi = 7;
  const double t_min = -60.0;
  const double t_max = 60.0;
  const size_t nt = 2;
  FILE *fp_data, *fp_corr;
  double *dp, *x, *y, *r;
  size_t n = 0;
  double sigma;

  fp_data = fopen(data_file, "w");
  fp_corr = fopen(corr_file, "w");

  bin2d_p = bin2d_alloc(phi_min, phi_max, nphi, t_min, t_max, nt);
  bin_p = bin_alloc(phi_min, phi_max, nphi);
  dp = malloc(data1->n * sizeof(double));
  x = malloc(data1->n * sizeof(double));
  y = malloc(data1->n * sizeof(double));
  r = malloc(data1->n * sizeof(double));

  i = 1;
  fprintf(fp_data, "# Field %zu: timestamp of satellite 1 crossing (UT)\n", i++);
  fprintf(fp_data, "# Field %zu: delta t (minutes)\n", i++);
  fprintf(fp_data, "# Field %zu: delta phi (degrees)\n", i++);
  fprintf(fp_data, "# Field %zu: satellite 1 peak current (mA/m)\n", i++);
  fprintf(fp_data, "# Field %zu: satellite 2 peak current (mA/m)\n", i++);

  for (i = 0; i < data1->n; ++i)
    {
      curr_profile *prof = &(data1->profiles[i]);
      curr_profile *prof2;
      time_t dt;
      double dphi, dt_min;
      size_t idx;
      double lt = get_localtime(prof->t, prof->phi * M_PI / 180.0);

      s = find_profile_t(prof->t, data2, &idx);
      if (s)
        continue;

      prof2 = &(data2->profiles[idx]);

      if (prof->npeak != 1 || prof2->npeak != 1)
        continue;

#if 0
      if (prof->peak_J > 20.0)
        continue; /* a couple outliers */
#endif

#if 0
      if (lt > 20.0)
        continue;
#endif

      dt = prof2->t - prof->t;
      dphi = wrap180(prof2->phi - prof->phi);
      dt_min = (double)dt / 60.0;

#if 0
      bin2d_add_element_corr(dphi, dt_min, prof->peak_J, prof2->peak_J, bin2d_p);
#endif

      if (fabs(dt_min) <= 30.0 && dphi >= phi_min && dphi <= phi_max)
        {
          /* store this point */
          x[n] = prof->peak_J;
          y[n] = prof2->peak_J;
          dp[n] = dphi;
          ++n;

          fprintf(fp_data, "%ld %f %f %f %f\n",
                  prof->t,
                  dt_min,
                  dphi,
                  prof->peak_J,
                  prof2->peak_J);
        }
    }

  /* robust fit straight line to (x,y) */
  {
    curvefit_workspace *curvefit_p = curvefit_alloc(gsl_multifit_robust_bisquare, curvefit_poly, n, 2);
    curvefit(0, x, y, curvefit_p);
    curvefit_residuals(x, y, r, curvefit_p);
    curvefit_free(curvefit_p);
    const double alpha = 2.5;

    sigma = gsl_stats_sd(r, 1, n);
    fprintf(stderr, "sigma = %f\n", sigma);

    /* loop through again and add (J1,J2) points which aren't outliers */
    for (i = 0; i < n; ++i)
      {
        if (fabs(r[i]) > alpha*sigma)
          continue;

        bin_add_element_corr(dp[i], x[i], y[i], bin_p);
      }
  }

#if 0
  i = 1;
  fprintf(fp_corr, "# Field %zu: delta longitude (degrees)\n", i++);
  fprintf(fp_corr, "# Field %zu: delta t (minutes)\n", i++);
  fprintf(fp_corr, "# Field %zu: correlation\n", i++);
  fprintf(fp_corr, "# Field %zu: number of points in bin\n", i++);

  for (i = 0; i < nphi; ++i)
    {
      for (j = 0; j < nt; ++j)
        {
          double dt, dlon;
          size_t n;
          double r;

          bin2d_xyval(i, j, &dlon, &dt, bin2d_p);
          n = bin2d_n(dlon, dt, bin2d_p);
          r = bin2d_correlation(dlon, dt, bin2d_p);

          fprintf(fp_corr, "%f %f %f %zu\n",
                  dlon,
                  dt,
                  r,
                  n);
        }

      fprintf(fp_corr, "\n");
    }
#else
  i = 1;
  fprintf(fp_corr, "# Field %zu: delta longitude (degrees)\n", i++);
  fprintf(fp_corr, "# Field %zu: correlation\n", i++);
  fprintf(fp_corr, "# Field %zu: number of points in bin\n", i++);

  for (i = 0; i < nphi; ++i)
    {
      double dlon;
      size_t n;
      double r;

      bin_xval(i, &dlon, bin_p);
      n = bin_n(dlon, bin_p);
      r = bin_correlation(dlon, bin_p);

      fprintf(fp_corr, "%f %f %zu\n",
              dlon,
              r,
              n);
    }
#endif

  fclose(fp_data);
  fclose(fp_corr);
  free(x);
  free(y);
  free(r);

  return s;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --curr_file  | -j line_current_file         - Line current file for satellite 1\n");
  fprintf(stderr, "\t --curr_file2 | -k line_current_file2        - Line current file for satellite 2\n");
  fprintf(stderr, "\t --swarm_file | -s swarm_index_file          - Swarm magnetic index file\n");
  fprintf(stderr, "\t --all | -a                                  - print all tracks (no filtering)\n");
}

int
main(int argc, char *argv[])
{
  const char *peak_file1 = "peak1.dat";
  const char *prof_file1 = "prof1.dat";
  const char *peak_file2 = "peak2.dat";
  const char *prof_file2 = "prof2.dat";
  satdata_mag *data = NULL;
  current_data sat1, sat2;
  peak_workspace *peak_workspace_p;
  struct timeval tv0, tv1;

  sat1.n = 0;
  sat2.n = 0;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "curr_file", required_argument, NULL, 'j' },
          { "curr_file2", required_argument, NULL, 'k' },
          { "swarm_file", required_argument, NULL, 's' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "j:k:s:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'j':
            fprintf(stderr, "main: reading %s...", optarg);
            read_lc(optarg, &sat1);
            fprintf(stderr, "done (%zu profiles read)\n", sat1.n);
            break;

          case 'k':
            fprintf(stderr, "main: reading %s...", optarg);
            read_lc(optarg, &sat2);
            fprintf(stderr, "done (%zu profiles read)\n", sat2.n);
            break;

          case 's':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = satdata_swarm_read_idx(optarg, 0);
            gettimeofday(&tv1, NULL);
            if (!data)
              exit(1);
            fprintf(stderr, "done (%zu points read, %g seconds)\n",
                    data->n, time_diff(tv0, tv1));

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

  peak_workspace_p = peak_alloc(NCURR);

  if (data)
    {
      track_workspace *track_p = track_alloc();

      fprintf(stderr, "main: separating into tracks...");
      track_init(data, NULL, track_p);
      fprintf(stderr, "done\n");

      /* north hemisphere peak finding */
      analyze_hemisphere(-3.0, 23.0, peak_file1, prof_file1, data, track_p, &sat1);

      /* south hemisphere peak finding */
      analyze_hemisphere(-23.0, 0.0, peak_file2, prof_file2, data, track_p, &sat1);

      satdata_mag_free(data);
      track_free(track_p);
    }
  else if (sat2.n > 0)
    {
      const char *data_file1 = "data.dat.north";
      const char *corr_file1 = "corr.dat.north";
      const char *data_file2 = "data.dat.south";
      const char *corr_file2 = "corr.dat.south";

      fprintf(stderr, "main: searching for northern J peaks in satellite 1...");
      find_J_peaks(5.0, 23.0, -100.0, peak_file1, prof_file1, NULL, NULL, &sat1);
      fprintf(stderr, "done\n");

      fprintf(stderr, "main: searching for northern J peaks in satellite 2...");
      find_J_peaks(5.0, 23.0, -100.0, peak_file2, prof_file2, NULL, NULL, &sat2);
      fprintf(stderr, "done\n");

      correlateJ(data_file1, corr_file1, &sat1, &sat2);

      fprintf(stderr, "main: data printed to %s\n", data_file1);
      fprintf(stderr, "main: correlation data printed to %s\n", corr_file1);

      fprintf(stderr, "main: searching for northern J peaks in satellite 2...");
      find_J_peaks(-23.0, -5.0, -100.0, peak_file1, prof_file1, NULL, NULL, &sat1);
      fprintf(stderr, "done\n");

      fprintf(stderr, "main: searching for northern J peaks in satellite 2...");
      find_J_peaks(-23.0, -5.0, -100.0, peak_file2, prof_file2, NULL, NULL, &sat2);
      fprintf(stderr, "done\n");

      correlateJ(data_file2, corr_file2, &sat1, &sat2);

      fprintf(stderr, "main: data printed to %s\n", data_file2);
      fprintf(stderr, "main: correlation data printed to %s\n", corr_file2);
    }

  peak_free(peak_workspace_p);

  return 0;
} /* main() */
