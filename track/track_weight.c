/*
 * track_weight.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>

#include <satdata/satdata.h>

#include <common/common.h>

#include "track_weight.h"

static double matrix_sum(const gsl_matrix *A);

track_weight_workspace *
track_weight_alloc(const size_t ntheta, const size_t nphi)
{
  track_weight_workspace *w;
  const double eps = 1.0e-8;

  w = calloc(1, sizeof(track_weight_workspace));
  if (!w)
    return 0;

  w->nphi = nphi;
  w->ntheta = ntheta;

  w->weight = gsl_matrix_alloc(nphi, ntheta);

  w->hist_p = gsl_histogram2d_alloc(nphi, ntheta);
  if (!w->hist_p)
    {
      track_weight_free(w);
      return 0;
    }

  /* adding eps to the upper bound of longitude will allow phi = pi values */
  gsl_histogram2d_set_ranges_uniform(w->hist_p,
                                     -M_PI, M_PI + eps,
                                     0.0, M_PI + eps);

  return w;
} /* track_weight_alloc() */

void
track_weight_free(track_weight_workspace *w)
{
  if (w->weight)
    gsl_matrix_free(w->weight);

  if (w->hist_p)
    gsl_histogram2d_free(w->hist_p);

  free(w);
} /* track_weight_free() */

int
track_weight_reset(track_weight_workspace *w)
{
  gsl_histogram2d_reset(w->hist_p);
  return 0;
}

/*
track_weight_add_data()
  Add data to running histogram

Inputs: theta - colatitude (radians)
        phi   - longitude (radians)
        w     - workspace
*/

int
track_weight_add_data(const double theta, const double phi, track_weight_workspace *w)
{
  int s = 0;

  gsl_histogram2d_increment(w->hist_p, phi, theta);

  return s;
} /* track_weight_add_data() */

/*
track_weight_calc()
  Compute weights according to:

w_ij = sqrt(a_ij / (A n_ij))

where a_ij is the area of bin (i,j), A is the total area, and
n_ij is the number of data points in bin (i,j)

The idea is that for regions which are sampled more frequently
(ie: polar regions), we assign smaller weights and less
sampled regions (ie: equator) we assign larger weights. Dividing
by total area helps make sum_ij w_ij of order unity
*/

int
track_weight_calc(track_weight_workspace *w)
{
  int s = 0;
  size_t i, j;

  gsl_matrix_set_zero(w->weight);

  /* calculate total surface area = sum_{ij} dS_{ij} and final weights */
  for (i = 0; i < w->nphi; ++i)
    {
      double phi0, phi1, dphi;

      gsl_histogram2d_get_xrange(w->hist_p, i, &phi0, &phi1);
      dphi = phi1 - phi0;
      assert(dphi > 0.0);

      for (j = 0; j < w->ntheta; ++j)
        {
          double npts = gsl_histogram2d_get(w->hist_p, i, j);

          /* if no data points in this bin, leave weight at 0 */
          if (npts > 0.0)
            {
              double theta0, theta1, theta, dtheta;
              double dS; /* area of this bin */
              double weight = 0.0;

              gsl_histogram2d_get_yrange(w->hist_p, j, &theta0, &theta1);
              dtheta = theta1 - theta0;
              assert(dtheta > 0.0);

              theta = 0.5 * (theta0 + theta1);

              /* compute area of bin on unit sphere */
              dS = sin(theta) * dtheta * dphi;

              weight = sqrt(dS / npts);
              gsl_matrix_set(w->weight, i, j, weight);
            }
        }
    }

#if 0
  /* divide by sum(W) to make sum of weights equal to unity */
  {
    double sum = matrix_sum(w->weight);
    gsl_matrix_scale(w->weight, 1.0 / sum);
    fprintf(stderr, "\n\t sum(W) before = %f\n", sum);
    fprintf(stderr, "\t sum(W) after = %f\n", matrix_sum(w->weight));
  }
#else
  /*
   * divide weight matrix by (10*stddev(weights)) to try to make each
   * weight of order unity
   */
  {
    double sum = matrix_sum(w->weight);
    double sd = gsl_stats_sd(w->weight->data, 1, w->weight->size1 * w->weight->size2);
    gsl_matrix_scale(w->weight, 1.0 / (10.0 * sd));
  }
#endif

  return s;
} /* track_weight_calc() */

int
track_weight_get(const double phi, const double theta, double *weight, track_weight_workspace *w)
{
  int s = 0;
  size_t i, j;

  s = gsl_histogram2d_find(w->hist_p, phi, theta, &i, &j);
  if (s)
    {
      fprintf(stderr, "track_weight_get: error in gsl_histogram2d_find: %d\n", s);
      return s;
    }

  *weight = gsl_matrix_get(w->weight, i, j);

  return s;
} /* track_weight_get() */

int
track_weight_n(const double phi, const double theta, size_t *n, track_weight_workspace *w)
{
  int s = 0;
  size_t i, j;

  s = gsl_histogram2d_find(w->hist_p, phi, theta, &i, &j);
  if (s)
    {
      fprintf(stderr, "track_weight_n: error in gsl_histogram2d_find: %d\n", s);
      return s;
    }

  *n = (size_t) gsl_histogram2d_get(w->hist_p, i, j);

  return s;
} /* track_weight_n() */

int
track_weight_write(const char *filename, track_weight_workspace *w)
{
  int s = 0;
  FILE *fp = fopen(filename, "w");

  if (!fp)
    {
      fprintf(stderr, "track_weight_write: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  fwrite(&(w->nphi), sizeof(size_t), 1, fp);
  fwrite(&(w->ntheta), sizeof(size_t), 1, fp);
  gsl_matrix_fwrite(fp, w->weight);

  fclose(fp);

  return s;
} /* track_weight_write() */

static double
matrix_sum(const gsl_matrix *A)
{
  size_t i, j;
  double sum = 0.0;

  for (i = 0; i < A->size1; ++i)
    {
      for (j = 0; j < A->size2; ++j)
        sum += gsl_matrix_get(A, i, j);
    }

  return sum;
}
