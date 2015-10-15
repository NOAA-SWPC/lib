/*
 * bin2d.c
 *
 * Bin a dataset in 2 dimensions
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <gsl/gsl_rstat.h>
#include <gsl/gsl_statistics.h>

#include "bin2d.h"

static size_t bin2d_findbin(double val, double a, double b, size_t n);
static double bin2d_interp1d(double a, double b, double fa, double fb, double x);

bin2d_workspace *
bin2d_alloc(double xmin, double xmax, size_t nx, double ymin, double ymax,
            size_t ny)
{
  bin2d_workspace *w;
  size_t i;

  w = calloc(1, sizeof(bin2d_workspace));
  if (!w)
    return 0;

  w->nbins = nx * ny;
  
  w->bins = malloc(w->nbins * sizeof(gsl_rstat_workspace *));
  if (!w->bins)
    {
      fprintf(stderr, "bin2d_alloc: malloc failed: %s\n", strerror(errno));
      bin2d_free(w);
      return 0;
    }

  w->z1 = malloc(w->nbins * sizeof(double *));
  w->z2 = malloc(w->nbins * sizeof(double *));
  w->n = calloc(1, w->nbins * sizeof(size_t));

  for (i = 0; i < w->nbins; ++i)
    {
      w->bins[i] = gsl_rstat_alloc();
      if (!w->bins[i])
        {
          fprintf(stderr, "bin2d_alloc: gsl_rstat_alloc failed\n");
          bin2d_free(w);
          return 0;
        }

      w->z1[i] = malloc(MAX_DATA_PER_BIN * sizeof(double));
      w->z2[i] = malloc(MAX_DATA_PER_BIN * sizeof(double));
    }

  w->nx = nx;
  w->ny = ny;

  w->xmin = xmin;
  w->xmax = xmax;
  w->ymin = ymin;
  w->ymax = ymax;

  return w;
} /* bin2d_alloc() */

void
bin2d_free(bin2d_workspace *w)
{
  size_t i;

  for (i = 0; i < w->nbins; ++i)
    {
      if (w->bins[i])
        gsl_rstat_free(w->bins[i]);

      if (w->z1[i])
        free(w->z1[i]);

      if (w->z2[i])
        free(w->z2[i]);
    }

  if (w->bins)
    free(w->bins);

  if (w->z1)
    free(w->z1);

  if (w->z2)
    free(w->z2);

  if (w->n)
    free(w->n);

  free(w);
} /* bin2d_free() */

int
bin2d_add_element(double x, double y, double data, bin2d_workspace *w)
{
  int s = 0;
  size_t binx, biny;

  if (x < w->xmin || x > w->xmax)
    {
      fprintf(stderr, "bin_add_element: error: x outside allowed range: %f\n", x);
      return -1;
    }
  else if (y < w->ymin || y > w->ymax)
    {
      fprintf(stderr, "bin_add_element: error: y outside allowed range: %f\n", y);
      return -1;
    }

  binx = bin2d_findbin(x, w->xmin, w->xmax, w->nx);
  biny = bin2d_findbin(y, w->ymin, w->ymax, w->ny);

  gsl_rstat_add(data, w->bins[BIN2D_IDX(binx, biny, w)]);

  return s;
} /* bin2d_add_element() */

int
bin2d_add_element_corr(double x, double y, double data1, double data2, bin2d_workspace *w)
{
  int s = 0;
  size_t binx, biny;
  double *z1, *z2;
  size_t n;

  if (x < w->xmin || x > w->xmax)
    {
      fprintf(stderr, "bin_add_element_corr: error: x outside allowed range: %f\n", x);
      return -1;
    }
  else if (y < w->ymin || y > w->ymax)
    {
      fprintf(stderr, "bin_add_element_corr: error: y outside allowed range: %f\n", y);
      return -1;
    }

  binx = bin2d_findbin(x, w->xmin, w->xmax, w->nx);
  biny = bin2d_findbin(y, w->ymin, w->ymax, w->ny);

  z1 = w->z1[BIN2D_IDX(binx, biny, w)];
  z2 = w->z2[BIN2D_IDX(binx, biny, w)];
  n = w->n[BIN2D_IDX(binx, biny, w)];

  if (n >= MAX_DATA_PER_BIN)
    {
      fprintf(stderr, "bin_add_element_corr: MAX_DATA_PER_BIN too small\n");
      return -1;
    }

  z1[n] = data1;
  z2[n] = data2;

  w->n[BIN2D_IDX(binx, biny, w)] = ++n;

  /* to update 'n' count */
  gsl_rstat_add(data1, w->bins[BIN2D_IDX(binx, biny, w)]);

  return s;
} /* bin2d_add_element() */

double
bin2d_correlation(const double x, const double y, const bin2d_workspace *w)
{
  size_t binx = bin2d_findbin(x, w->xmin, w->xmax, w->nx);
  size_t biny = bin2d_findbin(y, w->ymin, w->ymax, w->ny);
  double *z1 = w->z1[BIN2D_IDX(binx, biny, w)];
  double *z2 = w->z2[BIN2D_IDX(binx, biny, w)];
  size_t n = w->n[BIN2D_IDX(binx, biny, w)];
  double r = 0.0;
  
  if (n > 1)
    r = gsl_stats_correlation(z1, 1, z2, 1, n);

  return r;
} /* bin2d_correlation() */

double
bin2d_mean(const double x, const double y, const bin2d_workspace *w)
{
  size_t binx = bin2d_findbin(x, w->xmin, w->xmax, w->nx);
  size_t biny = bin2d_findbin(y, w->ymin, w->ymax, w->ny);
  double mean;

  mean = gsl_rstat_mean(w->bins[BIN2D_IDX(binx, biny, w)]);

  return mean;
} /* bin2d_mean() */

double
bin2d_sd(const double x, const double y, const bin2d_workspace *w)
{
  size_t binx = bin2d_findbin(x, w->xmin, w->xmax, w->nx);
  size_t biny = bin2d_findbin(y, w->ymin, w->ymax, w->ny);
  double sd;

  sd = gsl_rstat_sd(w->bins[BIN2D_IDX(binx, biny, w)]);

  return sd;
} /* bin2d_sd() */

double
bin2d_median(const double x, const double y, const bin2d_workspace *w)
{
  size_t binx = bin2d_findbin(x, w->xmin, w->xmax, w->nx);
  size_t biny = bin2d_findbin(y, w->ymin, w->ymax, w->ny);
  double median;

  median = gsl_rstat_median(w->bins[BIN2D_IDX(binx, biny, w)]);

  return median;
} /* bin2d_median() */

size_t
bin2d_n(const double x, const double y, const bin2d_workspace *w)
{
  size_t binx = bin2d_findbin(x, w->xmin, w->xmax, w->nx);
  size_t biny = bin2d_findbin(y, w->ymin, w->ymax, w->ny);
  size_t n;

  n = gsl_rstat_n(w->bins[BIN2D_IDX(binx, biny, w)]);

  return n;
} /* bin2d_n() */

/*
bin2d_xyval()
  Determine x and y bin values for a given bin

Inputs: i - x bin (i \in [0,nx-1])
        j - y bin (j \in [0,ny-1])
        x - (output) x = xmin + i * dx
        y - (output) y = ymin + j * dy
        w - workspace

Return: success or error
*/

int
bin2d_xyval(const size_t i, const size_t j,
            double *x, double *y, const bin2d_workspace *w)
{
  int s = 0;
  const double dx = (w->xmax - w->xmin) / (w->nx - 1.0);
  const double dy = (w->ymax - w->ymin) / (w->ny - 1.0);

  *x = w->xmin + i * dx;
  *y = w->ymin + j * dy;

  return s;
} /* bin2d_xyval() */

/* print 2d grid to file */
int
bin2d_print(const char *filename, const bin2d_workspace *w)
{
  int s = 0;
  size_t i, j;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "bin2d_print: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: x\n", i++);
  fprintf(fp, "# Field %zu: y\n", i++);
  fprintf(fp, "# Field %zu: mean\n", i++);
  fprintf(fp, "# Field %zu: median\n", i++);
  fprintf(fp, "# Field %zu: stddev\n", i++);
  fprintf(fp, "# Field %zu: number of data points in bin\n", i++);

  for (i = 0; i < w->nx; ++i)
    {
      for (j = 0; j < w->ny; ++j)
        {
          double x, y;
          double mean, sd, median;
          size_t n;

          bin2d_xyval(i, j, &x, &y, w);

          mean = bin2d_mean(x, y, w);
          sd = bin2d_sd(x, y, w);
          median = bin2d_median(x, y, w);
          n = bin2d_n(x, y, w);

          fprintf(fp, "%f %f %.12e %.12e %.12e %zu\n", x, y, mean, median, sd, n);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

/*
bin2d_findbin()
  Determine bin number corresponding to val \in [a, b] where n is
the number of bins using a linear map separately in x and y
directions
*/

static size_t
bin2d_findbin(double val, double a, double b, size_t n)
{
  size_t binno = (size_t) round(bin2d_interp1d(a, b, 0.0, n - 1.0, val));

  assert(binno < n);

  return binno;
} /* bin2d_findbin() */

/*
bin2d_interp1d()
  1D interpolation: find f(x) using endpoints f(a) and f(b)
*/

static double
bin2d_interp1d(double a, double b, double fa, double fb, double x)
{
  double fx;

  fx = (b - x) / (b - a) * fa + (x - a) / (b - a) * fb;

  return fx;
} /* bin2d_interp1d() */
