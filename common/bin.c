/*
 * bin.c
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

#include "bin.h"

static size_t bin_findbin(const double val, const double a, const double b, const size_t n);
static double bin_interp1d(const double a, const double b, const double fa,
                           const double fb, const double x);

bin_workspace *
bin_alloc(const double xmin, const double xmax, const size_t nx)
{
  bin_workspace *w;
  size_t i;

  w = calloc(1, sizeof(bin_workspace));
  if (!w)
    return 0;

  w->nbins = nx;
  
  w->bins = malloc(w->nbins * sizeof(gsl_rstat_workspace *));
  if (!w->bins)
    {
      fprintf(stderr, "bin_alloc: malloc failed: %s\n", strerror(errno));
      bin_free(w);
      return 0;
    }

  for (i = 0; i < w->nbins; ++i)
    {
      w->bins[i] = gsl_rstat_alloc();
      if (!w->bins[i])
        {
          fprintf(stderr, "bin_alloc: gsl_rstat_alloc failed\n");
          bin_free(w);
          return 0;
        }
    }

  w->nx = nx;
  w->xmin = xmin;
  w->xmax = xmax;

  return w;
} /* bin_alloc() */

void
bin_free(bin_workspace *w)
{
  size_t i;

  for (i = 0; i < w->nbins; ++i)
    {
      if (w->bins[i])
        gsl_rstat_free(w->bins[i]);
    }

  if (w->bins)
    free(w->bins);

  free(w);
} /* bin_free() */

int
bin_add_element(const double x, const double y, bin_workspace *w)
{
  int s = 0;
  size_t bin;

  if (x < w->xmin || x > w->xmax)
    {
      fprintf(stderr, "bin_add_element: error: x outside allowed range: %f\n", x);
      return -1;
    }

  bin = bin_findbin(x, w->xmin, w->xmax, w->nx);

  gsl_rstat_add(y, w->bins[bin]);

  return s;
} /* bin_add_element() */

double
bin_mean(const double x, bin_workspace *w)
{
  size_t bin = bin_findbin(x, w->xmin, w->xmax, w->nx);
  double mean;

  mean = gsl_rstat_mean(w->bins[bin]);

  return mean;
} /* bin_mean() */

double
bin_sd(const double x, bin_workspace *w)
{
  size_t bin = bin_findbin(x, w->xmin, w->xmax, w->nx);
  double sd;

  sd = gsl_rstat_sd(w->bins[bin]);

  return sd;
} /* bin_sd() */

size_t
bin_n(const double x, bin_workspace *w)
{
  size_t bin = bin_findbin(x, w->xmin, w->xmax, w->nx);
  size_t n;

  n = gsl_rstat_n(w->bins[bin]);

  return n;
} /* bin_n() */

/*
bin_xval()
  Determine x value for a given bin

Inputs: i - x bin (i \in [0,nx-1])
        x - (output) x = xmin + i * dx
        w - workspace

Return: success or error
*/

int
bin_xval(const size_t i, double *x, bin_workspace *w)
{
  int s = 0;
  const double dx = (w->xmax - w->xmin) / (w->nx - 1.0);

  *x = w->xmin + i * dx;

  return s;
} /* bin_xval() */

int
bin_reset(bin_workspace *w)
{
  int s = 0;
  size_t i;

  for (i = 0; i < w->nbins; ++i)
    gsl_rstat_reset(w->bins[i]);

  return s;
} /* bin_reset() */

/*
bin_findbin()
  Determine bin number corresponding to val \in [a, b] where n is
the number of bins using a linear map in x direction
*/

static size_t
bin_findbin(const double val, const double a, const double b, const size_t n)
{
  size_t binno = (size_t) round(bin_interp1d(a, b, 0.0, n - 1.0, val));

  assert(binno < n);

  return binno;
} /* bin_findbin() */

/*
bin_interp1d()
  1D interpolation: find f(x) using endpoints f(a) and f(b)
*/

static double
bin_interp1d(const double a, const double b, const double fa,
             const double fb, const double x)
{
  double fx;

  fx = (b - x) / (b - a) * fa + (x - a) / (b - a) * fb;

  return fx;
} /* bin_interp1d() */
