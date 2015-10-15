/*
 * ccorr.c
 *
 * Bin two datasets for cross-correlation
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>

#include "ccorr.h"

static size_t ccorr_findbin(const double val, const double a, const double b, const size_t n);
static double ccorr_interp1d(const double a, const double b, const double fa,
                           const double fb, const double x);

ccorr_workspace *
ccorr_alloc(const double xmin, const double xmax, const size_t nx)
{
  ccorr_workspace *w;

  w = calloc(1, sizeof(ccorr_workspace));
  if (!w)
    return 0;

  w->nbins = nx;
  
  w->bins = calloc(1, w->nbins * sizeof(ccorr_bin));
  if (!w->bins)
    {
      fprintf(stderr, "ccorr_alloc: malloc failed: %s\n", strerror(errno));
      ccorr_free(w);
      return 0;
    }

  w->nx = nx;
  w->xmin = xmin;
  w->xmax = xmax;

  return w;
} /* ccorr_alloc() */

void
ccorr_free(ccorr_workspace *w)
{
  if (w->bins)
    free(w->bins);

  free(w);
} /* ccorr_free() */

int
ccorr_add(const double x, const double u, const double v, ccorr_workspace *w)
{
  int s = 0;
  size_t bin;
  ccorr_bin *cbin;

  if (x < w->xmin || x > w->xmax)
    {
      fprintf(stderr, "ccorr_add_element: error: x outside allowed range: %f\n", x);
      return -1;
    }

  bin = ccorr_findbin(x, w->xmin, w->xmax, w->nx);

  cbin = &(w->bins[bin]);
  cbin->u[cbin->n] = u;
  cbin->v[cbin->n] = v;

  if (++(cbin->n) >= CCORR_MAX)
    {
      fprintf(stderr, "ccorr_add: CCORR_MAX not large enough\n");
      s = -1;
    }

  return s;
} /* ccorr_add() */

double
ccorr_r(const double x, ccorr_workspace *w)
{
  size_t bin = ccorr_findbin(x, w->xmin, w->xmax, w->nx);
  size_t n = w->bins[bin].n;
  double r;

  if (n < 2)
    return 0.0;

  r = gsl_stats_correlation(w->bins[bin].u, 1, w->bins[bin].v, 1, n);

  return r;
} /* ccorr_r() */

size_t
ccorr_n(const double x, ccorr_workspace *w)
{
  size_t bin = ccorr_findbin(x, w->xmin, w->xmax, w->nx);

  return w->bins[bin].n;
} /* ccorr_n() */

/*
ccorr_xval()
  Determine x value for a given bin

Inputs: i - x bin (i \in [0,nx-1])
        x - (output) x = xmin + i * dx
        w - workspace

Return: success or error
*/

int
ccorr_xval(const size_t i, double *x, ccorr_workspace *w)
{
  int s = 0;
  const double dx = (w->xmax - w->xmin) / (w->nx - 1.0);

  *x = w->xmin + i * dx;

  return s;
} /* ccorr_xval() */

int
ccorr_reset(ccorr_workspace *w)
{
  int s = 0;
  size_t i;

  for (i = 0; i < w->nbins; ++i)
    w->bins[i].n = 0;

  return s;
} /* ccorr_reset() */

/*
ccorr_findbin()
  Determine bin number corresponding to val \in [a, b] where n is
the number of bins using a linear map in x direction
*/

static size_t
ccorr_findbin(const double val, const double a, const double b, const size_t n)
{
  size_t binno = (size_t) round(ccorr_interp1d(a, b, 0.0, n - 1.0, val));

  assert(binno < n);

  return binno;
} /* ccorr_findbin() */

/*
ccorr_interp1d()
  1D interpolation: find f(x) using endpoints f(a) and f(b)
*/

static double
ccorr_interp1d(const double a, const double b, const double fa,
             const double fb, const double x)
{
  double fx;

  fx = (b - x) / (b - a) * fa + (x - a) / (b - a) * fb;

  return fx;
} /* ccorr_interp1d() */
