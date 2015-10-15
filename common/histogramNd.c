/*
 * histogramNd.c
 * Copyright (C) 2006 Patrick Alken
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_errno.h>

#include "histogramNd.h"

static int get_bin_by_x(struct histogramNd_workspace *h, double *x,
                        size_t *bin);
static int find (const size_t n, const double range[], 
                 const double x, size_t * i);
static int findNd(struct histogramNd_workspace *h, double *x);

/*
histogramNd_alloc()
  Allocate a histogramNd workspace

Inputs: n - number of dimensions
        N - array of length n containing number of bins in each
            dimension
*/

struct histogramNd_workspace *
histogramNd_alloc(size_t n, size_t *N)
{
  struct histogramNd_workspace *h;
  size_t i;
  size_t totbins;

  h = (struct histogramNd_workspace *)
      malloc(sizeof(struct histogramNd_workspace));
  if (!h)
    {
      fprintf(stderr, "histogramNd_alloc: malloc failed: %s\n",
              strerror(errno));
      return (0);
    }

  memset(h, '\0', sizeof(struct histogramNd_workspace));

  h->N = (size_t *) malloc(sizeof(size_t) * n);
  if (!h->N)
    {
      fprintf(stderr, "histogramNd_alloc: malloc failed: %s\n",
              strerror(errno));
      histogramNd_free(h);
      return (0);
    }

  totbins = 1;
  for (i = 0; i < n; ++i)
    {
      totbins *= N[i];
      h->N[i] = N[i];
    }

  h->bins = (double *) malloc(sizeof(double) * totbins);
  h->range = (double **) malloc(sizeof(double *) * n);
  h->loc = (size_t *) malloc(sizeof(size_t) * n);
  if (!h->bins || !h->range)
    {
      fprintf(stderr, "histogramNd_alloc: malloc failed: %s\n",
              strerror(errno));
      histogramNd_free(h);
      return (0);
    }

  for (i = 0; i < n; ++i)
    h->range[i] = (double *) malloc(sizeof(double) * (N[i] + 1));

  h->n = n;
  h->total_bins = totbins;

  return (h);
} /* histogramNd_alloc() */

/*
histogramNd_free()
  Free a histogramNd workspace
*/

void
histogramNd_free(struct histogramNd_workspace *h)
{
  if (!h)
    return;

  if (h->N)
    free(h->N);

  if (h->bins)
    free(h->bins);

  if (h->range)
    {
      size_t i;

      for (i = 0; i < h->n; ++i)
        {
          if (h->range[i])
            free(h->range[i]);
        }

      free(h->range);
    }

  if (h->loc)
    free(h->loc);

  free(h);
} /* histogramNd_free() */

void
histogramNd_set_ranges_uniform(struct histogramNd_workspace *h,
                               double *min, double *max)
{
  size_t i, j;

  for (i = 0; i < h->n; ++i)
    {
      size_t Ni = h->N[i];

      for (j = 0; j <= Ni; ++j)
        {
          double f1 = ((double) (Ni - j) / (double) Ni);
          double f2 = ((double) j / (double) Ni);

          h->range[i][j] = f1 * min[i] + f2 * max[i];
        }
    }

  /* initialize bins */
  for (i = 0; i < h->total_bins; ++i)
    h->bins[i] = 0.0;
} /* histogramNd_set_ranges_uniform() */

/*
histogramNd_increment()
  Add 1 to the bin containing 'x'
*/

void
histogramNd_increment(struct histogramNd_workspace *h,
                      double *x)
{
  histogramNd_accumulate(h, x, 1.0);
} /* histogramNd_increment() */

/*
histogramNd_accumulate()
*/

void
histogramNd_accumulate(struct histogramNd_workspace *h,
                      double *x, double weight)
{
  size_t bin;
  int status;

  status = get_bin_by_x(h, x, &bin);
  if (status)
    return;

  h->bins[bin] += weight;
} /* histogramNd_accumulate() */

int
histogramNd_get_by_x(struct histogramNd_workspace *h,
                     double *x, double *w)
{
  size_t bin;
  int status;

  status = get_bin_by_x(h, x, &bin);
  if (status)
    return status;

  *w = h->bins[bin];

  return 0;
} /* histogramNd_get_by_x() */

static int
get_bin_by_x(struct histogramNd_workspace *h, double *x, size_t *bin)
{
  size_t i;
  int status;

  status = findNd(h, x);
  if (status)
    return status;

  *bin = h->loc[h->n - 1];
  for (i = 0; i < h->n - 1; ++i)
    {
      *bin += h->loc[i] * h->N[i + 1];
    }

  return 0;
} /* get_bin_by_x() */

static int
findNd(struct histogramNd_workspace *h, double *x)
{
  size_t i, idx;

  for (i = 0; i < h->n; ++i)
    {
      int status = find(h->N[i], h->range[i], x[i], &idx);

      if (status)
        return (status);

      h->loc[i] = idx;
    }

  return 0;
} /* findNd() */

#define LINEAR_OPT 1

static int
find (const size_t n, const double range[], const double x, size_t * i)
{
  size_t i_linear, lower, upper, mid;

  if (x < range[0])
    {
      return -1;
    }

  if (x >= range[n])
    {
      return +1;
    }

  /* optimize for linear case */

#ifdef LINEAR_OPT
  {
    double u =  (x - range[0]) / (range[n] - range[0]);
    i_linear = (size_t) (u * n);
  }

  if (x >= range[i_linear] && x < range[i_linear + 1])
    {
      *i = i_linear;
      return 0;
    }
#endif

  /* perform binary search */

  upper = n ;
  lower = 0 ;

  while (upper - lower > 1)
    {
      mid = (upper + lower) / 2 ;
      
      if (x >= range[mid])
        {
          lower = mid ;
        }
      else
        {
          upper = mid ;
        }
    }

  *i = lower ;

  /* sanity check the result */

  if (x < range[lower] || x >= range[lower + 1])
    {
      GSL_ERROR ("x not found in range", GSL_ESANITY);
    }

  return 0;
}
