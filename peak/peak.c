/*
 * peak.c
 *
 * Routines for finding peaks in time series data. Calling sequence:
 *
 * 1. peak_alloc    - allocate peak workspace
 * 2. peak_init     - initialize by computing and smoothing first derivative
 * 3. peak_find     - find peaks by looking for zero crossings of first derivative
 * 4. peak_gaussian - fit gaussian to detected peak
 * 5. peak_free     - free memory
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include "gaussfit.h"
#include "peak.h"

static int deriv_calc(const gsl_vector *x, gsl_vector *d);
static int peak_smooth(const size_t window, const size_t n,
                       const double *x, double *smooth_x);
static int fit_gaussian(const size_t fit_width, const size_t i, const gsl_vector *x,
                        const gsl_vector *y, peak_workspace *w);

peak_workspace *
peak_alloc(const size_t n)
{
  peak_workspace *w;

  w = calloc(1, sizeof(peak_workspace));
  if (!w)
    return 0;

  w->deriv = gsl_vector_alloc(n);
  w->sderiv = gsl_vector_alloc(n);

  w->n = n;
  w->p = 3;
  w->idx = 0;
  w->pidx = 0;

  w->c = gsl_vector_alloc(w->p);

  return w;
}

void
peak_free(peak_workspace *w)
{
  if (w->deriv)
    gsl_vector_free(w->deriv);

  if (w->sderiv)
    gsl_vector_free(w->sderiv);

  if (w->c)
    gsl_vector_free(w->c);

  free(w);
}

/*
peak_init()
  Initialize peak finding algorithm for given (x,y) data. This
involves computing the first derivative of y and smoothing that
first derivative

Inputs: smooth_window - window size for smoothing (number of samples)
        x             - x data
        y             - y data
        w             - workspace

Return: success/error

Notes:
1) On output, w->deriv contains the first derivative of y, computed
as central differences

2) On output, w->sderiv contains the smoothed first derivative of y
using smooth_window
*/

int
peak_init(const size_t smooth_window, const gsl_vector *x,
          const gsl_vector *y, peak_workspace *w)
{
  if (x->size != w->n)
    {
      GSL_ERROR("x vector does not match workspace", GSL_EBADLEN);
    }
  else if (y->size != w->n)
    {
      GSL_ERROR("y vector does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int s;

      w->idx = 0;

      /* compute first derivative of y */
      s = deriv_calc(y, w->deriv);
      if (s)
        return s;

      /* smooth first derivative of y */
      s = peak_smooth(smooth_window, w->n, w->deriv->data, w->sderiv->data);
      if (s)
        return s;

      return s;
    }
}

/*
peak_find()

Inputs: minmax       - +1 to search for maximum (peaks)
                       -1 to search for minimum (valleys)
                        0 to search for both
        minslope     - minimum slope of first derivative needed
        minheight    - minimum height needed to detect peak; set to
                       0.0 to disable this test
        x            - x data
        y            - y data
        w            - workspace

Return:
GSL_CONTINUE if peak found
GSL_SUCCESS if no more peaks found

Notes:
1) When a peak is found, w->pidx stores the index into x/y of peak location
*/

int
peak_find(const int minmax, const double minslope, const double minheight,
          const gsl_vector *x, const gsl_vector *y, peak_workspace *w)
{
  const size_t n = w->n;

  if (x->size != n)
    {
      GSL_ERROR("x vector does not match workspace", GSL_EBADLEN);
    }
  else if (y->size != n)
    {
      GSL_ERROR("y vector does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int s = GSL_SUCCESS;

      if (w->idx >= n)
        return s; /* no more data to search */

      for ( ; w->idx < n - 1; ++(w->idx))
        {
          double di = gsl_vector_get(w->sderiv, w->idx);
          double dip1 = gsl_vector_get(w->sderiv, w->idx + 1);
          double yi, yip1;

          /* look for zero crossing of first derivative */
          if (di * dip1 > 0.0)
            continue;

          /*
           * for peaks, derivative changes from positive to negative and
           * vice versa for valleys
           */
          if (minmax == 1 && !(di >= 0.0 && dip1 < 0.0))
            continue;

          if (minmax == -1 && !(di <= 0.0 && dip1 > 0.0))
            continue;

          yi = gsl_vector_get(y, w->idx);
          yip1 = gsl_vector_get(y, w->idx + 1);

          /*
           * check if slope of derivative is less than minslope; this
           * test eliminates broad features as peaks
           */
          if (fabs(di - dip1) < minslope * yi)
            continue;

          /* check if peak height exceeds minimum threshold */
          if (minheight != 0.0 && yi < minheight && yip1 < minheight)
            continue;

          /* store peak location */
          w->pidx = w->idx;

          /* update idx for next call */
          ++(w->idx);

          s = GSL_CONTINUE;
          break;
        }

      /* found peak or done searching data */
      return s;
    }
} /* peak_find() */

/*
peak_gaussian()
  Fit a gaussian to time series around a located peak

Inputs: fit_width    - size of window for gaussian fitting
        x            - x data
        y            - y data
        w            - workspace

Return: success/error

Notes:
1) peak_find() must first be called to locate a peak and set w->pidx
to peak location
*/

int
peak_gaussian(const size_t fit_width, const gsl_vector *x, const gsl_vector *y, peak_workspace *w)
{
  const size_t n = w->n;

  if (x->size != n)
    {
      GSL_ERROR("x vector does not match workspace", GSL_EBADLEN);
    }
  else if (y->size != n)
    {
      GSL_ERROR("y vector does not match workspace", GSL_EBADLEN);
    }
  else
    {
      int s;
      double xi = gsl_vector_get(x, w->pidx);
      double yi = gsl_vector_get(y, w->pidx);

      /* initial guess vector */
      gsl_vector_set(w->c, 0, yi);  /* height */
      gsl_vector_set(w->c, 1, xi);  /* position */
      gsl_vector_set(w->c, 2, 5.0); /* stddev */

      /* fit gaussian */
      s = fit_gaussian(fit_width, w->pidx, x, y, w);

      return s;
    }
}

double
peak_eval(const gsl_vector *c, const double x)
{
  return gaussfit_eval(c, x);
}

double
peak_deriv(const size_t i, const peak_workspace *w)
{
  return gsl_vector_get(w->deriv, i);
}

double
peak_sderiv(const size_t i, const peak_workspace *w)
{
  return gsl_vector_get(w->sderiv, i);
}

static int
deriv_calc(const gsl_vector *x, gsl_vector *d)
{
  size_t i;
  const size_t n = x->size;

  gsl_vector_set(d, 0, gsl_vector_get(x, 1) - gsl_vector_get(x, 0));
  gsl_vector_set(d, n - 1, gsl_vector_get(x, n - 1) - gsl_vector_get(x, n - 2));

  for (i = 1; i < n - 1; ++i)
    {
      double xm1 = gsl_vector_get(x, i - 1);
      double xp1 = gsl_vector_get(x, i + 1);
      
      gsl_vector_set(d, i, 0.5 * (xp1 - xm1));
    }

  return GSL_SUCCESS;
} /* deriv_calc() */

static int
peak_smooth(const size_t window, const size_t n, const double *x, double *smooth_x)
{
  const size_t half_window = window / 2;
  gsl_vector_view sv = gsl_vector_view_array(smooth_x, n);
  double sum = 0.0;
  size_t i;

  gsl_vector_set_zero(&sv.vector);

  for (i = 0; i < window; ++i)
    sum += x[i];

  for (i = 0; i < n - window; ++i)
    {
      smooth_x[i + half_window - 1] = sum;
      sum = (sum + x[i + window]) - x[i];
    }

  for (i = n - window; i < n; ++i)
    smooth_x[n - half_window - 1] += x[i];

  gsl_vector_scale(&sv.vector, 1.0 / (double) window);

  return GSL_SUCCESS;
} /* peak_smooth() */

static int
fit_gaussian(const size_t fit_width, const size_t i, const gsl_vector *x,
             const gsl_vector *y, peak_workspace *w)
{
  int s = GSL_SUCCESS;
  const size_t half_width = fit_width / 2;
  size_t idx0, idx1, n;

  if (i < half_width)
    idx0 = 0;
  else
    idx0 = i - half_width;

  if (i + half_width >= x->size)
    idx1 = x->size - 1;
  else
    idx1 = i + half_width;

  n = idx1 - idx0 + 1;

  {
    gsl_vector_const_view xv = gsl_vector_const_subvector(x, idx0, n);
    gsl_vector_const_view yv = gsl_vector_const_subvector(y, idx0, n);
    gaussfit_workspace *gauss_p = gaussfit_alloc(n, w->p);

    /* initial guess */
    gaussfit_init(w->c, gauss_p);

    /* fit gaussian */
    s = gaussfit(xv.vector.data, yv.vector.data, gauss_p);
    
    /* save parameters */
    gsl_vector_memcpy(w->c, gauss_p->c);

    gaussfit_free(gauss_p);
  }

  w->idx0 = idx0;
  w->idx1 = idx1;

  return s;
} /* fit_gaussian() */
