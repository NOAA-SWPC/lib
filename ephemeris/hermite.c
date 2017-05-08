/*
 * hermite.c
 * 
 * This module contains routines for piecewise cubic Hermite
 * interpolation, following closely the following reference:
 *
 * [1] Burden and Faires, Numerical Analysis, 9th ed
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_poly.h>

#include "hermite.h"

static int hermite_calc_coeffs(const double xa[], const double ya[],
                               const double dya[], const double x,
                               gsl_interp_accel *acc, hermite_workspace *w);

hermite_workspace *
hermite_alloc(const size_t n, size_t degree)
{
  hermite_workspace *w;
  const size_t max_degree = 2 * n - 1;

  w = malloc(sizeof(hermite_workspace));
  if (!w)
    return 0;

  if (degree == 0)
    degree = max_degree; /* use maximum possible degree */

  w->n = n;

  w->degree = degree;
  w->ncoeff = degree + 1;
  w->max_degree = max_degree;

  /* piecewise interpolation window */
  w->npts = (degree - 1) / 2 + 1;

  w->z = malloc(w->ncoeff * sizeof(double));
  w->q = malloc(w->ncoeff * sizeof(double));
  w->coeff = malloc(w->ncoeff * sizeof(double));
  w->work = malloc(w->ncoeff * sizeof(double));

  return w;
}

void
hermite_free(hermite_workspace *w)
{
  if (w->z)
    free(w->z);

  if (w->q)
    free(w->q);

  if (w->coeff)
    free(w->coeff);

  if (w->work)
    free(w->work);

  free(w);
} /* hermite_free() */

/*
hermite_init()
  Initialize Hermite polynomial interpolator

If we are using the maximum degree polynomial (degree 2*n - 1),
then the coefficients only need to be computed once for each
dataset and this is done here. For piecewise cubic interpolation
the coefficients need to be computed for each interpolated
point x, and that is done in the evaluation functions.
*/

int
hermite_init(const double xa[], const double ya[],
             const double dya[], const size_t size,
             hermite_workspace *w)
{
  int status = GSL_SUCCESS;

  /*
   * If we're fitting the maximum degree polynomial, we
   * only need to calculate the divided difference representation
   * once. For piecewise cubic Hermite polynomials we need to
   * compute it for each interpolated value in the _eval functions
   */
  if (w->degree == w->max_degree)
    {
      status = gsl_poly_dd_hermite_init(w->q, w->z, xa, ya, dya, size);
    }

  return status;
} /* hermite_init() */

double
hermite_eval(const double xa[], const double ya[], 
             const double dya[], double x,
             gsl_interp_accel *acc,
             hermite_workspace *w)
{
  double yp;

  if (w->degree != w->max_degree)
    {
      /* calculate cubic Hermite coefficients for this x */
      hermite_calc_coeffs(xa, ya, dya, x, acc, w);
    }

  /* now evaluate polynomial at the point x */
  yp = gsl_poly_dd_eval(w->q, w->z, w->ncoeff, x);

  return yp;
} /* hermite_eval() */

double
hermite_eval_deriv(const double xa[], const double ya[], 
                   const double dya[], double x,
                   gsl_interp_accel *acc,
                   hermite_workspace *w)
{
  double dyp;

  if (w->degree != w->max_degree)
    {
      /* calculate cubic Hermite coefficients for this x */
      hermite_calc_coeffs(xa, ya, dya, x, acc, w);
    }

  /* now evaluate polynomial derivative at the point x */
  gsl_poly_dd_taylor(w->coeff, x, w->q, w->z, w->ncoeff, w->work);
  dyp = w->coeff[1];

  return dyp;
} /* hermite_eval_deriv() */

double
hermite_eval_deriv2(const double xa[], const double ya[], 
                    const double dya[], double x,
                    gsl_interp_accel *acc,
                    hermite_workspace *w)
{
  double dyp;

  if (w->degree != w->max_degree)
    {
      /* calculate cubic Hermite coefficients for this x */
      hermite_calc_coeffs(xa, ya, dya, x, acc, w);
    }

  /* now evaluate polynomial derivative at the point x */
  gsl_poly_dd_taylor(w->coeff, x, w->q, w->z, w->ncoeff, w->work);
  dyp = 2.0 * w->coeff[2];

  return dyp;
} /* hermite_eval_deriv2() */

/*
hermite_calc_coeffs()
  This function is called for the case of piecewise cubic
Hermite interpolation. First locate the neighboring points
to the desired point x, then construct coefficients of
a cubic Hermite polynomial matching function values and
derivatives at grid points

Inputs: xa  - x array data
        ya  - y array data
        dya - dy/dx array data
        x   - desired interpolation point
        acc - accelerator
        w   - workspace

Notes:

1) On output, the array w->q contains the divided difference
polynomial coefficients, and w->z contains the z array values
(see [1])

2) only cubic interpolation is currently supported
(using the two neighboring grid points surrounding the point x).
Higher order polynomials could easily be implemented by including
more grid points around the point x, but this is not currently
done to keep the code cleaner.
*/

static int
hermite_calc_coeffs(const double xa[], const double ya[],
                    const double dya[], const double x,
                    gsl_interp_accel *acc,
                    hermite_workspace *w)
{
  int s;
  size_t idx;

  /* find idx so that xa[idx] <= x < xa[idx+1] */
  idx = gsl_interp_accel_find(acc, xa, w->n, x);
  assert((xa[idx] <= x && x < xa[idx + 1]) || idx == w->n - 2);

  /* compute cubic Hermite polynomial coefficients */
  s = gsl_poly_dd_hermite_init(w->q, w->z, &xa[idx], &ya[idx], &dya[idx], w->npts);

  return s;
} /* hermite_calc_coeffs() */
