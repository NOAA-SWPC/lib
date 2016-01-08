/*
 * test.c
 */

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "gaussfit.h"

void
test_gaussfit(void)
{
  int s;
  const size_t n = 100;
  const size_t p = 6;
  const double a0 = 25.1;  /* amplitude */
  const double a1 = -2.3; /* position */
  const double a2 = 2.3;  /* stddev */
  const double a3 = -7.3; /* offset */
  const double a4 = 1.1;  /* linear term */
  const double a5 = 0.01; /* quadratic term */
  const double xmin = -15.0;
  const double xmax = 10.0;
  const double dx = (xmax - xmin) / (n - 1.0);
  size_t i;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  double *x = malloc(n * sizeof(double));
  double *y = malloc(n * sizeof(double));
  gaussfit_workspace *gauss_p = gaussfit_alloc(n, p);

  for (i = 0; i < n; ++i)
    {
      double xi = xmin + i * dx;
      double zi = (xi - a1) / a2;
      double yi = a0 * exp(-0.5 * zi * zi) + a3 + a4*xi + a5*xi*xi;
      double ei = gsl_ran_gaussian(r, 0.9);

      x[i] = xi;
      y[i] = yi + ei;
    }

  /* insert outliers */
  /*y[n/2] = 10000.0;
  y[n/4] = 500.0;
  y[n/3] = 800.0;*/

  s = gaussfit(x, y, gauss_p);
  if (s != GSL_SUCCESS)
    fprintf(stderr, "test_gaussfit: failed to converge\n");

  for (i = 0; i < n; ++i)
    {
      printf("%f %f %f\n",
             x[i],
             y[i],
             gaussfit_eval(gauss_p->c, x[i]));
    }

  fprintf(stderr, "test_gaussfit: amplitude = %f [%f]\n", gsl_vector_get(gauss_p->c, 0), a0);
  fprintf(stderr, "test_gaussfit: position  = %f [%f]\n", gsl_vector_get(gauss_p->c, 1), a1);
  fprintf(stderr, "test_gaussfit: stddev    = %f [%f]\n", gsl_vector_get(gauss_p->c, 2), a2);
  if (p > 3)
    fprintf(stderr, "test_gaussfit: constant  = %f [%f]\n", gsl_vector_get(gauss_p->c, 3), a3);
  if (p > 4)
    fprintf(stderr, "test_gaussfit: linear    = %f [%f]\n", gsl_vector_get(gauss_p->c, 4), a4);
  if (p > 5)
    fprintf(stderr, "test_gaussfit: quadratic = %f [%f]\n", gsl_vector_get(gauss_p->c, 5), a5);

  fprintf(stderr, "test_gaussfit: initial cost = %.12e\n", gauss_p->chisq0);
  fprintf(stderr, "test_gaussfit: final cost   = %.12e\n", gauss_p->chisq);

  gsl_rng_free(r);
  free(x);
  free(y);
  gaussfit_free(gauss_p);
}

int
main()
{
  test_gaussfit();

  exit (gsl_test_summary());
}
