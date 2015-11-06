/*
 * curvefit.c
 *
 * Routines for curve fitting of data
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_errno.h>

#include "curvefit.h"

#include "curvefit_bspline.c"
#include "curvefit_poly.c"

/*
curvefit_alloc()
  Allocate a curvefit workspace

Inputs: n - number of data observations
        p - number of fit coefficients
*/

curvefit_workspace *
curvefit_alloc(const gsl_multifit_robust_type *robust_t,
               const curvefit_type *T,
               const size_t n, const size_t p)
{
  curvefit_workspace *w;

  w = calloc(1, sizeof(curvefit_workspace));
  if (!w)
    return 0;

  w->p = p;
  w->n = n;
  w->type = T;

  w->state = malloc(T->size);
  (w->type->alloc)(w->state, n, p);

  w->robust_p = gsl_multifit_robust_alloc(robust_t, n, p);
  w->A = gsl_matrix_alloc(n, p);
  w->f = gsl_vector_alloc(n);
  w->c = gsl_vector_alloc(p);
  w->cov = gsl_matrix_alloc(p, p);
  w->x = gsl_vector_alloc(n);
  w->r = gsl_vector_alloc(n);

  if (w->robust_p == 0)
    {
      curvefit_free(w);
      return 0;
    }

  gsl_multifit_robust_maxiter(1000, w->robust_p);

  w->mean = 0.0;
  w->sigma = 1.0;

  return w;
} /* curvefit_alloc() */

void
curvefit_free(curvefit_workspace *w)
{
  if (w->state)
    {
      (w->type->free)(w->state);
      free(w->state);
    }

  if (w->robust_p)
    gsl_multifit_robust_free(w->robust_p);

  if (w->A)
    gsl_matrix_free(w->A);

  if (w->f)
    gsl_vector_free(w->f);

  if (w->c)
    gsl_vector_free(w->c);

  if (w->cov)
    gsl_matrix_free(w->cov);

  if (w->x)
    gsl_vector_free(w->x);

  if (w->r)
    gsl_vector_free(w->r);

  free(w);
} /* curvefit_free() */

/*
curvefit()
  Perform least-squares fit to data

Inputs: standardize - standardize x data: x <- (x - mean) / stddev
        x           - x data (length w->n)
        y           - y data (length w->n)
        w           - workspace

Return: success or error
*/

int
curvefit(const int standardize, const double *x, const double *y,
         curvefit_workspace *w)
{
  int s = 0;

  /* standardize and initialize least squares system */
  s = curvefit_init(standardize, x, y, w);
  if (s)
    return s;

  /* solve least-squares system */
  s = curvefit_solve(w);

  return s;
} /* curvefit() */

/*
curvefit_init()
  Standardize and initialize least squares system

Inputs: standardize - standardize x data: x <- (x - mean) / stddev
        x           - x data (length w->n)
        y           - y data (length w->n)
        w           - workspace

Return: success or error
*/

int
curvefit_init(const int standardize, const double *x, const double *y,
              curvefit_workspace *w)
{
  int s = 0;
  const size_t n = w->n;
  size_t i;

  w->mean = 0.0;
  w->sigma = 1.0;

  if (standardize)
    {
      w->mean = gsl_stats_mean(x, 1, n);
      w->sigma = gsl_stats_sd(x, 1, n);

      /* don't scale if sigma = 0 */
      if (w->sigma == 0.0)
        w->sigma = 1.0;
    }

  /* center and scale x data */
  for (i = 0; i < n; ++i)
    {
      double xi = (x[i] - w->mean) / w->sigma;
      gsl_vector_set(w->x, i, xi);
    }

  (w->type->init)(w->state, w->x->data, y);

  /* construct design matrix A and rhs vector f */
  for (i = 0; i < w->n; ++i)
    {
      gsl_vector_view v = gsl_matrix_row(w->A, i);
      double xi = gsl_vector_get(w->x, i);

      /* fill in row i of A: A(i,j) = x_i^j */
      (w->type->design_row)(w->state, xi, &v.vector);

      /* fill in RHS vector */
      gsl_vector_set(w->f, i, y[i]);
    }

  return s;
} /* curvefit_init() */

/*
curvefit_solve()
  Solve least squares system

Inputs: w           - workspace

Return: success or error
*/

int
curvefit_solve(curvefit_workspace *w)
{
  int s;
  gsl_error_handler_t *err_handler;

  /*
   * solve least-squares system; turn error handler off so a failure
   * to converge won't exit the program
   */
  err_handler = gsl_set_error_handler_off();
  s = gsl_multifit_robust(w->A, w->f, w->c, w->cov, w->robust_p);
  gsl_set_error_handler(err_handler);

  /* compute residuals */
  s += gsl_multifit_robust_residuals(w->A, w->f, w->c, w->r, w->robust_p);

  return s;
} /* curvefit_solve() */

/*
curvefit_eval()
  Evaluate model at a given point

Inputs: x - point at which to evaluate model
        w - workspace

Return: p(x) using previously computed model coefficients
from curvefit()
*/

double
curvefit_eval(const double x, curvefit_workspace *w)
{
  gsl_vector_view v;
  double p; /* p(x) */
  double xstd;

  /* center/scale x value */
  xstd = (x - w->mean) / w->sigma;

  /* use first row of A as temporary storage for the basis functions */
  v = gsl_matrix_row(w->A, 0);
  (w->type->design_row)(w->state, xstd, &v.vector);

  /* multiply basis functions by model coefficients */
  gsl_blas_ddot(w->c, &v.vector, &p);

  return p;
} /* curvefit_eval() */

double
curvefit_residual(const curvefit_workspace *w)
{
  double rnorm = gsl_blas_dnrm2(w->r);

  return rnorm;
}

int
curvefit_residuals(const double *x, const double *y, double *r,
                   curvefit_workspace *w)
{
  size_t i;

  for (i = 0; i < w->n; ++i)
    {
      double yest = curvefit_eval(x[i], w);
      r[i] = y[i] - yest;
    }

  return GSL_SUCCESS;
}
