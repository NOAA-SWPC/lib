/*
 * gaussfit.c
 *
 * The model is:
 *
 * f(t) = A0*exp(-z^2/2) + A3 + A4*t + A5*t^2
 *
 * with
 *
 * z = (t - A1) / A2
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_errno.h>

#include "gaussfit.h"

static int gaussfit_f(const gsl_vector *x, void *datap, gsl_vector *f);
static int gaussfit_df(const gsl_vector *x, void *datap, gsl_matrix *J);
static int gaussfit_fvv(const gsl_vector *x, const gsl_vector *v,
                        void *datap, gsl_vector *fvv);

/*
gaussfit_alloc()
  Allocate a gaussfit workspace

Inputs: n - number of data observations
        p - number of model parameters (between 3 and 6)
*/

gaussfit_workspace *
gaussfit_alloc(const size_t n, const size_t p)
{
  gaussfit_workspace *w;
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();

  w = calloc(1, sizeof(gaussfit_workspace));
  if (!w)
    return 0;

  w->n = n;
  w->p = p;

  fdf_params.accel = 1;
  w->nlinear_workspace_p = gsl_multifit_nlinear_alloc(T, &fdf_params, n, p);

  w->c = gsl_vector_calloc(p);

  /* initial guess */
  gsl_vector_set(w->c, 0, 1.0);
  gsl_vector_set(w->c, 2, 1.0);

  return w;
} /* gaussfit_alloc() */

void
gaussfit_free(gaussfit_workspace *w)
{
  if (w->c)
    gsl_vector_free(w->c);

  if (w->nlinear_workspace_p)
    gsl_multifit_nlinear_free(w->nlinear_workspace_p);

  free(w);
} /* gaussfit_free() */

/*
gaussfit_init()
  Provide initial guess to gaussian fit

Inputs: x - vector of initial parameter values
        w - workspace
*/

int
gaussfit_init(const gsl_vector *x, gaussfit_workspace *w)
{
  gsl_vector_memcpy(w->c, x);
  return GSL_SUCCESS;
}

/*
gaussfit()
  Perform least-squares fit to data

Inputs: t           - t data (length w->n)
        y           - y data (length w->n)
        w           - workspace

Return: success or error

Notes:
1) On output, w->c contains model coefficients
*/

int
gaussfit(const double *t, const double *y,
         gaussfit_workspace *w)
{
  int s = 0;
  gsl_multifit_nlinear_fdf fdf;
  gaussfit_data data;
  const double xtol = 1.0e-8;
  const double gtol = 1.0e-8;
  const double ftol = 0.0;
  int info;
  gsl_vector *f = gsl_multifit_nlinear_residual(w->nlinear_workspace_p);
  gsl_vector *c = gsl_multifit_nlinear_position(w->nlinear_workspace_p);

  data.t = (double *) t;
  data.y = (double *) y;
  data.w = w;

  fdf.f = gaussfit_f;
  fdf.df = gaussfit_df;
  fdf.fvv = gaussfit_fvv;
  fdf.n = w->n;
  fdf.p = w->p;
  fdf.params = &data;

  s = gsl_multifit_nlinear_init(w->c, &fdf, w->nlinear_workspace_p);
  if (s)
    return s;

  /* initial cost */
  gsl_blas_ddot(f, f, &(w->chisq0));

  s = gsl_multifit_nlinear_driver(100, xtol, gtol, ftol, NULL, NULL,
                                  &info, w->nlinear_workspace_p);

  /* save model parameters */
  gsl_vector_memcpy(w->c, c);

  /* final cost */
  gsl_blas_ddot(f, f, &(w->chisq));

  return s;
} /* gaussfit() */

/*
gaussfit_eval()
  Evaluate model at a given point

Inputs: x - model coefficient vector
        t - point at which to evaluate model

Return: f(t) using previously computed model coefficients
from gaussfit()
*/

double
gaussfit_eval(const gsl_vector *x, const double t)
{
  const size_t p = x->size;
  double a[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  size_t i;
  double z, f;

  /* fill in model coefficients */
  for (i = 0; i < p; ++i)
    a[i] = gsl_vector_get(x, i);

  z = (t - a[1]) / a[2];
  f = a[0] * exp(-0.5 * z * z) + a[3] + a[4]*t + a[5]*t*t;

  return f;
} /* gaussfit_eval() */

static int
gaussfit_f(const gsl_vector *x, void *datap, gsl_vector *f)
{
  gaussfit_data *data = (gaussfit_data *) datap;
  gaussfit_workspace *w = data->w;
  size_t i;

  for (i = 0; i < w->n; ++i)
    {
      double ti = data->t[i];
      double yi = data->y[i];
      double Yi = gaussfit_eval(x, ti);

      gsl_vector_set(f, i, Yi - yi);
    }

  return GSL_SUCCESS;
} /* gaussfit_f() */

static int
gaussfit_df(const gsl_vector *x, void *datap, gsl_matrix *J)
{
  gaussfit_data *data = (gaussfit_data *) datap;
  gaussfit_workspace *w = data->w;
  const size_t p = w->p;
  size_t i;
  double a0 = gsl_vector_get(x, 0);
  double a1 = gsl_vector_get(x, 1);
  double a2 = gsl_vector_get(x, 2);

  for (i = 0; i < w->n; ++i)
    {
      double ti = data->t[i];
      double zi = (ti - a1) / a2;
      double expterm = exp(-0.5 * zi * zi);
      
      gsl_matrix_set(J, i, 0, expterm);
      gsl_matrix_set(J, i, 1, (a0 / a2) * zi * expterm);
      gsl_matrix_set(J, i, 2, (a0 / a2) * zi * zi * expterm);

      if (p > 3)
        gsl_matrix_set(J, i, 3, 1.0);

      if (p > 4)
        gsl_matrix_set(J, i, 4, ti);

      if (p > 5)
        gsl_matrix_set(J, i, 5, ti * ti);
    }

  return GSL_SUCCESS;
} /* gaussfit_df() */

static int
gaussfit_fvv(const gsl_vector *x, const gsl_vector *v,
             void *datap, gsl_vector *fvv)
{
  gaussfit_data *data = (gaussfit_data *) datap;
  gaussfit_workspace *w = data->w;
  size_t i;
  double A0 = gsl_vector_get(x, 0);
  double A1 = gsl_vector_get(x, 1);
  double A2 = gsl_vector_get(x, 2);
  double v0 = gsl_vector_get(v, 0);
  double v1 = gsl_vector_get(v, 1);
  double v2 = gsl_vector_get(v, 2);

  for (i = 0; i < w->n; ++i)
    {
      double ti = data->t[i];
      double zi = (ti - A1) / A2;
      double e = exp(-0.5 * zi * zi);
      double f01 = (zi / A2) * e;
      double f02 = zi * f01;
      double f11 = (A0 / A2) * (zi*zi - 1.0) / A2 * e;
      double f12 = (A0 / A2) * (zi / A2) * (zi*zi - 2.0) * e;
      double f22 = A0 * (zi / A2) * (zi / A2) * (zi*zi - 3.0) * e;
      double fvvi = 2.0 * f01 * v0 * v1 +
                    2.0 * f02 * v0 * v2 +
                    2.0 * f12 * v1 * v2 +
                          f11 * v1 * v1 +
                          f22 * v2 * v2;

      gsl_vector_set(fvv, i, fvvi);
    }

  return GSL_SUCCESS;
}
