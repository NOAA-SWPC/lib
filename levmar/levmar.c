/*
 * levmar.c
 * wrapper for the C levmar library
 */

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>

#include "levmar.h"

static void levmar_f(double *p, double *hx, int m, int n, void *adata);
static void levmar_df(double *p, double *jac, int m, int n, void *params);

levmar_workspace *
levmar_alloc(const size_t n, const size_t p)
{
  levmar_workspace *w;

  w = calloc(1, sizeof(levmar_workspace));
  if (!w)
    return 0;

  w->n = n;
  w->p = p;

  w->rhs = gsl_vector_calloc(n); /* set to 0 */
  w->covar = gsl_matrix_alloc(p, p);

  w->opts[0] = LM_INIT_MU;
  w->opts[1] = 1.0e-15;
  w->opts[2] = 1.0e-15;
  w->opts[3] = 1.0e-15;
  /*w->opts[3] = 1.e-20;*/
  w->opts[4] = LM_DIFF_DELTA;

  w->maxit = 100;

  return w;
} /* levmar_alloc() */

void
levmar_free(levmar_workspace *w)
{
  if (w->rhs)
    gsl_vector_free(w->rhs);

  free(w);
}

/*
levmar_proc

Inputs: x     - (input/output) vector of model parameters
        fdf   - gsl struct containing f,df,fdf
        w     - workspace

Notes:
1) On output, w->covar contains the parameter covariance matrix
*/

int
levmar_proc(gsl_vector *x, gsl_multifit_function_fdf *fdf, levmar_workspace *w)
{
  int s = 0;
  int n = (int) w->n;
  int m = (int) w->p;
  levmar_params params;

  params.fdf = fdf;

  if (fdf->df == NULL)
    {
      dlevmar_dif(&levmar_f,
                  x->data,
                  w->rhs->data,
                  m,
                  n,
                  w->maxit,
                  w->opts,
                  w->info,
                  NULL,
                  w->covar->data,
                  &params);
    }
  else
    {
      dlevmar_der(&levmar_f,
                  &levmar_df,
                  x->data,
                  w->rhs->data,
                  m,
                  n,
                  w->maxit,
                  w->opts,
                  w->info,
                  NULL,
                  w->covar->data,
                  &params);
    }

  if (w->info[6] > 2.0)
    s = (int) w->info[6];

  return s;
}

/*
levmar_bc_proc()
  Compute levmar solution with box constraints

Inputs: x     - (input/output) vector of model parameters
        fdf   - gsl struct containing f,df,fdf
        lb    - vector of lower bounds (size p)
        ub    - vector of upper bounds (size p)
        w     - workspace

Notes:
1) On output, w->covar contains the parameter covariance matrix
*/

int
levmar_bc_proc(gsl_vector *x, gsl_multifit_function_fdf *fdf,
               gsl_vector *lb, gsl_vector *ub, levmar_workspace *w)
{
  int s = 0;
  int n = (int) w->n;
  int m = (int) w->p;
  levmar_params params;

  params.fdf = fdf;

  if (fdf->df == NULL)
    {
      dlevmar_bc_dif(&levmar_f,
                     x->data,
                     w->rhs->data,
                     m,
                     n,
                     lb->data,
                     ub->data,
                     NULL,
                     w->maxit,
                     w->opts,
                     w->info,
                     NULL,
                     w->covar->data,
                     &params);
    }
  else
    {
      dlevmar_bc_der(&levmar_f,
                     &levmar_df,
                     x->data,
                     w->rhs->data,
                     m,
                     n,
                     lb->data,
                     ub->data,
                     NULL,
                     w->maxit,
                     w->opts,
                     w->info,
                     NULL,
                     w->covar->data,
                     &params);
    }

  if (w->info[6] > 2.0)
    s = (int) w->info[6];

  return s;
}

/* final chi^2 */
double
levmar_chisq(const levmar_workspace *w)
{
  return w->info[1];
}

/* initial chi^2 */
double
levmar_chisq0(const levmar_workspace *w)
{
  return w->info[0];
}

size_t
levmar_niter(const levmar_workspace *w)
{
  return (size_t) w->info[5];
}

static void
levmar_f(double *p, double *x, int m, int n, void *params)
{
  levmar_params *lp = (levmar_params *) params;
  gsl_multifit_function_fdf *fdf = lp->fdf;
  gsl_vector_view pv = gsl_vector_view_array(p, m);
  gsl_vector_view xv = gsl_vector_view_array(x, n);

  gsl_multifit_eval_wf(fdf, &pv.vector, NULL, &xv.vector);
} /* levmar_f() */

static void
levmar_df(double *p, double *jac, int m, int n, void *params)
{
  levmar_params *lp = (levmar_params *) params;
  gsl_multifit_function_fdf *fdf = lp->fdf;
  gsl_vector_view pv = gsl_vector_view_array(p, m);
  gsl_matrix_view J = gsl_matrix_view_array(jac, n, m);

  gsl_multifit_eval_wdf(fdf, &pv.vector, NULL, &J.matrix);
}
