#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "levmar.h"
#include "oct.h"

#define N               400

#define USE_BC          0

struct data {
  size_t n;
  double * y;
  double * sigma;
};

int
expb_f (const gsl_vector * x, void *data, 
        gsl_vector * f)
{
  size_t n = ((struct data *)data)->n;
  double *y = ((struct data *)data)->y;
  double *sigma = ((struct data *) data)->sigma;

  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);
  double b = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Model Yi = A * exp(-lambda * i) + b */
      double t = i;
      double Yi = A * exp (-lambda * t) + b;
      gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
    }

  return GSL_SUCCESS;
}

int
expb_df (const gsl_vector * x, void *data, 
         gsl_matrix * J)
{
  size_t n = ((struct data *)data)->n;
  double *sigma = ((struct data *) data)->sigma;

  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * exp(-lambda * i) + b  */
      /* and the xj are the parameters (A,lambda,b) */
      double t = i;
      double s = sigma[i];
      double e = exp(-lambda * t);
      gsl_matrix_set (J, i, 0, e/s); 
      gsl_matrix_set (J, i, 1, -t * A * e/s);
      gsl_matrix_set (J, i, 2, 1/s);
    }
  return GSL_SUCCESS;
}

int
main (void)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int i, iter = 0;
  const size_t n = N;
  const size_t p = 3;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  double y[N], sigma[N];
  struct data d = { n, y, sigma};
  gsl_multifit_function_fdf f;
  double x_init[3] = { 1.0, 0.0, 0.0 };
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  const gsl_rng_type * type;
  gsl_rng * r;

  gsl_rng_env_setup();

  type = gsl_rng_default;
  r = gsl_rng_alloc (type);

  f.f = &expb_f;
  f.df = &expb_df;
  f.n = n;
  f.p = p;
  f.params = &d;

  /* This is the data to be fitted */

  for (i = 0; i < n; i++)
    {
      double t = i;
      y[i] = 1.0 + 5 * exp (-0.1 * t) 
                 + gsl_ran_gaussian (r, 0.1);
      sigma[i] = 0.1;
    };

  {
    levmar_workspace *lm_workspace_p = levmar_alloc(n, p);
    gsl_vector *x_lev = gsl_vector_alloc(p);
    gsl_matrix *covar_lm = lm_workspace_p->covar;
    double lb_data[3] = { -10.0, 0.0, -10.0 };
    double ub_data[3] = { 10.0, 1.0, 10.0 };
    gsl_vector_view lb = gsl_vector_view_array(lb_data, 3);
    gsl_vector_view ub = gsl_vector_view_array(ub_data, 3);

    f.df = NULL;
    gsl_vector_memcpy(x_lev, &x.vector);

#if USE_BC
    levmar_bc_proc(x_lev, &f, &lb.vector, &ub.vector, lm_workspace_p);
#else
    levmar_proc(x_lev, &f, lm_workspace_p);
#endif

#define ERR_LM(i) sqrt(gsl_matrix_get(covar_lm,i,i))

    printf("=== LEVMAR DIF SOLUTION (%zu iterations) ===\n", levmar_niter(lm_workspace_p));
    printf ("A      = %.12f +/- %.12f\n", gsl_vector_get(x_lev, 0), ERR_LM(0));
    printf ("lambda = %.12f +/- %.12f\n", gsl_vector_get(x_lev, 1), ERR_LM(1));
    printf ("b      = %.12f +/- %.12f\n", gsl_vector_get(x_lev, 2), ERR_LM(2));

    f.df = &expb_df;
    gsl_vector_memcpy(x_lev, &x.vector);

#if USE_BC
    levmar_bc_proc(x_lev, &f, &lb.vector, &ub.vector, lm_workspace_p);
#else
    levmar_proc(x_lev, &f, lm_workspace_p);
#endif

    printf("=== LEVMAR DER SOLUTION (%zu iterations) ===\n", levmar_niter(lm_workspace_p));
    printf ("A      = %.12f +/- %.12f\n", gsl_vector_get(x_lev, 0), ERR_LM(0));
    printf ("lambda = %.12f +/- %.12f\n", gsl_vector_get(x_lev, 1), ERR_LM(1));
    printf ("b      = %.12f +/- %.12f\n", gsl_vector_get(x_lev, 2), ERR_LM(2));

    fprintf(stderr, "printing covar.levmar\n");
    print_octave(covar_lm, "covar.levmar");

    gsl_vector_free(x_lev);
    levmar_free(lm_workspace_p);
  }

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x,
                                        1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 500);

  printf("=== GSL SOLUTION (%zu iterations) ===\n", iter);

  /*gsl_multifit_covar (s->J, 0.0, covar);*/

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  { 
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

    printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

    printf ("A      = %.12f +/- %.12f\n", FIT(0), c*ERR(0));
    printf ("lambda = %.12f +/- %.12f\n", FIT(1), c*ERR(1));
    printf ("b      = %.12f +/- %.12f\n", FIT(2), c*ERR(2));
  }

  printf ("status = %s\n", gsl_strerror (status));

  fprintf(stderr, "printing covar.gsl\n");
  print_octave(covar, "covar.gsl");

  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  gsl_rng_free (r);
  return 0;
}
