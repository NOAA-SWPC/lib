#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "levmar.h"
#include "oct.h"

#include "test_powell.c"

void print_state (size_t iter, gsl_multifit_fdfsolver * s);

int
main (void)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int i, iter = 0;
  const size_t n = powell_N;
  const size_t p = powell_P;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  gsl_multifit_function_fdf f;
  gsl_vector_view x = gsl_vector_view_array (powell_x0, p);
  gsl_rng * r;

  f.f = &powell_f;
  f.df = &powell_df;
  f.fdf = &powell_fdf;
  f.n = n;
  f.p = p;
  f.params = NULL;

  {
    levmar_workspace *lm_workspace_p = levmar_alloc(n, p);
    gsl_vector *x_lev = gsl_vector_alloc(p);
    gsl_matrix *covar_lm = lm_workspace_p->covar;

    f.df = NULL;
    f.fdf = NULL;
    gsl_vector_memcpy(x_lev, &x.vector);
    levmar_proc(x_lev, &f, lm_workspace_p);

#define ERR_LM(i) sqrt(gsl_matrix_get(covar_lm,i,i))

    printf("=== LEVMAR DIF SOLUTION ===\n");
    printf ("x0 = %.12f +/- %.12f\n", gsl_vector_get(x_lev, 0), ERR_LM(0));
    printf ("x1 = %.12f +/- %.12f\n", gsl_vector_get(x_lev, 1), ERR_LM(1));

    f.df = &powell_df;
    f.fdf = &powell_fdf;
    gsl_vector_memcpy(x_lev, &x.vector);
    levmar_proc(x_lev, &f, lm_workspace_p);

    printf("=== LEVMAR DER SOLUTION ===\n");
    printf ("x0 = %.12f +/- %.12f\n", gsl_vector_get(x_lev, 0), ERR_LM(0));
    printf ("x1 = %.12f +/- %.12f\n", gsl_vector_get(x_lev, 1), ERR_LM(1));

    fprintf(stderr, "printing covar.levmar\n");
    print_octave(covar_lm, "covar.levmar");

    gsl_vector_free(x_lev);
    levmar_free(lm_workspace_p);
  }

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  printf("=== GSL SOLUTION ===\n");

#if 0
  print_state (iter, s);
#endif

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

#if 0
      printf ("status = %s\n", gsl_strerror (status));
      print_state (iter, s);
#endif

      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x,
                                        1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 500);

  /*gsl_multifit_covar (s->J, 0.0, covar);*/

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  { 
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

    printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

    printf ("x0 = %.12f +/- %.12f\n", FIT(0), c*ERR(0));
    printf ("x1 = %.12f +/- %.12f\n", FIT(1), c*ERR(1));
  }

  printf ("status = %s\n", gsl_strerror (status));

  fprintf(stderr, "printing covar.gsl\n");
  print_octave(covar, "covar.gsl");

  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  return 0;
}

void
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3zu % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_blas_dnrm2 (s->f));
}
