/*
 * levmar.h
 */

#ifndef INCLUDED_levmar_h
#define INCLUDED_levmar_h

#include <levmar_lib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

typedef struct
{
  size_t n;      /* number of data */
  size_t p;      /* number of parameters */
  double opts[LM_OPTS_SZ];
  double info[LM_INFO_SZ];
  int maxit;     /* maximum iterations */

  gsl_vector *rhs; /* rhs vector */

  gsl_matrix *covar; /* covariance matrix */
} levmar_workspace;

typedef struct
{
  gsl_multifit_function_fdf *fdf;
} levmar_params;

/*
 * Prototypes
 */

levmar_workspace *levmar_alloc(const size_t n, const size_t p);
void levmar_free(levmar_workspace *w);
int levmar_proc(gsl_vector *x, gsl_multifit_function_fdf *fdf,
                levmar_workspace *w);
int levmar_bc_proc(gsl_vector *x, gsl_multifit_function_fdf *fdf,
                   gsl_vector *lb, gsl_vector *ub, levmar_workspace *w);
double levmar_chisq(const levmar_workspace *w);
double levmar_chisq0(const levmar_workspace *w);
size_t levmar_niter(const levmar_workspace *w);

#endif /* INCLUDED_levmar_h */
