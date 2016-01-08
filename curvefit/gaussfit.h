/*
 * gaussfit.h
 */

#ifndef INCLUDED_gaussfit_h
#define INCLUDED_gaussfit_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit_nlinear.h>

typedef struct
{
  size_t n;       /* number of observations */
  size_t p;       /* number of model coefficients */

  double chisq0;  /* initial chi^2 */
  double chisq;   /* final chi^2 */

  /*
   * coefficient vector is organized as:
   * c = [ a0 a1 a2 a3 a4 a5 ]
   */
  gsl_vector *c;  /* model coefficients */

  gsl_multifit_nlinear_workspace *nlinear_workspace_p;
} gaussfit_workspace;

typedef struct
{
  double *t;
  double *y;
  gaussfit_workspace *w;
} gaussfit_data;

/*
 * Prototypes
 */

gaussfit_workspace *gaussfit_alloc(const size_t n, const size_t p);
void gaussfit_free(gaussfit_workspace *w);
int gaussfit_init(const gsl_vector *x, gaussfit_workspace *w);
int gaussfit(const double *t, const double *y, gaussfit_workspace *w);
double gaussfit_eval(const gsl_vector *x, const double t);

#endif /* INCLUDED_gaussfit_h */
