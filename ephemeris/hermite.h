/*
 * hermite.h
 */

#ifndef INCLUDED_hermite_h
#define INCLUDED_hermite_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>

typedef struct
{
  size_t n;      /* total number of points to interpolate */
  double *q;     /* Hermite polynomial coefficients */
  double *z;     /* Hermite polynomial z values */
  double *coeff; /* Taylor coefficients */
  double *work;  /* additional workspace */
  size_t degree; /* degree of polynomial */
  size_t ncoeff; /* number of polynomial coefficients */
  size_t npts;   /* number of points needed for piecewise interpolation */
  size_t max_degree; /* maximum possible polynomial degree */
} hermite_workspace;

hermite_workspace *hermite_alloc(const size_t n, size_t degree);
void hermite_free(hermite_workspace *w);
int hermite_init(const double xa[], const double ya[],
                 const double dya[], const size_t size,
                 hermite_workspace *w);
double hermite_eval(const double xa[], const double ya[], 
                    const double dya[], double x,
                    gsl_interp_accel *acc,
                    hermite_workspace *w);
double hermite_eval_deriv(const double xa[], const double ya[], 
                          const double dya[], double x,
                          gsl_interp_accel *acc,
                          hermite_workspace *w);
double hermite_eval_deriv2(const double xa[], const double ya[], 
                           const double dya[], double x,
                           gsl_interp_accel *acc,
                           hermite_workspace *w);

#endif /* INCLUDED_hermite_h */
