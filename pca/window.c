#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_bessel.h>

int
apply_ps1(gsl_vector *v)
{
  const double p[11] = { 5.3476939016920851e-11, 2.2654256220146656e-9,
                         7.8075102004229667e-8, 2.1373409644281953e-6,
                         4.5094847544714943e-5, 7.0498957221483167e-4,
                         7.7412693304064753e-3, 5.5280627452077586e-2,
                         2.2753754228751827e-1, 4.3433904277546202e-1,
                         2.2902051859068017e-1 };
  size_t n = v->size;
  double d = sqrt(2.0 / (double)n);
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double xt = (2.0*i-1.0)/(double)n - 1.0;
      double u = (1.0 - xt) * (1.0 + xt);
      double wi = ((((((((((p[0]*u+p[1])*u+p[2])*u+p[3])*u+p[4])*u+p[5])*u+p[6])*u+p[7])*u+p[8])*u+p[9])*u+p[10])*d;
      double vi;

      vi = gsl_vector_get(v, i);
      gsl_vector_set(v, i, vi * wi);
    }

  return 0;
}

int
apply_hamming(gsl_vector *v)
{
  const size_t N = v->size;
  const double alpha = 0.53836;
  const double beta = 0.46164;
  size_t i;

  for (i = 0; i < N; ++i)
    {
      double ratio = 2.0 * M_PI * i / ((double)N - 1.0);
      double wi = alpha - beta*cos(ratio);
      double *ptr = gsl_vector_ptr(v, i);

      *ptr *= wi;
    }

  return 0;
}

int
apply_kaiser(gsl_vector *v)
{
  const size_t n = v->size;
  const double alpha = 3.0;
  const double beta = M_PI * alpha;
  const double d = gsl_sf_bessel_I0(beta);
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double term1 = 2.0*i / (n - 1.0) - 1.0;
      double term2 = 1.0 - term1*term1;
      double wi = gsl_sf_bessel_I0(beta * sqrt(term2)) / d;
      double *ptr = gsl_vector_ptr(v, i);

      *ptr *= wi;
    }

  return 0;
}
