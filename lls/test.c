/*
 * test.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_test.h>

#include "lls.h"
#include "tsqr.h"

static void random_vector(gsl_vector *v, gsl_rng *r,
                          const double lower, const double upper);

#include "test_complex.c"
#include "test_shaw.c"
#include "test_tsqr.c"

static void
random_vector(gsl_vector *v, gsl_rng *r,
              const double lower, const double upper)

{
  const size_t N = v->size;
  size_t i;

  for (i = 0; i < N; ++i)
    {
      double vi = gsl_rng_uniform(r) * (upper - lower) + lower;

      gsl_vector_set(v, i, vi);
    }
} /* random_vector() */

int
main(int argc, char *argv[])
{
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

#if 0
  test_shaw(r);
  test_complex(r);
#endif
  test_tsqr(r);

  gsl_rng_free(r);

  exit (gsl_test_summary());
}
