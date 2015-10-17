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

#include "test_complex.c"
#include "test_shaw.c"

int
main(int argc, char *argv[])
{
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

#if 0
  test_shaw(r);
#endif
  test_complex(r);

  gsl_rng_free(r);

  exit (gsl_test_summary());
}
