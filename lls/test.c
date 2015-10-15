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

#include "test_shaw.c"

static void
random_complex_matrix(gsl_matrix_complex *m, gsl_rng *r,
                      const double lower, const double upper)

{
  const size_t N = m->size1;
  const size_t M = m->size2;
  size_t i, j;

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < M; ++j)
        {
          double aij = gsl_rng_uniform(r) * (upper - lower) + lower;
          double bij = gsl_rng_uniform(r) * (upper - lower) + lower;
          gsl_complex mij = gsl_complex_rect(aij, bij);

          gsl_matrix_complex_set(m, i, j, mij);
        }
    }
} /* random_complex_matrix() */

static void
random_complex_vector(gsl_vector_complex *v, gsl_rng *r,
                      const double lower, const double upper)

{
  const size_t N = v->size;
  size_t i;

  for (i = 0; i < N; ++i)
    {
      double ai = gsl_rng_uniform(r) * (upper - lower) + lower;
      double bi = gsl_rng_uniform(r) * (upper - lower) + lower;
      gsl_complex vi = gsl_complex_rect(ai, bi);

      gsl_vector_complex_set(v, i, vi);
    }
} /* random_complex_vector() */

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

/*
test_complex()
  Test complex LS system using a random LS matrix and RHS vector.
The system is first partioned into two parts, folded into the LLS
system and solved. Then, the system is solved whole without partitioning,
and the two solutions are compared.
*/

int
test_complex(const size_t n, const size_t p, gsl_rng *r)
{
  int s = 0;
  const double tol = 1.0e-10;
  const size_t nblock = (size_t) (0.75 * n);
  gsl_matrix_complex *A = gsl_matrix_complex_alloc(n, p);
  gsl_vector_complex *b = gsl_vector_complex_alloc(n);
  gsl_vector_complex *c1 = gsl_vector_complex_alloc(p);
  gsl_vector_complex *c2 = gsl_vector_complex_alloc(p);
  gsl_matrix_complex *work_A = gsl_matrix_complex_alloc(n, p);
  gsl_vector_complex *work_b = gsl_vector_complex_alloc(n);
  gsl_vector *w = gsl_vector_alloc(n);
  lls_complex_workspace *lls1 = lls_complex_alloc(nblock, p);
  lls_complex_workspace *lls2 = lls_complex_alloc(n, p);
  gsl_matrix_complex_view Av;
  gsl_vector_complex_view bv;
  gsl_vector_view wv;
  size_t i;

  random_complex_matrix(A, r, -1.0, 1.0);
  random_complex_vector(b, r, -1.0, 1.0);
  random_vector(w, r, 0.0, 1.0);

  gsl_vector_complex_memcpy(work_b, b);
  gsl_matrix_complex_memcpy(work_A, A);

  /* first part */
  Av = gsl_matrix_complex_submatrix(work_A, 0, 0, nblock, p);
  bv = gsl_vector_complex_subvector(work_b, 0, nblock);
  wv = gsl_vector_subvector(w, 0, nblock);
  lls_complex_fold(&Av.matrix, &bv.vector, &wv.vector, lls1);

  /* second part */
  Av = gsl_matrix_complex_submatrix(work_A, nblock, 0, n - nblock, p);
  bv = gsl_vector_complex_subvector(work_b, nblock, n - nblock);
  wv = gsl_vector_subvector(w, nblock, n - nblock);
  lls_complex_fold(&Av.matrix, &bv.vector, &wv.vector, lls1);

  /* solve partitioned system */
  lls_complex_solve(c1, lls1);

  /* build and solve complete system */
  gsl_vector_complex_memcpy(work_b, b);
  gsl_matrix_complex_memcpy(work_A, A);
  lls_complex_fold(work_A, work_b, w, lls2);
  lls_complex_solve(c2, lls2);

  /* check that chi^2 matches */
  gsl_test_rel(lls1->chisq, lls2->chisq, tol, "chisq");

  /* check coefficient vectors */
  for (i = 0; i < p; ++i)
    {
      gsl_complex c1i = gsl_vector_complex_get(c1, i);
      gsl_complex c2i = gsl_vector_complex_get(c2, i);

      gsl_test_rel(GSL_REAL(c1i), GSL_REAL(c2i), tol, "real i=%zu", i);
      gsl_test_rel(GSL_IMAG(c1i), GSL_IMAG(c2i), tol, "imag i=%zu", i);
    }

  gsl_matrix_complex_free(A);
  gsl_vector_complex_free(b);
  gsl_vector_complex_free(c1);
  gsl_vector_complex_free(c2);
  gsl_vector_complex_free(work_b);
  gsl_matrix_complex_free(work_A);
  gsl_vector_free(w);
  lls_complex_free(lls1);
  lls_complex_free(lls2);

  return s;
} /* test_complex() */

int
main(int argc, char *argv[])
{
  int s = 0;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  test_shaw(r);

  s += test_complex(100, 50, r);

  gsl_rng_free(r);

  exit (gsl_test_summary());
}
