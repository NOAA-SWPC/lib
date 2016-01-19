/*
 * test.c
 *
 * Test lis module by creating random sparse matrices and rhs
 * vectors, solving the systems and comparing with GSL output
 */

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_spmatrix.h>

#include "lisw.h"

void
test_vectors(double *x, double *x_exact, int n)
{
  int i;

  for (i = 0; i < n; ++i)
    gsl_test_rel(x[i], x_exact[i], 1.0e-9, "n = %d, i = %d", n, i);
}

void
create_random_sparse_matrix(gsl_matrix *m, gsl_rng *r, double lower,
                            double upper)
{
  size_t N = m->size1;
  size_t M = m->size2;
  size_t i, j;
  int numel; /* number of non-zero elements in a row */
  double x;

  if (N == 1)
    {
      x = gsl_rng_uniform(r) * (upper - lower) + lower;
      gsl_matrix_set(m, 0, 0, x);
      return;
    }

  for (i = 0; i < N; ++i)
    {
      /* pick a random number between 1 and M/2 - this is how many
       * nonzero elements are in this row
       */
      numel = (int) (gsl_rng_uniform(r) * (M / 2 - 1) + 1);
      for (j = 0; j < numel; ++j)
        {
          int k = (int) (gsl_rng_uniform(r) * (M - 2));
          x = gsl_rng_uniform(r) * (upper - lower) + lower;
          gsl_matrix_set(m, i, k, x);
        }

      /* always set the diagonal element */
      x = gsl_rng_uniform(r) * (upper - lower) + lower;
      gsl_matrix_set(m, i, i, x);
    }
} /* create_random_sparse_matrix() */

void
create_random_vector(gsl_vector *v, gsl_rng *r, double lower, double upper)
{
  size_t i;

  for (i = 0; i < v->size; ++i)
    {
      double x = gsl_rng_uniform(r) * (upper - lower) + lower;
      gsl_vector_set(v, i, x);
    }
}

void
test_lis()
{
#if 0
  const int N_max = 50;
  int n, i;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  for (n = 1; n <= N_max; ++n)
    {
      gsl_matrix *A = gsl_matrix_alloc(n, n);
      gsl_matrix *A_copy = gsl_matrix_alloc(n, n);
      gsl_vector *rhs = gsl_vector_alloc(n);
      gsl_permutation *p = gsl_permutation_alloc(n);
      gsl_vector *x_gsl = gsl_vector_alloc(n);
      gsl_vector *x_lis = gsl_vector_alloc(n);
      lis_workspace *w = lis_alloc(n, n);
      double scale = gsl_rng_uniform(r);
      int s;

      gsl_matrix_set_zero(A);

      for (i = 0; i < 1; ++i)
        {
          create_random_sparse_matrix(A, r, -10.0, 10.0);
          create_random_vector(rhs, r, -10.0, 10.0);

          gsl_matrix_memcpy(A_copy, A);

          gsl_linalg_LU_decomp(A, p, &s);
          gsl_linalg_LU_solve(A, p, rhs, x_gsl);

          gsl_matrix_memcpy(A, A_copy);

          lis_reinit(w);
          test_cs_add(A, w);
          /*lis_scale(scale, w);*/
          lis_proc(rhs->data, x_lis->data, w);
          /*gsl_vector_scale(x_lis, scale);*/

          test_vectors(x_lis->data, x_gsl->data, n);
        }

      gsl_matrix_free(A);
      gsl_matrix_free(A_copy);
      gsl_vector_free(rhs);
      gsl_vector_free(x_gsl);
      gsl_vector_free(x_lis);
      gsl_permutation_free(p);
      lis_free(w);
    }

  gsl_rng_free(r);
#endif
} /* test_lis() */

int
main(int argc, char *argv[])
{
  lis_workspace *w;
  size_t i;
  const int m = 5;
  const int n = 5;
  gsl_spmatrix *A = gsl_spmatrix_alloc(m, n);
  double rhs[m];
  double sol[n];
  int nprocs;
  double s, u, p, e, r, l;
  double x[] = { -0.0312500000000000, 0.0654761904761905, 0.0133928571428571,
                  0.0625000000000000, 0.0327380952380952 };

  nprocs = 1;
  s = 19.0;
  u = 21.0;
  p = 16.0;
  e = 5.0;
  r = 18.0;
  l = 12.0;

  for (i = 0; i < n; ++i)
    rhs[i] = 1.0;

  w = lis_alloc(m, n);

  gsl_spmatrix_set(A, 0, 0, s);
  gsl_spmatrix_set(A, 1, 0, l);
  gsl_spmatrix_set(A, 4, 0, l);
  gsl_spmatrix_set(A, 1, 1, u);
  gsl_spmatrix_set(A, 2, 1, l);
  gsl_spmatrix_set(A, 4, 1, l);
  gsl_spmatrix_set(A, 0, 2, u);
  gsl_spmatrix_set(A, 2, 2, p);
  gsl_spmatrix_set(A, 0, 3, u);
  gsl_spmatrix_set(A, 3, 3, e);
  gsl_spmatrix_set(A, 3, 4, u);
  gsl_spmatrix_set(A, 4, 4, r);

  lis_proc(A, rhs, 1.0e-12, sol, w);

  test_vectors(sol, x, n);

  printf("sol = [\n");
  for (i = 0; i < n; ++i)
    printf("%.12e\n", sol[i]);
  printf("]\n");

  printf("residual = %.12e\n", w->residual);

  test_lis();

  mylis_free(w);

  gsl_spmatrix_free(A);

  exit(gsl_test_summary());

  return 0;
}
