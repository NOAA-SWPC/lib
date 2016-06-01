/*
 * test.c
 * Patrick Alken
 *
 * Test superlu module by creating random sparse matrices and rhs
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

#include "superlu.h"

void
test_vectors(double *x, double *x_exact, int n)
{
  int i;

  for (i = 0; i < n; ++i)
    gsl_test_rel(x[i], x_exact[i], 1.0e-7, "n = %d, i = %d", n, i);
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

  gsl_matrix_set_zero(m);

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
      for (j = 0; j < (size_t) numel; ++j)
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
test_superlu()
{
  const int N_max = 100;
  int n, i;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  for (n = 1; n <= N_max; ++n)
    {
      gsl_spmatrix *S = gsl_spmatrix_alloc(n, n);
      gsl_spmatrix *C;
      gsl_matrix *A = gsl_matrix_alloc(n, n);
      gsl_matrix *A_copy = gsl_matrix_alloc(n, n);
      gsl_vector *rhs = gsl_vector_alloc(n);
      gsl_permutation *p = gsl_permutation_alloc(n);
      gsl_vector *x_gsl = gsl_vector_alloc(n);
      gsl_vector *x_slu = gsl_vector_alloc(n);
      slu_workspace *w = slu_alloc(n, n, 1);
      double scale = gsl_rng_uniform(r);
      int s;

      for (i = 0; i < 50; ++i)
        {
          create_random_sparse_matrix(A, r, -10.0, 10.0);
          create_random_vector(rhs, r, -10.0, 10.0);

          gsl_matrix_memcpy(A_copy, A);

          gsl_linalg_LU_decomp(A, p, &s);
          gsl_linalg_LU_solve(A, p, rhs, x_gsl);

          gsl_matrix_memcpy(A, A_copy);

          /* convert dense matrix to sparse (triplet) format */
          gsl_spmatrix_d2sp(S, A);

          /* convert to compressed column format */
          C = gsl_spmatrix_compcol(S);

          gsl_spmatrix_scale(C, scale);
          slu_proc(C, rhs->data, x_slu->data, w);
          gsl_vector_scale(x_slu, scale);

          gsl_test_abs(w->rnorm, 0.0, 1.0e-8, "residual n = %zu, i = %zu cond = %.12e",
                       n, i, 1.0 / w->rcond);

          test_vectors(x_slu->data, x_gsl->data, n);

          gsl_spmatrix_free(C);
        }

      gsl_spmatrix_free(S);
      gsl_matrix_free(A);
      gsl_matrix_free(A_copy);
      gsl_vector_free(rhs);
      gsl_vector_free(x_gsl);
      gsl_vector_free(x_slu);
      gsl_permutation_free(p);
      slu_free(w);
    }

  gsl_rng_free(r);
} /* test_superlu() */

int
main()
{
  slu_workspace *sw;
  gsl_spmatrix *A, *C;
  size_t i;
  const int m = 5;
  const int n = 5;
  double rhs[m];
  double sol[n];
  int nprocs;
  double s, u, p, e, r, l;
  double min, max;
  double x[] = { -0.0312500000000000, 0.0654761904761905, 0.0133928571428571,
                  0.0625000000000000, 0.0327380952380952 };

  nprocs = 1;
  s = 19.0;
  u = 21.0;
  p = 16.0;
  e = 5.0;
  r = 18.0;
  l = 12.0;

  for (i = 0; i < (size_t) n; ++i)
    rhs[i] = 1.0;

  sw = slu_alloc(m, n, nprocs);
  A = gsl_spmatrix_alloc(m, n);

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

  gsl_spmatrix_minmax(A, &min, &max);
  fprintf(stderr, "min = %f, max = %f\n", min, max);

  C = gsl_spmatrix_compcol(A);

  slu_proc(C, rhs, sol, sw);

  test_vectors(sol, x, n);

  printf("sol = [\n");
  for (i = 0; i < (size_t) n; ++i)
    printf("%.12e\n", sol[i]);
  printf("]\n");

  printf("residual = %.12e\n", slu_residual(sw));
  printf("cond(A)  = %.12e\n", 1.0 / sw->rcond);

  test_superlu();

  slu_free(sw);
  gsl_spmatrix_free(A);
  gsl_spmatrix_free(C);

  exit (gsl_test_summary());

  return 0;
}
