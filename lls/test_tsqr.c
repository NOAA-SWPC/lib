/*
 * test_tsqr.c
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>

#include "oct.h"

static void
random_matrix(gsl_matrix *m, const gsl_rng *r,
              const double lower, const double upper)

{
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double mij = gsl_rng_uniform(r) * (upper - lower) + lower;

          gsl_matrix_set(m, i, j, mij);
        }
    }
}

/* generate random square orthogonal matrix via QR decomposition */
static void
random_matrix_orth(gsl_matrix *m, const gsl_rng *r)
{
  const size_t M = m->size1;
  gsl_matrix *A = gsl_matrix_alloc(M, M);
  gsl_vector *tau = gsl_vector_alloc(M);
  gsl_matrix *R = gsl_matrix_alloc(M, M);

  random_matrix(A, r, -1.0, 1.0);
  gsl_linalg_QR_decomp(A, tau);
  gsl_linalg_QR_unpack(A, tau, m, R);

  gsl_matrix_free(A);
  gsl_matrix_free(R);
  gsl_vector_free(tau);
}

/* construct ill-conditioned matrix via SVD */
static void
random_matrix_ill(gsl_matrix *m, const gsl_rng *r)
{
  const size_t M = m->size1;
  const size_t N = m->size2;
  gsl_matrix *U = gsl_matrix_alloc(M, M);
  gsl_matrix *V = gsl_matrix_alloc(N, N);
  gsl_vector *S = gsl_vector_alloc(N);
  gsl_matrix_view Uv = gsl_matrix_submatrix(U, 0, 0, M, N);
  const double smin = 16.0 * GSL_DBL_EPSILON;
  const double smax = 10.0;
  const double ratio = pow(smin / smax, 1.0 / (N - 1.0));
  double s;
  size_t j;

  random_matrix_orth(U, r);
  random_matrix_orth(V, r);

  /* compute U * S */

  s = smax;
  for (j = 0; j < N; ++j)
    {
      gsl_vector_view uj = gsl_matrix_column(U, j);

      gsl_vector_scale(&uj.vector, s);
      s *= ratio;
    }

  /* compute m = (U * S) * V' */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &Uv.matrix, V, 0.0, m);

  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(S);
}

static void
random_noise(gsl_vector *y, const gsl_rng *r)
{
  size_t i;

  for (i = 0; i < y->size; ++i)
    {
      double *ptr = gsl_vector_ptr(y, i);
      *ptr += 1.0e-3 * gsl_rng_uniform(r);
    }
}

static void
test_solve_tsqr(const double lambda, const gsl_matrix * X, const gsl_vector * y,
                gsl_vector * c)
{
  const size_t nblock = 5;
  const size_t n = X->size1;
  const size_t p = X->size2;
  const size_t nrows = n / nblock; /* number of rows per block */
  tsqr_workspace *w = tsqr_alloc(nrows, p);
  size_t rowidx = 0; /* index of start of current block */
  double rnorm, snorm;

  while (rowidx < n)
    {
      size_t nleft = n - rowidx;
      size_t nr = GSL_MIN(nrows, nleft);
      gsl_matrix_const_view Xv = gsl_matrix_const_submatrix(X, rowidx, 0, nr, p);
      gsl_vector_const_view yv = gsl_vector_const_subvector(y, rowidx, nr);

      tsqr_accumulate(&Xv.matrix, &yv.vector, w);

      rowidx += nr;
    }

  /* compute SVD of R */
  tsqr_svd(w);

  /* solve LS system */
  tsqr_solve(lambda, c, &rnorm, &snorm, w);

  tsqr_free(w);
}

static void
test_solve_gsl(const double lambda, const gsl_matrix * X, const gsl_vector * y,
               gsl_vector * c)
{
  const size_t n = X->size1;
  const size_t p = X->size2;
  gsl_multifit_linear_workspace *w = gsl_multifit_linear_alloc(n, p);
  double rnorm, snorm;

  gsl_multifit_linear_svd(X, w);
  gsl_multifit_linear_solve(lambda, X, y, c, &rnorm, &snorm, w);

  gsl_multifit_linear_free(w);
}

static int
test_system(const size_t n, const size_t p, gsl_rng *r)
{
  int s = 0;
  const double tol = 1.0e-5;
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_vector *c = gsl_vector_alloc(p);
  gsl_vector *c0 = gsl_vector_alloc(p);
  gsl_vector *c1 = gsl_vector_alloc(p);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *w = gsl_vector_alloc(n);
  gsl_vector *L = gsl_vector_alloc(p);
  double lambda;
  size_t i, j;

  /* generate ill-conditioned random matrix */
  random_matrix_ill(X, r);

  /* generate random solution vector */
  random_vector(c, r, -1.0, 1.0);

  /* compute y = X*c and add random noise */
  gsl_blas_dgemv(CblasNoTrans, 1.0, X, c, 0.0, y);
  random_noise(y, r);

  for (i = 0; i < 5; ++i)
    {
      lambda = pow(10.0, -(double) i);

      test_solve_tsqr(lambda, X, y, c0);
      test_solve_gsl(lambda, X, y, c1);

      /* test c0 = c1 */
      for (j = 0; j < p; ++j)
        {
          double c0j = gsl_vector_get(c0, j);
          double c1j = gsl_vector_get(c1, j);

          gsl_test_rel(c0j, c1j, tol, "tsqr lambda=%g n=%zu p=%zu j=%zu", lambda, n, p, j);
        }
    }

  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_vector_free(c0);
  gsl_vector_free(c1);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(L);

  return s;
}

void
test_tsqr(gsl_rng *r)
{
  test_system(163, 87, r);
  test_system(200, 70, r);
  test_system(503, 452, r);
}
