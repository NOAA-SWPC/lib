/*
 * test_complex.c
 */

#include <lapacke/lapacke.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "oct.h"

static void
random_vector_complex(gsl_vector_complex *v, gsl_rng *r,
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
} /* random_vector_complex() */

static void
random_vector_complex_noise(gsl_vector_complex *v, gsl_rng *r)

{
  const size_t N = v->size;
  size_t i;

  for (i = 0; i < N; ++i)
    {
      gsl_complex vi = gsl_vector_complex_get(v, i);
      double sigma1 = 0.1 * GSL_REAL(vi);
      double sigma2 = 0.1 * GSL_IMAG(vi);
      double ai = gsl_ran_gaussian(r, sigma1);
      double bi = gsl_ran_gaussian(r, sigma2);
      gsl_complex val = gsl_complex_rect(GSL_REAL(vi) + ai, GSL_IMAG(vi) + bi);

      gsl_vector_complex_set(v, i, val);
    }
} /* random_vector_complex_noise() */

static void
random_matrix_complex(gsl_matrix_complex *m, gsl_rng *r,
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
} /* random_matrix_complex() */

static void
random_unitary(gsl_matrix_complex *m, gsl_rng *r)
{
  const size_t M = m->size1;
  const size_t N = m->size2;
  lapack_int K = GSL_MIN(M, N);
  lapack_int lda = M;
  lapack_int ldc = M;
  gsl_vector_complex *tau = gsl_vector_complex_alloc(GSL_MIN(M, N));
  gsl_matrix_complex *A = gsl_matrix_complex_alloc(M, N);
  gsl_matrix_complex *m2 = gsl_matrix_complex_alloc(N, M);

  random_matrix_complex(A, r, -1.0, 1.0);
  gsl_matrix_complex_set_identity(m);

  gsl_matrix_complex_transpose_memcpy(m2, m);

  LAPACKE_zgeqrf(LAPACK_COL_MAJOR,
                 M,
                 N,
                 (lapack_complex_double *) A->data,
                 lda,
                 (lapack_complex_double *) tau->data);

  LAPACKE_zunmqr(LAPACK_COL_MAJOR,
                 'L',
                 'N',
                 M,
                 N,
                 K,
                 (lapack_complex_double *) A->data,
                 lda,
                 (lapack_complex_double *) tau->data,
                 (lapack_complex_double *) m2->data,
                 ldc);

  gsl_matrix_complex_transpose_memcpy(m, m2);

  gsl_vector_complex_free(tau);
  gsl_matrix_complex_free(A);
  gsl_matrix_complex_free(m2);
} /* random_unitary() */

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

void
random_matrix_complex_ill(gsl_matrix_complex * m, gsl_rng *r)
{
  const size_t M = m->size1;
  const size_t N = m->size2;
  gsl_matrix_complex *U = gsl_matrix_complex_alloc(M, N);
  gsl_matrix_complex *V = gsl_matrix_complex_alloc(N, N);
  size_t j;
  const double smin = 16.0 * GSL_DBL_EPSILON;
  const double smax = 10.0;
  const double ratio = pow(smin / smax, 1.0 / (N - 1.0));
  double s;

  random_unitary(U, r);
  random_unitary(V, r);
  
  /* calculate U * S */

  s = smax;
  for (j = 0; j < N; ++j)
    {
      gsl_vector_complex_view uj = gsl_matrix_complex_column(U, j);
      gsl_complex val = gsl_complex_rect(s, 0.0);

      gsl_vector_complex_scale(&uj.vector, val);

      s *= ratio;
    }

  /* compute (U * S) * V' */
  gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, U, V, GSL_COMPLEX_ZERO, m);

  gsl_matrix_complex_free(U);
  gsl_matrix_complex_free(V);
}

void
test_complex_solve(const double lambda, const gsl_matrix_complex * X, const gsl_vector_complex * y,
                   const gsl_vector * w, const gsl_vector * L, const char * desc)
{
  const size_t n = X->size1;
  const size_t p = X->size2;
  size_t i;
  const double tol = 1.0e-7;

  gsl_vector_complex *c1 = gsl_vector_complex_alloc(p);
  gsl_vector_complex *c2 = gsl_vector_complex_alloc(p);

  gsl_matrix_complex *Xtmp = gsl_matrix_complex_alloc(n, p);
  gsl_vector_complex *ytmp = gsl_vector_complex_alloc(n);

  /* solve with LLS module */
  {
    lls_complex_workspace *lls_p = lls_complex_alloc(10, p);

    gsl_matrix_complex_memcpy(Xtmp, X);
    gsl_vector_complex_memcpy(ytmp, y);

    lls_complex_stdform(Xtmp, ytmp, w, L, lls_p);
    lls_complex_fold(Xtmp, ytmp, lls_p);
    lls_complex_solve(lambda, c1, lls_p);
    lls_complex_btransform(L, c1, lls_p);

    lls_complex_free(lls_p);
  }

  /* solve by forming (X^T W X + \lambda^2 I)^{-1} X^T W y */
  {
    gsl_matrix_complex *XHX = gsl_matrix_complex_alloc(p, p); /* X^H W X */
    gsl_vector_complex *XHy = gsl_vector_complex_alloc(p);    /* X^H W y */

    gsl_matrix_complex_memcpy(Xtmp, X);
    gsl_vector_complex_memcpy(ytmp, y);

    /* compute Xtmp = W X and ytmp = W y */
    for (i = 0; i < n; ++i)
      {
        gsl_vector_complex_view Xi = gsl_matrix_complex_row(Xtmp, i);
        gsl_complex yi = gsl_vector_complex_get(ytmp, i);
        double wi = gsl_vector_get(w, i);
        gsl_complex val = gsl_complex_rect(wi, 0.0);

        gsl_vector_complex_scale(&Xi.vector, val);
        gsl_vector_complex_set(ytmp, i, gsl_complex_mul(yi, val));
      }

    /* compute XHy = X^H W y */
    gsl_blas_zgemv(CblasConjTrans, GSL_COMPLEX_ONE, X, ytmp, GSL_COMPLEX_ZERO, XHy);

    /* compute XHX = X^H W X */
    gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, GSL_COMPLEX_ONE, X, Xtmp, GSL_COMPLEX_ZERO, XHX);

    /* add lambda^2 L^H L to diagonal */
    for (i = 0; i < p; ++i)
      {
        gsl_complex Xii = gsl_matrix_complex_get(XHX, i, i);
        double Li = gsl_vector_get(L, i);
        double val = pow(lambda * Li, 2.0);

        gsl_matrix_complex_set(XHX, i, i, gsl_complex_add_real(Xii, val));
      }

    /* solve XHX c2 = ytmp */
    {
      gsl_matrix_complex_transpose(XHX);
      gsl_vector_complex_memcpy(c2, XHy);
      LAPACKE_zposv(LAPACK_COL_MAJOR,
                    'U',
                    p,
                    1,
                    (lapack_complex_double *) XHX->data,
                    p,
                    (lapack_complex_double *) c2->data,
                    p);
    }

    gsl_matrix_complex_free(XHX);
    gsl_vector_complex_free(XHy);
  }

  /* now test if c1 = c2 */
  for (i = 0; i < p; ++i)
    {
      gsl_complex c1i = gsl_vector_complex_get(c1, i);
      gsl_complex c2i = gsl_vector_complex_get(c2, i);

      gsl_test_rel(GSL_REAL(c1i), GSL_REAL(c2i), tol, "%s, real lambda=%g n=%zu p=%zu i=%zu", desc, lambda, n, p, i);
      gsl_test_rel(GSL_IMAG(c1i), GSL_IMAG(c2i), tol, "%s, imag lambda=%g n=%zu p=%zu i=%zu", desc, lambda, n, p, i);
    }

  printcv_octave(c1, "c1");
  printcv_octave(c2, "c2");

  gsl_vector_complex_free(c1);
  gsl_vector_complex_free(c2);
  gsl_matrix_complex_free(Xtmp);
  gsl_vector_complex_free(ytmp);
}

int
test_complex_system(const size_t n, const size_t p, gsl_rng *r)
{
  int s = 0;
  gsl_matrix_complex *X = gsl_matrix_complex_alloc(n, p);
  gsl_vector_complex *c = gsl_vector_complex_alloc(p);
  gsl_vector_complex *y = gsl_vector_complex_alloc(n);
  gsl_vector *w = gsl_vector_alloc(n);
  gsl_vector *L = gsl_vector_alloc(p);

  /* generate ill-conditioned random matrix */
  random_matrix_complex_ill(X, r);

  /* generate random solution vector */
  random_vector_complex(c, r, -1.0, 1.0);

  /* compute y = X*c and add random noise */
  gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, X, c, GSL_COMPLEX_ZERO, y);
  random_vector_complex_noise(y, r);

  /* test with unit weights */
  gsl_vector_set_all(w, 1.0);

  /* test with L = I */
  gsl_vector_set_all(L, 1.0);

  test_complex_solve(1.0, X, y, w, L, "unit w, unit L");
  test_complex_solve(0.1, X, y, w, L, "unit w, unit L");
  test_complex_solve(0.01, X, y, w, L, "unit w, unit L");

  /* test with random L */
  random_vector(L, r, -5.0, 5.0);

  test_complex_solve(1.0, X, y, w, L, "unit w, random L");
  test_complex_solve(0.1, X, y, w, L, "unit w, random L");
  test_complex_solve(0.01, X, y, w, L, "unit w, random L");

  /* test with random weights */
  random_vector(w, r, 0.0, 1.0);

  /* test with L = I */
  gsl_vector_set_all(L, 1.0);

  test_complex_solve(1.0, X, y, w, L, "random w, unit L");
  test_complex_solve(0.1, X, y, w, L, "random w, unit L");
  test_complex_solve(0.01, X, y, w, L, "random w, unit L");

  /* test with random L */
  random_vector(L, r, -5.0, 5.0);

  test_complex_solve(1.0, X, y, w, L, "random w, random L");
  test_complex_solve(0.1, X, y, w, L, "random w, random L");
  test_complex_solve(0.01, X, y, w, L, "random w, random L");

  printc_octave(X, "X");
  printcv_octave(c, "c");
  printcv_octave(y, "y");
  printv_octave(w, "w");
  printv_octave(L, "L");

  gsl_matrix_complex_free(X);
  gsl_vector_complex_free(c);
  gsl_vector_complex_free(y);
  gsl_vector_free(w);
  gsl_vector_free(L);

  return s;
}

void
test_complex(gsl_rng *r)
{
  test_complex_system(100, 50, r);
}
