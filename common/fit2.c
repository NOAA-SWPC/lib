#include <stdio.h>
#include <math.h>

#include "fit2.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

/*
fit2_proc()
  Perform the fit by solving the system

A x = b

with

x = (A^t A)^{-1} A^t b

Inputs: A - design matrix
        b - rhs side
        x - where to store unknown coefficients

Return: 1 on success, 0 on failure
*/

int
fit2_proc(gsl_matrix *A, gsl_vector *b, gsl_vector *x)

{
  int s;
  gsl_matrix *Af;
  gsl_vector *bf;
  gsl_permutation *p;

  Af = gsl_matrix_alloc(A->size2, A->size2);
  bf = gsl_vector_alloc(A->size2);
  p = gsl_permutation_alloc(A->size2);

  /* compute Af = A^t A */
  gsl_blas_dgemm(CblasTrans,
                 CblasNoTrans,
                 1.0,
                 A,
                 A,
                 0.0,
                 Af);

  /* compute bf = A^t b */
  gsl_blas_dgemv(CblasTrans,
                 1.0,
                 A,
                 b,
                 0.0,
                 bf);

  /* we now wish to solve the equation Af x = bf */

  gsl_linalg_LU_decomp(Af, p, &s);
  gsl_linalg_LU_solve(Af, p, bf, x);

  gsl_matrix_free(Af);
  gsl_vector_free(bf);
  gsl_permutation_free(p);

  return (1);
} /* fit2_proc() */

/*
Solve X c = y using weighted data in w

so that

X' W X c = X' W y

and

c = (X' W X)^{-1} X' W y
*/

int
fit2_wproc(const gsl_matrix *X, const gsl_vector *y, const gsl_vector *w,
           gsl_vector *c)
{
  size_t i;
  int s = 0;
  gsl_matrix *Xf, *Xbar;
  gsl_vector *yf, *ybar;

  Xbar = gsl_matrix_alloc(X->size1, X->size2);
  ybar = gsl_vector_alloc(y->size);

  Xf = gsl_matrix_alloc(X->size2, X->size2);
  yf = gsl_vector_alloc(X->size2);

  gsl_matrix_memcpy(Xbar, X);

  /* Xbar = sqrt(W) X, ybar = sqrt(W) y */
  for (i = 0; i < X->size1; ++i)
    {
      double wi = gsl_vector_get(w, i);
      gsl_vector_view v = gsl_matrix_row(Xbar, i);

      gsl_vector_scale(&v.vector, sqrt(wi));
      gsl_vector_set(ybar, i, gsl_vector_get(y, i) * sqrt(wi));
    }

  /* Xf = X^T W X */
  gsl_blas_dgemm(CblasTrans,
                 CblasNoTrans,
                 1.0,
                 Xbar,
                 Xbar,
                 0.0,
                 Xf);

  /* yf = X^T W y */
  gsl_blas_dgemv(CblasTrans,
                 1.0,
                 Xbar,
                 ybar,
                 0.0,
                 yf);

  /*
   * Now solve: Xf c = yf
   *
   * Xf is symmetric and positive (semi) definite
   */

#if 0
  {
    gsl_eigen_symmv_workspace *symm = gsl_eigen_symmv_alloc(Xf->size1);
    gsl_vector *eval = gsl_vector_alloc(Xf->size1);
    gsl_matrix *evec = gsl_matrix_alloc(Xf->size1, Xf->size2);

    gsl_eigen_symmv(Xf, eval, evec, symm);

    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

    for (i = 0; i < eval->size; ++i)
      {
        printf("%zu %f\n", i, log(gsl_vector_get(eval, i)));
      }

    gsl_eigen_symmv_free(symm);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);

    exit(1);
  }
#endif

  s += gsl_linalg_cholesky_decomp(Xf);
  s += gsl_linalg_cholesky_solve(Xf, yf, c);

  gsl_matrix_free(Xbar);
  gsl_matrix_free(Xf);
  gsl_vector_free(ybar);
  gsl_vector_free(yf);

  return s;
} /* fit2_wproc() */
