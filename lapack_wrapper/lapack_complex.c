/*
 * lapack_complex.c
 *
 * Wrappers for complex LAPACK routines
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <lapacke/lapacke.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

#include "lapack_wrapper.h"

/*
lapack_complex_cholesky_decomp()
  Compute Cholesky decomposition of Hermitian square matrix A.

Inputs: A - on input, Hermitian matrix, stored in upper triangle
            on output, Cholesky factor stored in upper triangle

Return: success/error
*/

int
lapack_complex_cholesky_decomp(gsl_matrix_complex * A)
{
  if (A->size1 != A->size2)
    {
      GSL_ERROR("matrix must be square", GSL_ENOTSQR);
    }
  else
    {
      int s;
      lapack_int N = A->size1;
      lapack_int lda = A->size1;

      gsl_matrix_complex_transpose(A);

      /* compute Cholesky decomposition of A */
      s = LAPACKE_zpotrf(LAPACK_COL_MAJOR,
                         'U',
                         N,
                         (lapack_complex_double *) A->data,
                         lda);
      if (s)
        {
          fprintf(stderr, "lapack_complex_cholesky_decomp: error in ZPOTRF: %d\n", s);
          return s;
        }

      gsl_matrix_complex_transpose(A);

      return s;
    }
}

/*
lapack_complex_cholesky_invert()
  Compute inverse of Hermitian matrix A, given its Cholesky decomposition

Inputs: A - on input, Cholesky factor, stored in upper triangle
            on output, inverse matrix A^{-1} stored in upper triangle

Return: success/error
*/

int
lapack_complex_cholesky_invert(gsl_matrix_complex * A)
{
  if (A->size1 != A->size2)
    {
      GSL_ERROR("matrix must be square", GSL_ENOTSQR);
    }
  else
    {
      int s;
      lapack_int N = A->size1;
      lapack_int lda = A->size1;

      gsl_matrix_complex_transpose(A);

      /* compute Cholesky decomposition of A */
      s = LAPACKE_zpotri(LAPACK_COL_MAJOR,
                         'U',
                         N,
                         (lapack_complex_double *) A->data,
                         lda);
      if (s)
        {
          fprintf(stderr, "lapack_complex_cholesky_decomp: error in ZPOTRF: %d\n", s);
          return s;
        }

      gsl_matrix_complex_transpose(A);

      return s;
    }
}

/*
lapack_complex_zposv()
  Solve A x = b for square Hermitian system and compute condition number

Inputs: b     - right hand side vector
        A     - on input, Hermitian matrix stored in upper triangle
                on output, Cholesky factor stored in upper triangle
        x     - (output) solution vector
        rcond - (output) reciprocal condition number

Return: success/error
*/

int
lapack_complex_zposv(const gsl_vector_complex * b, gsl_matrix_complex * A,
                     gsl_vector_complex *x, double * rcond)
{
  if (A->size1 != A->size2)
    {
      GSL_ERROR("matrix must be square", GSL_ENOTSQR);
    }
  else
    {
      int s = 0;
      lapack_int N = A->size1;
      lapack_int nrhs = 1;
      lapack_int lda = A->size1;
      lapack_int ldb = b->size;
      lapack_int ldx = x->size;
      gsl_matrix_complex *AF = gsl_matrix_complex_alloc(N, N);
      gsl_vector_complex *work_b = gsl_vector_complex_alloc(N);
      gsl_vector *S = gsl_vector_alloc(N);
      double ferr, berr;
      lapack_int ldaf = N;
      char equed;
      size_t j;

      gsl_matrix_complex_transpose(A);
      gsl_vector_complex_memcpy(work_b, b);

      /* use expert driver to get condition number estimate */
      s = LAPACKE_zposvx(LAPACK_COL_MAJOR,
                         'E',
                         'U',
                         N,
                         nrhs,
                         (lapack_complex_double *) A->data,
                         lda,
                         (lapack_complex_double *) AF->data,
                         ldaf,
                         &equed,
                         S->data,
                         (lapack_complex_double *) work_b->data,
                         ldb,
                         (lapack_complex_double *) x->data,
                         ldx,
                         rcond,
                         &ferr,
                         &berr);

      /* save Cholesky factor (of equilibrated matrix S A S) in A */
      gsl_matrix_complex_transpose_memcpy(A, AF);

      if (equed == 'Y')
        {
          /* undo scaling transformation to recover Cholesky factor of A,
           * R = R~ S^{-1} */
          for (j = 0; j < (size_t) N; ++j)
            {
              double sj = gsl_vector_get(S, j);
              gsl_vector_complex_view v = gsl_matrix_complex_subcolumn(A, j, 0, j + 1);
              gsl_blas_zdscal(1.0 / sj, &v.vector);
            }
        }

      gsl_matrix_complex_free(AF);
      gsl_vector_complex_free(work_b);
      gsl_vector_free(S);

      return s;
    }
}
