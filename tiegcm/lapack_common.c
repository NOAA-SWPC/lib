#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <lapacke/lapacke.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

int
lapack_lls(const gsl_matrix * A, const gsl_matrix * B, gsl_matrix * X,
           int *rank)
{
  int s;
  lapack_int m = A->size1;
  lapack_int n = A->size2;
  lapack_int nrhs = B->size2;
  lapack_int lda = A->size1;
  lapack_int ldb = B->size1;
  lapack_int lrank;
  lapack_int *jpvt = malloc(n * sizeof(lapack_int));
  gsl_matrix *work_A = gsl_matrix_alloc(A->size2, A->size1);
  gsl_matrix *work_B = gsl_matrix_alloc(B->size2, B->size1);
  double rcond = 1.0e-6;

  gsl_matrix_transpose_memcpy(work_A, A);
  gsl_matrix_transpose_memcpy(work_B, B);

  s = LAPACKE_dgelsy(LAPACK_COL_MAJOR,
                     m,
                     n,
                     nrhs,
                     work_A->data,
                     lda,
                     work_B->data,
                     ldb,
                     jpvt,
                     rcond,
                     &lrank);

  /* store solution in X */
  {
    gsl_matrix_view m = gsl_matrix_submatrix(work_B, 0, 0, X->size2, X->size1);
    gsl_matrix_transpose_memcpy(X, &m.matrix);
  }

  *rank = lrank;

  gsl_matrix_free(work_A);
  gsl_matrix_free(work_B);
  free(jpvt);

  return s;
}

int
lapack_complex_lls(const gsl_matrix_complex * A, const gsl_matrix_complex * B,
                   gsl_matrix_complex * X, int *rank)
{
  int s;
  lapack_int m = A->size1;
  lapack_int n = A->size2;
  lapack_int nrhs = B->size2;
  lapack_int lda = A->size1;
  lapack_int ldb = B->size1;
  lapack_int lrank;
  lapack_int *jpvt = malloc(n * sizeof(lapack_int));
  gsl_matrix_complex *work_A = gsl_matrix_complex_alloc(A->size2, A->size1);
  gsl_matrix_complex *work_B = gsl_matrix_complex_alloc(B->size2, B->size1);
  double rcond = 1.0e-6;

  gsl_matrix_complex_transpose_memcpy(work_A, A);
  gsl_matrix_complex_transpose_memcpy(work_B, B);

  s = LAPACKE_zgelsy(LAPACK_COL_MAJOR,
                     m,
                     n,
                     nrhs,
                     (lapack_complex_double *) work_A->data,
                     lda,
                     (lapack_complex_double *) work_B->data,
                     ldb,
                     jpvt,
                     rcond,
                     &lrank);

  /* store solution in X */
  {
    gsl_matrix_complex_view m = gsl_matrix_complex_submatrix(work_B, 0, 0, X->size2, X->size1);
    gsl_matrix_complex_transpose_memcpy(X, &m.matrix);
  }

  *rank = lrank;

  gsl_matrix_complex_free(work_A);
  gsl_matrix_complex_free(work_B);
  free(jpvt);

  return s;
}

int
lapack_eigen_symm(const gsl_matrix * m, gsl_vector *eval, gsl_matrix *evec,
                  int *eval_found)
{
  int s;
  const lapack_int N = m->size1;
  gsl_matrix *A = gsl_matrix_alloc(N, N);
  double vl = 0.0, vu = 0.0;
  lapack_int il = 0, iu = 0;
  lapack_int lda = A->size1;
  lapack_int ldz = evec->size1;
  double abstol = 0.0;
  lapack_int M = 0;
  lapack_int *isuppz = malloc(2*N*sizeof(lapack_int));

  gsl_matrix_transpose_memcpy(A, m);

  s = LAPACKE_dsyevr(LAPACK_COL_MAJOR,
                     'V',
                     'A',
                     'L',
                     N,
                     A->data,
                     lda,
                     vl,
                     vu,
                     il,
                     iu,
                     abstol,
                     &M,
                     eval->data,
                     evec->data,
                     ldz,
                     isuppz);

  gsl_matrix_transpose(evec);

  /* sort into descending order */
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  *eval_found = M;

  free(isuppz);
  gsl_matrix_free(A);

  return s;
}

int
lapack_eigen_herm(const gsl_matrix_complex * m, gsl_vector *eval, gsl_matrix_complex *evec,
                  int *eval_found)
{
  int s;
  const lapack_int N = m->size1;
  gsl_matrix_complex *A = gsl_matrix_complex_alloc(N, N);
  double vl = 0.0, vu = 0.0;
  lapack_int il = 0, iu = 0;
  lapack_int lda = A->size1;
  lapack_int ldz = evec->size1;
  double abstol = 0.0;
  lapack_int M = 0;
  lapack_int *isuppz = malloc(2*N*sizeof(lapack_int));

  gsl_matrix_complex_transpose_memcpy(A, m);

  s = LAPACKE_zheevr(LAPACK_COL_MAJOR,
                     'V',
                     'A',
                     'L',
                     N,
                     (lapack_complex_double *) A->data,
                     lda,
                     vl,
                     vu,
                     il,
                     iu,
                     abstol,
                     &M,
                     eval->data,
                     (lapack_complex_double *) evec->data,
                     ldz,
                     isuppz);

  gsl_matrix_complex_transpose(evec);

  *eval_found = M;

  free(isuppz);
  gsl_matrix_complex_free(A);

  return s;
}
