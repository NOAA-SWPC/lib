/*
 * lapack_inverse.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work,
             int *lwork, int *info);

/* on output, A is replaced by A^{-1} */
int
lapack_inverse(gsl_matrix *A)
{
  int s = 0;
  int M = A->size1;
  int N = A->size2;
  int lda = N;
  int *ipiv;
  int lwork;
  double *work;
  double q[1];

  ipiv = malloc(N * sizeof(int));

  dgetrf_(&M, &N, A->data, &lda, ipiv, &s);
  if (s != 0)
    {
      fprintf(stderr, "lapack_inverse: error: %d\n", s);
      return s;
    }

  lwork = -1;
  dgetri_(&N, A->data, &lda, ipiv, q, &lwork, &s);

  lwork = (int) q[0];
  work = malloc(lwork * sizeof(double));

  /* compute inverse */
  dgetri_(&N, A->data, &lda, ipiv, work, &lwork, &s);

  free(ipiv);
  free(work);

  return s;
}
