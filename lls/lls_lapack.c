/*
 * lls_lapack.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lapacke/lapacke.h>

#include "lls.h"

int
lls_lapack_dposv(const double lambda, gsl_vector *c, lls_workspace *w)
{
  int s = 0;
  lapack_int N = w->p;
  lapack_int nrhs = 1;
  lapack_int lda = w->work_A->size1;
  lapack_int ldb = w->work_b->size;
  lapack_int ldx = c->size;
  double rcond;
  double ferr, berr;
  lapack_int ldaf = N;
  char equed;

  gsl_matrix_transpose_memcpy(w->work_A, w->ATA);
  gsl_vector_memcpy(w->work_b, w->ATb);

  lls_regularize(lambda, w->work_A);

  /* use expert driver to get condition number estimate */
  s = LAPACKE_dposvx(LAPACK_COL_MAJOR,
                     'E',
                     'U',
                     N,
                     nrhs,
                     w->work_A->data,
                     lda,
                     w->AF->data,
                     ldaf,
                     &equed,
                     w->S->data,
                     w->work_b->data,
                     ldb,
                     c->data,
                     ldx,
                     &rcond,
                     &ferr,
                     &berr);

  if (s == 0)
    w->cond = 1.0 / rcond;

  return s;
} /* lls_lapack_dposv() */

int
lls_lapack_zposv(gsl_vector_complex *c, lls_complex_workspace *w)
{
  int s = 0;
  lapack_int n = w->work_A_complex->size1;
  lapack_int nrhs = 1;
  lapack_int lda = w->work_A_complex->size1;
  lapack_int ldb = w->work_b_complex->size;
  lapack_int ldx = c->size;
  double rcond;
  double ferr, berr;
  lapack_int ldaf = n;
  char equed = 'N';

  gsl_matrix_complex_transpose_memcpy(w->work_A_complex, w->AHA);
  gsl_vector_complex_memcpy(w->work_b_complex, w->AHb);

  /* use expert driver to get condition number estimate */
  s = LAPACKE_zposvx(LAPACK_COL_MAJOR,
                     'E',
                     'U',
                     n,
                     nrhs,
                     (lapack_complex_double *) w->work_A_complex->data,
                     lda,
                     (lapack_complex_double *) w->AF->data,
                     ldaf,
                     &equed,
                     w->S->data,
                     (lapack_complex_double *) w->work_b_complex->data,
                     ldb,
                     (lapack_complex_double *) c->data,
                     ldx,
                     &rcond,
                     &ferr,
                     &berr);

  if (s == 0)
    w->cond = 1.0 / rcond;

  return s;
} /* lls_lapack_zposv() */

/* compute B = (A^H A)^{-1} */
int
lls_lapack_zinvert(gsl_matrix_complex *B, const lls_complex_workspace *w)
{
  int s = 0;
  lapack_int n = w->work_A_complex->size1;
  lapack_int lda = w->work_A_complex->size1;

  gsl_matrix_complex_transpose_memcpy(w->work_A_complex, w->AHA);

  /* compute Cholesky decomposition of A^H A */
  s = LAPACKE_zpotrf(LAPACK_COL_MAJOR,
                     'U',
                     n,
                     (lapack_complex_double *) w->work_A_complex->data,
                     lda);
  if (s)
    {
      fprintf(stderr, "lls_lapack_zinvert: error in ZPOTRF: %d\n", s);
      return s;
    }

  /* compute matrix inverse */
  s = LAPACKE_zpotri(LAPACK_COL_MAJOR,
                     'U',
                     n,
                     (lapack_complex_double *) w->work_A_complex->data,
                     lda);
  if (s)
    {
      fprintf(stderr, "lls_lapack_zinvert: error in ZPOTRI: %d\n", s);
      return s;
    }

  /* store output in upper half of B */
  gsl_matrix_complex_transpose_memcpy(B, w->work_A_complex);

  return s;
} /* lls_lapack_zinvert() */
