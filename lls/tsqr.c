/*
 * tsqr.c
 *
 * This module implements the Tall Skinny QR (TSQR) algorithm
 * for least squares problems
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "tsqr.h"

#include "oct.h"

static int tsqr_zero_R(gsl_matrix *R);

tsqr_workspace *
tsqr_alloc(const size_t nmax, const size_t p)
{
  tsqr_workspace *w;

  w = calloc(1, sizeof(tsqr_workspace));
  if (!w)
    return 0;

  w->tau = gsl_vector_alloc(p);
  w->R = gsl_matrix_alloc(nmax + p, p);
  w->QTb = gsl_vector_alloc(nmax + p);

  w->init = 0;
  w->p = p;
  w->nmax = nmax;

  return w;
}

void
tsqr_free(tsqr_workspace *w)
{
  if (w->tau)
    gsl_vector_free(w->tau);

  if (w->R)
    gsl_matrix_free(w->R);

  if (w->QTb)
    gsl_vector_free(w->QTb);

  free(w);
}

int
tsqr_accumulate(const gsl_matrix * A, const gsl_vector * b, tsqr_workspace * w)
{
  const size_t n = A->size1;
  const size_t p = A->size2;

  if (p != w->p)
    {
      GSL_ERROR("A has wrong size2", GSL_EBADLEN);
    }
  else if (n != b->size)
    {
      GSL_ERROR("size mismatch between A and b", GSL_EBADLEN);
    }
  else if (w->init == 0)
    {
      int status;
      gsl_vector_view tau = gsl_vector_subvector(w->tau, 0, GSL_MIN(n, p));
      gsl_matrix_view A0 = gsl_matrix_submatrix(w->R, 0, 0, n, p);
      gsl_vector_view QTb = gsl_vector_subvector(w->QTb, 0, n);

      gsl_matrix_set_zero(w->R);
      gsl_vector_set_zero(w->QTb);

      /* copy A into the upper portion of w->R, so that R = [ A ; 0 ] */
      gsl_matrix_memcpy(&A0.matrix, A);

      /* compute QR decomposition of A */
      status = gsl_linalg_QR_decomp(&A0.matrix, &tau.vector);
      if (status)
        return status;

      /* compute Q^T b */
      gsl_vector_memcpy(&QTb.vector, b);
      gsl_linalg_QR_QTvec(&A0.matrix, &tau.vector, &QTb.vector);

      w->init = 1;

      return GSL_SUCCESS;
    }
  else
    {
      int status;
      gsl_vector_view tau = gsl_vector_subvector(w->tau, 0, GSL_MIN(n + p, p));
      gsl_matrix_view R = gsl_matrix_submatrix(w->R, 0, 0, n + p, p);
      gsl_vector_view QTb = gsl_vector_subvector(w->QTb, 0, n + p);
      gsl_matrix_view Ai = gsl_matrix_submatrix(w->R, p, 0, n, p);
      gsl_vector_view bi = gsl_vector_subvector(w->QTb, p, n);

      /* zero the matrix below the diagonal */
      tsqr_zero_R(w->R);

      /* form w->R = [ R_{i-1} ; A_i ] */
      gsl_matrix_memcpy(&Ai.matrix, A);

      status = gsl_linalg_QR_decomp(&R.matrix, &tau.vector);
      if (status)
        return status;

      gsl_vector_memcpy(&bi.vector, b);
      gsl_linalg_QR_QTvec(&R.matrix, &tau.vector, &QTb.vector);

      {
        R = gsl_matrix_submatrix(w->R, 0, 0, p, p);
        printtri_octave(&R.matrix, "R1");

        QTb = gsl_vector_subvector(w->QTb, 0, p);
        printv_octave(&QTb.vector, "QTb1");
      }

      return GSL_SUCCESS;
    }
}

int
tsqr_solve(const double lambda, gsl_vector * c, tsqr_workspace * w)
{
  const size_t p = w->p;

  if (p != c->size)
    {
      GSL_ERROR("c vector has wrong size", GSL_EBADLEN);
    }
  else
    {
      gsl_matrix_view R = gsl_matrix_submatrix(w->R, 0, 0, p, p);
      gsl_vector_view QTb = gsl_vector_subvector(w->QTb, 0, p);

      /* solve: R c = Q^T b */
      gsl_vector_memcpy(c, &QTb.vector);
      gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, &R.matrix, c);

      return GSL_SUCCESS;
    }
}

/* zero everything below the diagonal */
static int
tsqr_zero_R(gsl_matrix *R)
{
  const size_t n = R->size1;
  const size_t p = R->size2;
  size_t j;

  for (j = 0; j < p; ++j)
    {
      gsl_vector_view v = gsl_matrix_subcolumn(R, j, j + 1, n - j - 1);
      gsl_vector_set_zero(&v.vector);
    }

  return GSL_SUCCESS;
}
