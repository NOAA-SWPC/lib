/*
 * tsqr.c
 *
 * This module implements the Tall Skinny QR (TSQR) algorithm
 * for least squares problems
 *
 * 1) tsqr_accumulate: the tall matrix A is decomposed into QR, one block
 *    a a time.
 * 2) tsqr_svd: compute SVD of final R factor
 * 3) tsqr_lcurve: compute L-curve of system and determine lambda
 * 4) tsqr_solve: solve system using optimal lambda
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

  w->multifit_workspace_p = gsl_multifit_linear_alloc(p, p);

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

  if (w->multifit_workspace_p)
    gsl_multifit_linear_free(w->multifit_workspace_p);

  free(w);
}

/*
tsqr_accumulate()
  Accumulate a new matrix block into the QR system

Inputs: A - matrix n-by-p with n <= nmax
        b - right hand side vector n-by-1 with n <= nmax
        w - workspace

Return: success/error

Notes:
1) On output,
w->R(1:p,1:p) contains current R matrix and the rest of w->R is 0
w->QTb(1:p) contains current Q^T b vector
*/

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
      gsl_matrix_view R = gsl_matrix_submatrix(w->R, 0, 0, n, p);
      gsl_vector_view QTb = gsl_vector_subvector(w->QTb, 0, n);

      gsl_matrix_set_zero(w->R);
      gsl_vector_set_zero(w->QTb);

      /* copy A into the upper portion of w->R, so that R = [ A ; 0 ] */
      gsl_matrix_memcpy(&R.matrix, A);

      /* compute QR decomposition of A */
      status = gsl_linalg_QR_decomp(&R.matrix, &tau.vector);
      if (status)
        return status;

      /* compute Q^T b */
      gsl_vector_memcpy(&QTb.vector, b);
      gsl_linalg_QR_QTvec(&R.matrix, &tau.vector, &QTb.vector);

      /* XXX: zero the R matrix below the diagonal */
      tsqr_zero_R(w->R);

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

      /* form w->R = [ R_{i-1} ; A_i ] */
      gsl_matrix_memcpy(&Ai.matrix, A);

      /* XXX: this routine could be optimized to take advantage of the
       * sparse structure of R
       */
      status = gsl_linalg_QR_decomp(&R.matrix, &tau.vector);
      if (status)
        return status;

      gsl_vector_memcpy(&bi.vector, b);
      gsl_linalg_QR_QTvec(&R.matrix, &tau.vector, &QTb.vector);

      /* XXX: zero the R matrix below the diagonal */
      tsqr_zero_R(w->R);

      return GSL_SUCCESS;
    }
}

/* compute SVD of R factor */
int
tsqr_svd(tsqr_workspace * w)
{
  int status;
  const size_t p = w->p;
  gsl_matrix_view R = gsl_matrix_submatrix(w->R, 0, 0, p, p);

  status = gsl_multifit_linear_svd(&R.matrix, w->multifit_workspace_p);

  return status;
}

int
tsqr_solve(const double lambda, gsl_vector * c, double *rnorm, double *snorm,
           tsqr_workspace * w)
{
  const size_t p = w->p;

  if (p != c->size)
    {
      GSL_ERROR("c vector has wrong size", GSL_EBADLEN);
    }
  else
    {
      int status;
      gsl_matrix_view R = gsl_matrix_submatrix(w->R, 0, 0, p, p);
      gsl_vector_view QTb = gsl_vector_subvector(w->QTb, 0, p);

      status = gsl_multifit_linear_solve(lambda, &R.matrix, &QTb.vector, c, rnorm, snorm,
                                         w->multifit_workspace_p);

      return status;
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
