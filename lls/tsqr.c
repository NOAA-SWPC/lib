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

static int tsqr_zero_R(gsl_matrix *R);
static double tsqr_householder_transform (const size_t N, const size_t j,
                                          gsl_vector * v);
static int tsqr_householder_hv (const size_t N, const size_t colidx, const double tau,
                                const gsl_vector * v, gsl_vector * w);
static int tsqr_QR_decomp (gsl_matrix * A, gsl_vector * tau);

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

      w->init = 1;

      return GSL_SUCCESS;
    }
  else
    {
      int status;
      const size_t npp = n + p;
      gsl_vector_view tau = gsl_vector_subvector(w->tau, 0, p);
      gsl_matrix_view R = gsl_matrix_submatrix(w->R, 0, 0, npp, p);
      gsl_vector_view QTb = gsl_vector_subvector(w->QTb, 0, npp);
      gsl_matrix_view Ai = gsl_matrix_submatrix(w->R, p, 0, n, p);
      gsl_vector_view bi = gsl_vector_subvector(w->QTb, p, n);

      /* form w->R = [ R_{i-1} ; A_i ] */
      gsl_matrix_memcpy(&Ai.matrix, A);

      status = tsqr_QR_decomp(&R.matrix, &tau.vector);
      if (status)
        return status;

      {
        size_t i;

        /* compute Q^T [ QTb_{i - 1}; b_i ], accounting for the sparse
         * structure of the Householder reflectors */
        gsl_vector_memcpy(&bi.vector, b);
        for (i = 0; i < p; i++)
          {
            gsl_vector_const_view h = gsl_matrix_const_subcolumn (&R.matrix, i, i, npp - i);
            gsl_vector_view w = gsl_vector_subvector (&QTb.vector, i, npp - i);
            double ti = gsl_vector_get (&tau.vector, i);
            tsqr_householder_hv (p, i, ti, &(h.vector), &(w.vector));
          }
      }

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

  /* XXX: zero R below the diagonal for SVD routine - need to develop
   * SVD routine for triangular matrices? */
  tsqr_zero_R(w->R);
  status = gsl_multifit_linear_svd(&R.matrix, w->multifit_workspace_p);

  return status;
}

/*
tsqr_solve()
*/

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

/*
tsqr_householder_transform()
  This routine is an optimized version of
gsl_linalg_householder_transform(), designed for the QR
decomposition of M-by-N matrices of the form:

T = [ R ]
    [ A ]

where R is N-by-N upper triangular, and A is (M-N)-by-N dense.
This routine computes a householder transformation (tau,v) of a 
x so that P x = [ I - tau*v*v' ] x annihilates x(1:n-1). x will
be a subcolumn of the matrix T, and so its structure will be:

x = [ * ] <- 1 nonzero value
    [ 0 ] <- N - j - 1 zeros, where j is column of matrix in [0,N-1]
    [ * ] <- M-N nonzero values for the dense part A

Inputs: N      - number of columns in matrix T
        colidx - column number in [0, N-1]
        v      - on input, x vector
                 on output, householder vector v
*/

static double
tsqr_householder_transform (const size_t N, const size_t colidx, gsl_vector * v)
{
  /* replace v[0:M-1] with a householder vector (v[0:M-1]) and
     coefficient tau that annihilate v[1:M-1] */

  const size_t M = v->size ;

  if (M == 1)
    {
      return 0.0; /* tau = 0 */
    }
  else
    { 
      double alpha, beta, tau ;
      
      /* A portion of vector */
      gsl_vector_view x = gsl_vector_subvector (v, N - colidx, M - (N - colidx));

      /* compute xnorm = || v[1:M-1] ||, ignoring zero part of vector */
      double xnorm = gsl_blas_dnrm2(&x.vector);

      if (xnorm == 0) 
        {
          return 0.0; /* tau = 0 */
        }
      
      alpha = gsl_vector_get (v, 0) ;
      beta = - (alpha >= 0.0 ? +1.0 : -1.0) * hypot(alpha, xnorm) ;
      tau = (beta - alpha) / beta ;
      
      {
        double s = (alpha - beta);
        
        if (fabs(s) > GSL_DBL_MIN) 
          {
            gsl_blas_dscal (1.0 / s, &x.vector);
            gsl_vector_set (v, 0, beta) ;
          }
        else
          {
            gsl_blas_dscal (GSL_DBL_EPSILON / s, &x.vector);
            gsl_blas_dscal (1.0 / GSL_DBL_EPSILON, &x.vector);
            gsl_vector_set (v, 0, beta) ;
          }
      }
      
      return tau;
    }
}

static int
tsqr_householder_hv (const size_t N, const size_t colidx, const double tau,
                     const gsl_vector * v, gsl_vector * w)
{
  /* applies a householder transformation v to vector w */
  const size_t M = v->size;
 
  if (tau == 0)
    return GSL_SUCCESS ;

  {
    /* compute d = v'w */

    double w0 = gsl_vector_get(w,0);
    double d1, d;

    gsl_vector_const_view v1 = gsl_vector_const_subvector(v, N - colidx, M - (N - colidx));
    gsl_vector_view w1 = gsl_vector_subvector(w, N - colidx, M - (N - colidx));

    /* compute d1 = v(2:n)'w(2:n) */
    gsl_blas_ddot (&v1.vector, &w1.vector, &d1);

    /* compute d = v'w = w(1) + d1 since v(1) = 1 */
    d = w0 + d1;

    /* compute w = w - tau (v) (v'w) */

    gsl_vector_set (w, 0, w0 - tau * d);
    gsl_blas_daxpy (-tau * d, &v1.vector, &w1.vector);
  }
  
  return GSL_SUCCESS;
}

static int
tsqr_householder_hm (const size_t N, const size_t colidx, const double tau, const gsl_vector * v, gsl_matrix * A)
{
  /* applies a householder transformation v,tau to matrix m */

  if (tau == 0.0)
    {
      return GSL_SUCCESS;
    }
  else
    {
      gsl_vector_const_view v1 = gsl_vector_const_subvector (v, N - colidx, v->size - (N - colidx));
      gsl_matrix_view A1 = gsl_matrix_submatrix (A, N - colidx, 0, A->size1 - (N - colidx), A->size2);
      size_t j;

      for (j = 0; j < A->size2; j++)
        {
          double A0j = gsl_matrix_get (A, 0, j);
          double wj;
          gsl_vector_view A1j = gsl_matrix_column(&A1.matrix, j);

          gsl_blas_ddot (&A1j.vector, &v1.vector, &wj);
          wj += A0j;

          gsl_matrix_set (A, 0, j, A0j - tau *  wj);

          gsl_blas_daxpy (-tau * wj, &v1.vector, &A1j.vector);
        }

      return GSL_SUCCESS;
    }
}

/* Factorise a general M x N matrix A into
 *  
 *   A = Q R
 *
 * where Q is orthogonal (M x M) and R is upper triangular (M x N).
 *
 * Q is stored as a packed set of Householder transformations in the
 * strict lower triangular part of the input matrix.
 *
 * R is stored in the diagonal and upper triangle of the input matrix.
 *
 * The full matrix for Q can be obtained as the product
 *
 *       Q = Q_k .. Q_2 Q_1
 *
 * where k = MIN(M,N) and
 *
 *       Q_i = (I - tau_i * v_i * v_i')
 *
 * and where v_i is a Householder vector
 *
 *       v_i = [1, m(i+1,i), m(i+2,i), ... , m(M,i)]
 *
 * This storage scheme is the same as in LAPACK.  */

static int
tsqr_QR_decomp (gsl_matrix * A, gsl_vector * tau)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (tau->size != GSL_MIN (M, N))
    {
      GSL_ERROR ("size of tau must be MIN(M,N)", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      for (i = 0; i < GSL_MIN (M, N); i++)
        {
          /* Compute the Householder transformation to reduce the j-th
             column of the matrix to a multiple of the j-th unit vector,
             taking into account the sparse structure of A */

          gsl_vector_view c = gsl_matrix_subcolumn (A, i, i, M - i);
          double tau_i = tsqr_householder_transform(N, i, &c.vector);

          gsl_vector_set (tau, i, tau_i);

          /* Apply the transformation to the remaining columns and
             update the norms */

          if (i + 1 < N)
            {
              gsl_matrix_view m = gsl_matrix_submatrix (A, i, i + 1, M - i, N - (i + 1));
              tsqr_householder_hm (N, i, tau_i, &(c.vector), &(m.matrix));
            }
        }

      return GSL_SUCCESS;
    }
}
