/*
 * lls_complex.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>

#include "lls.h"

/*
lls_complex_alloc()
  Allocate lls workspace

Inputs: max_block - maximum number of observations which will be
                    added to LS system at a time
        p         - number of coefficients (columns of LS design matrix)
*/

lls_complex_workspace *
lls_complex_alloc(const size_t max_block, const size_t p)
{
  const gsl_multifit_robust_type *robust_t = gsl_multifit_robust_huber;
  const size_t nblock = GSL_MAX(max_block, p);
  lls_complex_workspace *w;

  w = calloc(1, sizeof(lls_complex_workspace));
  if (!w)
    return 0;

  w->AHA = gsl_matrix_complex_alloc(p, p);
  w->work_A_complex = gsl_matrix_complex_alloc(p, p);
  if (!w->AHA || !w->work_A_complex)
    {
      lls_complex_free(w);
      return 0;
    }

  w->c = gsl_vector_alloc(p);
  w->AHb = gsl_vector_complex_alloc(p);
  w->work_b_complex = gsl_vector_complex_alloc(p);
  w->AF = gsl_matrix_complex_alloc(p, p);
  w->S = gsl_vector_alloc(p);
  if (!w->AHb || !w->work_b_complex)
    {
      lls_complex_free(w);
      return 0;
    }

  w->r_complex = gsl_vector_complex_alloc(nblock);
  w->c_complex = gsl_vector_complex_alloc(p);
  w->r = gsl_vector_alloc(nblock);
  w->w_robust = gsl_vector_alloc(nblock);
  w->robust_workspace_p = gsl_multifit_robust_alloc(robust_t, nblock, p);
  if (!w->r_complex || !w->w_robust || !w->c_complex)
    {
      lls_complex_free(w);
      return 0;
    }

  w->p = p;
  w->max_block = nblock;

  w->residual = 0.0;
  w->chisq = 0.0;
  w->cond = -1.0;
  w->niter = 0;
  w->bHb = 0.0;

  /* initialize matrices/vectors to 0 */
  lls_complex_reset(w);

  return w;
} /* lls_complex_alloc() */

void
lls_complex_free(lls_complex_workspace *w)
{
  if (w->c)
    gsl_vector_free(w->c);

  if (w->AHA)
    gsl_matrix_complex_free(w->AHA);

  if (w->AHb)
    gsl_vector_complex_free(w->AHb);

  if (w->work_A_complex)
    gsl_matrix_complex_free(w->work_A_complex);

  if (w->work_b_complex)
    gsl_vector_complex_free(w->work_b_complex);

  if (w->AF)
    gsl_matrix_complex_free(w->AF);

  if (w->S)
    gsl_vector_free(w->S);

  if (w->c_complex)
    gsl_vector_complex_free(w->c_complex);

  if (w->r_complex)
    gsl_vector_complex_free(w->r_complex);

  if (w->r)
    gsl_vector_free(w->r);

  if (w->w_robust)
    gsl_vector_free(w->w_robust);

  if (w->robust_workspace_p)
    gsl_multifit_robust_free(w->robust_workspace_p);

  free(w);
} /* lls_complex_free() */

/*
lls_complex_reset()
  Re-initialize matrices and vectors to 0
*/

int
lls_complex_reset(lls_complex_workspace *w)
{
  gsl_matrix_complex_set_zero(w->AHA);
  gsl_vector_complex_set_zero(w->AHb);
  w->bHb = 0.0;

  return 0;
}

/*
lls_complex_fold()
  Add a set of n observations to the current A^H W A matrix and
A^H W b rhs vector. This routine is designed so it can be called multiple
times, so the caller can read in a set of observations, add them to the LS
system, and then repeat

Inputs: A   - matrix of basis functions (size n-by-p)
              A_{ij} = g_j(x_i)
        b   - right hand side vector (size n)
        wts - (input/output)
              On input, data weights = diag(W) (size n)
              On output, final weights (including robust weight multipliers)
        w   - workspace

Notes:
1) On output,
AHA <- AHA + A^H W A
AHb <- AHb + A^H W b
bHb <- bHb + b^H W b (for computing chi^2 later)

2) On output,
A <- sqrt(W) A
b <- W b

3) On output, w->ssq_res contains ||b - A c||^2 using the coefficient
vector from the previous iteration

4) For robust iterations, w->niter and w->c_complex must be initialized
prior to calling the function
*/

int
lls_complex_fold(gsl_matrix_complex *A, gsl_vector_complex *b,
                 gsl_vector *wts, lls_complex_workspace *w)
{
  const size_t n = A->size1;

  if (A->size2 != w->p)
    {
      fprintf(stderr, "lls_complex_fold: A has wrong size2\n");
      return GSL_EBADLEN;
    }
  else if (n != b->size)
    {
      fprintf(stderr, "lls_complex_fold: b has wrong size\n");
      return GSL_EBADLEN;
    }
  else if (n != wts->size)
    {
      fprintf(stderr, "lls_complex_fold: wts has wrong size\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;
      size_t i;
      gsl_vector_view wv = gsl_vector_subvector(w->w_robust, 0, n);

      if (w->niter > 0)
        {
          gsl_vector_complex_view rc = gsl_vector_complex_subvector(w->r_complex, 0, n);
          gsl_vector_view rv = gsl_vector_subvector(w->r, 0, n);

          /* calculate residuals with previously computed coefficients: r = b - A c */
          gsl_vector_complex_memcpy(&rc.vector, b);
          gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_NEGONE, A, w->c_complex, GSL_COMPLEX_ONE, &rc.vector);

          /* compute Re(r) */
          for (i = 0; i < n; ++i)
            {
              gsl_complex ri = gsl_vector_complex_get(&rc.vector, i);
              gsl_vector_set(&rv.vector, i, GSL_REAL(ri));
            }

          /* calculate weights with robust weighting function */
          gsl_multifit_robust_weights(&rv.vector, &wv.vector, w->robust_workspace_p);
        }
      else
        gsl_vector_set_all(&wv.vector, 1.0);

      /* compute final weights as product of input and robust weights */
      gsl_vector_mul(wts, &wv.vector);

      for (i = 0; i < n; ++i)
        {
          gsl_vector_complex_view rv = gsl_matrix_complex_row(A, i);
          gsl_complex bi = gsl_vector_complex_get(b, i);
          double wi = gsl_vector_get(wts, i);
          gsl_complex val;

          GSL_SET_COMPLEX(&val, sqrt(wi), 0.0);

          /* A <- sqrt(W) A */
          gsl_vector_complex_scale(&rv.vector, val);

          /* b <- W b */
          val = gsl_complex_mul_real(bi, wi);
          gsl_vector_complex_set(b, i, val);

          /* bHb += b^H W b */
          w->bHb += wi * gsl_complex_abs2(bi);
        }
 
      /* AHA += A^H W A, using only the upper half of the matrix */
      s = gsl_blas_zherk(CblasUpper, CblasConjTrans, 1.0, A, 1.0, w->AHA);
      if (s)
        return s;

      /* AHb += A^H W b */
      s = gsl_blas_zgemv(CblasConjTrans, GSL_COMPLEX_ONE, A, b, GSL_COMPLEX_ONE, w->AHb);
      if (s)
        return s;

      return s;
    }
} /* lls_complex_fold() */

/*
lls_complex_solve()
  Solve the least squares system:

A^H A c = A^H b

where A^H A and A^H b have been previously computed with
lls_complex_fold

Inputs: c - (output) coefficient vector
        w - workspace

Notes:

1) on output, the residual || A^H A c - A^H b || is stored in
w->residual

2) on output, the residual ssq is stored in w->chisq:
chi^2 = b^H W b - 2 c^H A^H W b + c^H A^H W A c

3) on output, the matrix condition number is stored in w->cond

4) on output, w->niter is incremented for possible future robust
iterations
*/

int
lls_complex_solve(gsl_vector_complex *c, lls_complex_workspace *w)
{
  if (c->size != w->p)
    {
      fprintf(stderr, "lls_complex_solve: coefficient vector has wrong size\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;

      /* solve AHA c = AHb and estimate condition number */
      s = lls_lapack_zposv(c, w);

      /* compute residual || AHA c - AHb || */
      gsl_vector_complex_memcpy(w->work_b_complex, w->AHb);
      gsl_blas_zhemv(CblasUpper, GSL_COMPLEX_ONE, w->AHA, c, GSL_COMPLEX_NEGONE, w->work_b_complex);
      w->residual = gsl_blas_dznrm2(w->work_b_complex);

      /* compute chi^2 = b^H W b - 2 c^H A^H W b + c^H A^H W A c */
      {
        gsl_complex negtwo = gsl_complex_rect(-2.0, 0.0);
        gsl_complex val;

        /* compute: AHA c - 2 AHb */
        gsl_vector_complex_memcpy(w->work_b_complex, w->AHb);
        gsl_blas_zhemv(CblasUpper, GSL_COMPLEX_ONE, w->AHA, c, negtwo, w->work_b_complex);

        /* compute: c^H ( AHA c - 2 AHb ) */
        gsl_blas_zdotc(c, w->work_b_complex, &val);

        w->chisq = w->bHb + GSL_REAL(val);
      }

      /* save coefficient vector for future robust iterations */
      gsl_vector_complex_memcpy(w->c_complex, c);

      ++(w->niter);

      return s;
    }
} /* lls_complex_solve() */

/*
lls_complex_regularize2()
  Perform Tikhonov regularization on the least-squares
problem by adding a positive constant to the diagonal
matrix terms

x = (A^H A + D^T D)^{-1) A^H b

Inputs: diag - diag(D); real valued
        w    - workspace
*/

int
lls_complex_regularize2(const gsl_vector *diag, lls_complex_workspace *w)
{
  int s = 0;
  size_t n = w->AHA->size1;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double di = gsl_vector_get(diag, i);
      gsl_complex z = gsl_matrix_complex_get(w->AHA, i, i);

      GSL_REAL(z) += di * di;
      gsl_matrix_complex_set(w->AHA, i, i, z);
    }

  return s;
} /* lls_complex_regularize2() */

/* compute B = (A^H A)^{-1} */
int
lls_complex_invert(gsl_matrix_complex *B, const lls_complex_workspace *w)
{
  size_t n = w->AHA->size1;

  if (B->size1 != n || B->size2 != n)
    {
      fprintf(stderr, "lls_complex_invert: B has wrong dimensions\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s;
      size_t i, j;

      s = lls_lapack_zinvert(B, w);
      if (s)
        fprintf(stderr, "lls_complex_invert: error inverting A^H A: %d\n", s);

      /* fill in lower half of B */
      for (j = 1; j < n; ++j)
        {
          for (i = 0; i < j; ++i)
            {
              gsl_complex z = gsl_matrix_complex_get(B, i, j);
              GSL_SET_IMAG(&z, -GSL_IMAG(z));
              gsl_matrix_complex_set(B, j, i, z); 
            }
        }

      return s;
    }
} /* lls_complex_invert() */

/*
lls_complex_correlation()
  Compute correlation matrix:
  
B = diag(C)^{-1/2} C diag(C)^{-1/2}

where C is the covariance matrix: C = (A^H A)^{-1}
*/

int
lls_complex_correlation(gsl_matrix_complex *B, const lls_complex_workspace *w)
{
  size_t n = w->AHA->size1;

  if (B->size1 != n || B->size2 != n)
    {
      fprintf(stderr, "lls_complex_correlation: B has wrong dimensions\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s;
      size_t i;
      gsl_vector_complex_view d = gsl_matrix_complex_diagonal(B);

      /* compute covariance matrix */
      s = lls_complex_invert(B, w);
      if (s)
        {
          fprintf(stderr, "lls_complex_correlation: error computing covariance matrix: %d\n", s);
          return s;
        }

      /* compute diag(C)^{-1/2} C diag(C)^{-1/2} */
      for (i = 0; i < n; ++i)
        {
          gsl_complex di = gsl_vector_complex_get(&d.vector, i);
          gsl_vector_complex_view ri = gsl_matrix_complex_row(B, i);
          gsl_vector_complex_view ci = gsl_matrix_complex_column(B, i);
          gsl_complex z;

          GSL_SET_COMPLEX(&z, 1.0 / sqrt(GSL_REAL(di)), 0.0);
          gsl_vector_complex_scale(&ri.vector, z);
          gsl_vector_complex_scale(&ci.vector, z);
        }

      return s;
    }
} /* lls_complex_correlation() */

/*
lls_complex_save()
  Save matrices and rhs vectors to a binary file
*/

int
lls_complex_save(const char *filename, lls_complex_workspace *w)
{
  int s = 0;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "lls_complex_save: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  fwrite(&(w->bHb), sizeof(double), 1, fp);
  s += gsl_matrix_complex_fwrite(fp, w->AHA);
  s += gsl_vector_complex_fwrite(fp, w->AHb);

  fclose(fp);

  return s;
} /* lls_complex_save() */

/*
lls_complex_load()
  Load matrices and rhs vectors from binary file
*/

int
lls_complex_load(const char *filename, lls_complex_workspace *w)
{
  int s = 0;
  FILE *fp;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "lls_complex_load: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  fread(&(w->bHb), sizeof(double), 1, fp);
  s += gsl_matrix_complex_fread(fp, w->AHA);
  s += gsl_vector_complex_fread(fp, w->AHb);

  fclose(fp);

  return s;
} /* lls_complex_load() */
