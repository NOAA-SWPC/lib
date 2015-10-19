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
#include <gsl/gsl_eigen.h>

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
  w->work_A = gsl_matrix_complex_alloc(p, p);
  if (!w->AHA || !w->work_A)
    {
      lls_complex_free(w);
      return 0;
    }

  w->AHb = gsl_vector_complex_alloc(p);
  w->work_b = gsl_vector_complex_alloc(p);
  w->AF = gsl_matrix_complex_alloc(p, p);
  w->S = gsl_vector_alloc(p);
  if (!w->AHb || !w->work_b)
    {
      lls_complex_free(w);
      return 0;
    }

  w->r_complex = gsl_vector_complex_alloc(nblock);
  w->c = gsl_vector_complex_alloc(p);
  w->r = gsl_vector_alloc(nblock);
  w->w_robust = gsl_vector_alloc(nblock);
  w->robust_workspace_p = gsl_multifit_robust_alloc(robust_t, nblock, p);
  if (!w->r_complex || !w->w_robust || !w->c)
    {
      lls_complex_free(w);
      return 0;
    }

  w->eigen_p = gsl_eigen_herm_alloc(p);
  w->eval = gsl_vector_alloc(p);

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
  if (w->AHA)
    gsl_matrix_complex_free(w->AHA);

  if (w->AHb)
    gsl_vector_complex_free(w->AHb);

  if (w->work_A)
    gsl_matrix_complex_free(w->work_A);

  if (w->work_b)
    gsl_vector_complex_free(w->work_b);

  if (w->AF)
    gsl_matrix_complex_free(w->AF);

  if (w->S)
    gsl_vector_free(w->S);

  if (w->c)
    gsl_vector_complex_free(w->c);

  if (w->r_complex)
    gsl_vector_complex_free(w->r_complex);

  if (w->r)
    gsl_vector_free(w->r);

  if (w->eigen_p)
    gsl_eigen_herm_free(w->eigen_p);

  if (w->eval)
    gsl_vector_free(w->eval);

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
lls_complex_stdform()
  Transform a regularized least squares system to standard form:

  chi^2 = ||b - A*c||_W^2 + \lambda^2 ||L c||^2

A~ = sqrt(W) A inv(L)
b~ = sqrt(W) b
c~ = L c

then:

  chi^2 = ||b~ - A~ c~||^2 + \lambda^2 || c~ ||^2

Inputs: A   - (input/output) On input, LS matrix A; on output,
              transformed matrix A~
        b   - (input/output) On input, right hand side vector b;
              on output, transformed vector b~
        wts - data weights = diag(W) (size n) or NULL for unit weights
        L   - regularization matrix L = diag(l_1,...,l_p) or NULL
        w   - workspace

Notes:
1) On output,
A <- sqrt(W) A inv(L)
b <- sqrt(W) b
*/

int
lls_complex_stdform(gsl_matrix_complex *A, gsl_vector_complex *b,
                    const gsl_vector *wts, const gsl_vector *L,
                    lls_complex_workspace *w)
{
  const size_t n = A->size1;
  const size_t p = A->size2;

  if (p != w->p)
    {
      fprintf(stderr, "lls_complex_stdform: A has wrong size2\n");
      return GSL_EBADLEN;
    }
  else if (n != b->size)
    {
      fprintf(stderr, "lls_complex_stdform: b has wrong size\n");
      return GSL_EBADLEN;
    }
  else if (wts != NULL && n != wts->size)
    {
      fprintf(stderr, "lls_complex_stdform: wts has wrong size\n");
      return GSL_EBADLEN;
    }
  else if (L != NULL && p != L->size)
    {
      fprintf(stderr, "lls_complex_stdform: L has wrong size\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;
      size_t i;

      if (wts != NULL)
        {
          for (i = 0; i < n; ++i)
            {
              gsl_vector_complex_view rv = gsl_matrix_complex_row(A, i);
              gsl_complex bi = gsl_vector_complex_get(b, i);
              double wi = gsl_vector_get(wts, i);
              double sqrtwi = sqrt(wi);
              gsl_complex val;

              GSL_SET_COMPLEX(&val, sqrtwi, 0.0);

              /* A <- sqrt(W) A */
              gsl_vector_complex_scale(&rv.vector, val);

              /* b <- sqrt(W) b */
              val = gsl_complex_mul_real(bi, sqrtwi);
              gsl_vector_complex_set(b, i, val);
            }
        }

      if (L != NULL)
        {
          /* A <- sqrt(W) A L^{-1} */
          for (i = 0; i < p; ++i)
            {
              gsl_vector_complex_view cv = gsl_matrix_complex_column(A, i);
              double Li = gsl_vector_get(L, i);
              gsl_complex val;

              GSL_SET_COMPLEX(&val, 1.0 / Li, 0.0);

              gsl_vector_complex_scale(&cv.vector, val);
            }
        }

      return s;
    }
} /* lls_complex_stdform() */

/*
lls_complex_fold()
  Add a set of n observations to the current A^H A matrix and
A^H b rhs vector. This routine is designed so it can be called multiple
times, so the caller can read in a set of observations, add them to the LS
system, and then repeat

Inputs: A   - matrix of basis functions in standard form (size n-by-p)
        b   - right hand side vector in standard form(size n)
        w   - workspace

Notes:
1) On output,
AHA <- AHA + A^H A
AHb <- AHb + A^H b
bHb <- bHb + b^H b (for computing chi^2 later)

2) On output, w->ssq_res contains ||b - A c||^2 using the coefficient
vector from the previous iteration

3) For robust iterations, w->niter and w->c must be initialized
prior to calling the function
*/

int
lls_complex_fold(const gsl_matrix_complex *A, const gsl_vector_complex *b,
                 lls_complex_workspace *w)
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
  else
    {
      int s = 0;
      double bnorm;
#if 0
      size_t i;

      gsl_vector_view wv = gsl_vector_subvector(w->w_robust, 0, n);

      if (w->niter > 0)
        {
          gsl_vector_complex_view rc = gsl_vector_complex_subvector(w->r_complex, 0, n);
          gsl_vector_view rv = gsl_vector_subvector(w->r, 0, n);

          /* calculate residuals with previously computed coefficients: r = b - A c */
          gsl_vector_complex_memcpy(&rc.vector, b);
          gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_NEGONE, A, w->c, GSL_COMPLEX_ONE, &rc.vector);

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

#endif
 
      /* AHA += A^H A, using only the upper half of the matrix */
      s = gsl_blas_zherk(CblasUpper, CblasConjTrans, 1.0, A, 1.0, w->AHA);
      if (s)
        return s;

      /* AHb += A^H b */
      s = gsl_blas_zgemv(CblasConjTrans, GSL_COMPLEX_ONE, A, b, GSL_COMPLEX_ONE, w->AHb);
      if (s)
        return s;

      /* bHb += b^H b */
      bnorm = gsl_blas_dznrm2(b);
      w->bHb += bnorm * bnorm;

      return s;
    }
} /* lls_complex_fold() */

/*
lls_complex_solve()
  Solve the least squares system:

A^H A c = A^H b

where A^H A and A^H b have been previously computed with
lls_complex_fold

Inputs: lambda - regularization parameter
        c      - (output) coefficient vector
        w      - workspace

Notes:

1) on output, the residual || A^H A c - A^H b || is stored in
w->residual

2) on output, the residual ssq is stored in w->chisq:
chi^2 = b^H b - 2 c^H A^H b + c^H A^H A c

3) on output, the matrix condition number is stored in w->cond

4) on output, w->niter is incremented for possible future robust
iterations
*/

int
lls_complex_solve(const double lambda, gsl_vector_complex *c, lls_complex_workspace *w)
{
  if (c->size != w->p)
    {
      fprintf(stderr, "lls_complex_solve: coefficient vector has wrong size\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;

      /* solve (AHA + lambda^2 I) c = AHb and estimate condition number */
      s = lls_lapack_zposv(lambda, c, w);

      /* compute residual || AHA c - AHb || */
      gsl_vector_complex_memcpy(w->work_b, w->AHb);
      gsl_blas_zhemv(CblasUpper, GSL_COMPLEX_ONE, w->AHA, c, GSL_COMPLEX_NEGONE, w->work_b);
      w->residual = gsl_blas_dznrm2(w->work_b);

      /* compute chi^2 = b^H b - 2 c^H A^H b + c^H A^H A c */
      {
        gsl_complex negtwo = gsl_complex_rect(-2.0, 0.0);
        gsl_complex val;

        /* compute: AHA c - 2 AHb */
        gsl_vector_complex_memcpy(w->work_b, w->AHb);
        gsl_blas_zhemv(CblasUpper, GSL_COMPLEX_ONE, w->AHA, c, negtwo, w->work_b);

        /* compute: c^H ( AHA c - 2 AHb ) */
        gsl_blas_zdotc(c, w->work_b, &val);

        w->chisq = w->bHb + GSL_REAL(val);
      }

      /* save coefficient vector for future robust iterations */
      gsl_vector_complex_memcpy(w->c, c);

      ++(w->niter);

      return s;
    }
} /* lls_complex_solve() */

/*
lls_complex_btransform()
  Backtransform a regularized solution in standard form to the solution
of the original problem

c = inv(L) c~

Inputs: L   - diag(L) regularization matrix
        c   - (input/output) on input, standard form solution vector
              output from lls_complex_solve(); on output, solution
              vector of original problem
        w   - workspace
*/

int
lls_complex_btransform(const gsl_vector *L, gsl_vector_complex *c,
                       lls_complex_workspace *w)
{
  size_t p = c->size;

  if (L->size != p)
    {
      GSL_ERROR("L and c vectors have different sizes", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      for (i = 0; i < p; ++i)
        {
          gsl_complex ci = gsl_vector_complex_get(c, i);
          double Li = gsl_vector_get(L, i);
          gsl_complex val = gsl_complex_div_real(ci, Li);

          gsl_vector_complex_set(c, i, val);
        }

      return GSL_SUCCESS;
    }
} /* lls_complex_btransform() */

int
lls_complex_lcurve(gsl_vector *reg_param, gsl_vector *rho, gsl_vector *eta,
                   lls_complex_workspace *w)
{
  const size_t N = rho->size; /* number of points on L-curve */

  if (N != reg_param->size)
    {
      GSL_ERROR("size of reg_param and rho do not match", GSL_EBADLEN);
    }
  else if (N != eta->size)
    {
      GSL_ERROR("size of eta and rho do not match", GSL_EBADLEN);
    }
  else
    {
      int s;
      const gsl_complex negtwo = gsl_complex_rect(-2.0, 0.0);

      /* smallest regularization parameter */
      const double smin_ratio = 16.0 * GSL_DBL_EPSILON;

      double s1, sp, ratio, tmp;
      size_t i;

      /* compute eigenvalues of A^H A */
      gsl_matrix_complex_transpose_memcpy(w->work_A, w->AHA);
      s = gsl_eigen_herm(w->work_A, w->eval, w->eigen_p);
      if (s)
        return s;

      /* find largest and smallest eigenvalues */
      gsl_vector_minmax(w->eval, &sp, &s1);

      /* singular values are square roots of eigenvalues */
      s1 = sqrt(s1);
      if (sp > GSL_DBL_EPSILON)
        sp = sqrt(fabs(sp));

      tmp = GSL_MAX(sp, s1*smin_ratio);
      gsl_vector_set(reg_param, N - 1, tmp);

      /* ratio so that reg_param(1) = s(1) */
      ratio = pow(s1 / tmp, 1.0 / (N - 1.0));

      /* calculate the regularization parameters */
      for (i = N - 1; i > 0 && i--; )
        {
          double rp1 = gsl_vector_get(reg_param, i + 1);
          gsl_vector_set(reg_param, i, ratio * rp1);
        }

      for (i = 0; i < N; ++i)
        {
          double r2;
          double lambda = gsl_vector_get(reg_param, i);
          gsl_complex val;

          lls_complex_solve(lambda, w->c, w);

          /* store ||c|| */
          gsl_vector_set(eta, i, gsl_blas_dznrm2(w->c));

          /* compute: A^H A c - 2 A^H b */
          gsl_vector_complex_memcpy(w->work_b, w->AHb);
          gsl_blas_zhemv(CblasUpper, GSL_COMPLEX_ONE, w->AHA, w->c, negtwo, w->work_b);

          /* compute: c^T A^T A c - 2 c^T A^T b */
          gsl_blas_zdotc(w->c, w->work_b, &val);
          r2 = GSL_REAL(val) + w->bHb;

          gsl_vector_set(rho, i, sqrt(r2));
        }

      return GSL_SUCCESS;
    }
} /* lls_complex_lcurve() */

/*
lls_complex_regularize()
  Perform Tikhonov regularization on the least-squares
problem by adding a positive constant to the diagonal
matrix terms

x = (A^H A + lambda^2*I)^{-1) A^H b

Inputs: lambda - regularization parameter
        AHA    - A^H A matrix
*/

int
lls_complex_regularize(const double lambda, gsl_matrix_complex *AHA)
{
  int s;
  gsl_vector_complex_view d = gsl_matrix_complex_diagonal(AHA);
  gsl_complex val = gsl_complex_rect(lambda * lambda, 0.0);

  s = gsl_vector_complex_add_constant(&d.vector, val);

  return s;
} /* lls_complex_regularize() */

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
