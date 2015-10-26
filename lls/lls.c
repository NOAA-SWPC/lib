/*
 * lls.c
 * Patrick Alken
 *
 * This module is designed to solve linear least squares problems with
 * huge numbers of observations (for example a satellite data set spanning
 * years with 1-second data). Such data-sets would lead to a least-squares
 * matrix A too large to fit into memory. Therefore, as the dataset is
 * read in a block at a time, a block of rows of the matrix A can be
 * computed, and the matrix A^T A can be calculated as a running total
 * of all observations read in thus far.
 *
 * For example, if 1000 observations are read in, 1000 rows of A can
 * be computed, and these rows can be folded in to the computation of
 * A^T A as follows: Let the matrix M = A^T A, and let R represent
 * an r-by-p submatrix of A for a block of r observations:
 *
 * M_0 = 0
 * M_{k+1} = M_k + R^T R
 *
 * Similarly let v = A^T b where b is the right hand side vector. Then:
 *
 * v_0 = 0
 * v_{k+1} = v_k + R^T b
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multifit.h>

#include "common.h"
#include "lls.h"

/*
lls_alloc()
  Allocate lls workspace

Inputs: max_block - maximum number of observations which will be
                    added to LS system at a time
        p         - number of coefficients (columns of LS design matrix)
*/

lls_workspace *
lls_alloc(const size_t max_block, const size_t p)
{
  const gsl_multifit_robust_type *robust_t = gsl_multifit_robust_huber;
  const size_t nblock = GSL_MAX(max_block, p);
  lls_workspace *w;

  w = calloc(1, sizeof(lls_workspace));
  if (!w)
    return 0;

  /* A^T A matrix */
  w->ATA = gsl_matrix_alloc(p, p);
  w->work_A = gsl_matrix_alloc(p, p);
  if (!w->ATA || !w->work_A)
    {
      lls_free(w);
      return 0;
    }

  w->c = gsl_vector_alloc(p);
  w->ATb = gsl_vector_alloc(p);
  w->work_b = gsl_vector_alloc(p);
  w->AF = gsl_matrix_alloc(p, p);
  w->S = gsl_vector_alloc(p);
  if (!w->ATb || !w->work_b)
    {
      lls_free(w);
      return 0;
    }

  w->r = gsl_vector_alloc(nblock);
  w->w_robust = gsl_vector_alloc(nblock);
  w->robust_workspace_p = gsl_multifit_robust_alloc(robust_t, nblock, p);
  if (!w->w_robust || !w->r)
    {
      lls_free(w);
      return 0;
    }

  w->eigen_p = gsl_eigen_symm_alloc(p);
  w->eval = gsl_vector_alloc(p);

  w->p = p;
  w->max_block = nblock;

  w->residual = 0.0;
  w->chisq = 0.0;
  w->cond = -1.0;
  w->niter = 0;
  w->bTb = 0.0;

  /* initialize matrices/vectors to 0 */
  lls_reset(w);

  return w;
} /* lls_alloc() */

void
lls_free(lls_workspace *w)
{
  if (w->ATA)
    gsl_matrix_free(w->ATA);

  if (w->ATb)
    gsl_vector_free(w->ATb);

  if (w->work_A)
    gsl_matrix_free(w->work_A);

  if (w->work_b)
    gsl_vector_free(w->work_b);

  if (w->c)
    gsl_vector_free(w->c);

  if (w->AF)
    gsl_matrix_free(w->AF);

  if (w->S)
    gsl_vector_free(w->S);

  if (w->r)
    gsl_vector_free(w->r);

  if (w->w_robust)
    gsl_vector_free(w->w_robust);

  if (w->eigen_p)
    gsl_eigen_symm_free(w->eigen_p);

  if (w->eval)
    gsl_vector_free(w->eval);

  if (w->robust_workspace_p)
    gsl_multifit_robust_free(w->robust_workspace_p);

  free(w);
} /* lls_free() */

/*
lls_reset()
  Re-initialize matrices and vectors to 0
*/

int
lls_reset(lls_workspace *w)
{
  gsl_matrix_set_zero(w->ATA);
  gsl_vector_set_zero(w->ATb);
  w->bTb = 0.0;

  return 0;
}

/*
lls_fold()
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
ATA <- ATA + A^T W A
ATb <- ATb + A^T W b
bTb <- bTb + b^T W b

2) On output,
A <- sqrt(W) A
b <- W b
*/

int
lls_fold(gsl_matrix *A, gsl_vector *b,
         gsl_vector *wts, lls_workspace *w)
{
  const size_t n = A->size1;

  if (A->size2 != w->p)
    {
      GSL_ERROR("A has wrong size2", GSL_EBADLEN);
    }
  else if (n != b->size)
    {
      GSL_ERROR("b has wrong size", GSL_EBADLEN);
    }
  else if (n != wts->size)
    {
      GSL_ERROR("wts has wrong size", GSL_EBADLEN);
    }
  else
    {
      int s = 0;
      size_t i;
      double bnorm;

      for (i = 0; i < n; ++i)
        {
          gsl_vector_view rv = gsl_matrix_row(A, i);
          double *bi = gsl_vector_ptr(b, i);
          double wi = gsl_vector_get(wts, i);
          double swi = sqrt(wi);

          /* A <- sqrt(W) A */
          gsl_vector_scale(&rv.vector, swi);

          /* b <- sqrt(W) b */
          *bi *= swi;
        }
 
      /* ATA += A^T W A, using only the upper half of the matrix */
      s = gsl_blas_dsyrk(CblasUpper, CblasTrans, 1.0, A, 1.0, w->ATA);
      if (s)
        return s;

      /* ATb += A^T W b */
      s = gsl_blas_dgemv(CblasTrans, 1.0, A, b, 1.0, w->ATb);
      if (s)
        return s;

      /* bTb += b^T W b */
      bnorm = gsl_blas_dnrm2(b);
      w->bTb += bnorm * bnorm;

      return s;
    }
} /* lls_fold() */

/*
lls_regularize()
  Perform Tikhonov regularization on the least-squares
problem by adding a positive constant to the diagonal
matrix terms

x = (A^T A + lambda^2*I)^{-1) A^T b

Inputs: lambda - regularization parameter
        ATA    - A^T W A matrix
*/

int
lls_regularize(const double lambda, gsl_matrix *ATA)
{
  int s;
  gsl_vector_view d = gsl_matrix_diagonal(ATA);

  s = gsl_vector_add_constant(&d.vector, lambda * lambda);

  return s;
} /* lls_regularize() */

/*
lls_regularize2()
  Perform Tikhonov regularization on the least-squares
problem by adding a positive constant to the diagonal
matrix terms

x = (A^T A + D^T D)^{-1) A^T b

Inputs: diag - diag(D)
        w    - workspace
*/

int
lls_regularize2(const gsl_vector *diag, lls_workspace *w)
{
  int s = 0;
  size_t n = w->ATA->size1;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double di = gsl_vector_get(diag, i);
      double *Aii = gsl_matrix_ptr(w->ATA, i, i);

      *Aii += di;
    }

  return s;
} /* lls_regularize2() */

/*
lls_solve()
  Solve the least squares system:

A^T A c = A^T b

where A^T A and A^T b have been previously computed with
lls_calc_ATA and lls_calc_ATb

Inputs: lambda - regularization parameter
        c      - (output) coefficient vector
        w      - workspace

Notes: on output, the residual || A^T A c - A^T b || is stored in
w->residual
*/

int
lls_solve(const double lambda, gsl_vector *c, lls_workspace *w)
{
  if (c->size != w->p)
    {
      fprintf(stderr, "lls_solve: coefficient vector has wrong size\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;

      s = lls_lapack_dposv(lambda, c, w); /* matrix solve */

      /* compute residual || ATA c - ATb || */
      gsl_vector_memcpy(w->work_b, w->ATb);
      gsl_blas_dgemv(CblasNoTrans, 1.0, w->ATA, c, -1.0, w->work_b);
      w->residual = gsl_blas_dnrm2(w->work_b);

      return s;
    }
} /* lls_solve() */

int
lls_lcurve(gsl_vector *reg_param, gsl_vector *rho, gsl_vector *eta,
           lls_workspace *w)
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
      double smax, smin;
      size_t i;

      /* compute eigenvalues of A^T W A */
      gsl_matrix_transpose_memcpy(w->work_A, w->ATA);
      s = gsl_eigen_symm(w->work_A, w->eval, w->eigen_p);
      if (s)
        return s;

      /* find largest and smallest eigenvalues */
      gsl_vector_minmax(w->eval, &smin, &smax);

      /* singular values are square roots of eigenvalues */
      smax = sqrt(smax);
      if (smin > GSL_DBL_EPSILON)
        smin = sqrt(fabs(smin));

      gsl_multifit_linear_lreg(smin, smax, reg_param);

      for (i = 0; i < N; ++i)
        {
          double r2;
          double lambda = gsl_vector_get(reg_param, i);

          lls_solve(lambda, w->c, w);

          /* store ||c|| */
          gsl_vector_set(eta, i, gsl_blas_dnrm2(w->c));

          /* compute: A^T A c - 2 A^T y */
          gsl_vector_memcpy(w->work_b, w->ATb);
          gsl_blas_dsymv(CblasUpper, 1.0, w->ATA, w->c, -2.0, w->work_b);

          /* compute: c^T A^T A c - 2 c^T A^T y */
          gsl_blas_ddot(w->c, w->work_b, &r2);

          r2 += w->bTb;
          gsl_vector_set(rho, i, sqrt(r2));
        }

      return GSL_SUCCESS;
    }
} /* lls_lcurve() */

/*
lls_save()
  Save matrices and rhs vectors to a binary file
*/

int
lls_save(const char *filename, lls_workspace *w)
{
  int s = 0;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "lls_save: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  fwrite(&(w->bTb), sizeof(double), 1, fp);
  s += gsl_matrix_fwrite(fp, w->ATA);
  s += gsl_vector_fwrite(fp, w->ATb);

  fclose(fp);

  return s;
} /* lls_save() */

/*
lls_load()
  Load matrices and rhs vectors from binary file
*/

int
lls_load(const char *filename, lls_workspace *w)
{
  int s = 0;
  FILE *fp;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "lls_load: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  fread(&(w->bTb), sizeof(double), 1, fp);
  s += gsl_matrix_fread(fp, w->ATA);
  s += gsl_vector_fread(fp, w->ATb);

  fclose(fp);

  return s;
} /* lls_load() */
