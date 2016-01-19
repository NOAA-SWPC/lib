/*
 * lls.c
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
 * an r-by-ncoeff submatrix of A for a block of r observations:
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
#include <gsl/gsl_eigen.h>

#include "common.h"
#include "lls.h"
#include "oct.h"

void dgesv_(int *N, int *nrhs, double *a, int *lda, int *ipiv,
            double *b, int *ldb, int *info);
void zgesv_(int *N, int *nrhs, double *a, int *lda, int *ipiv,
            double *b, int *ldb, int *info);
void dgels_(char *trans, int *m, int *n, int *nrhs, double *a,
            int *lda, double *b, int *ldb, double *work, int *lwk,
            int *info);

static int lls_lapack_dgesv(gsl_vector *c, lls_workspace *w);
static int lls_lapack_dgels(gsl_vector *c, lls_workspace *w);
static int lls_lapack_zgesv(gsl_vector_complex *c, lls_workspace *w);

lls_workspace *
lls_alloc(size_t ncoeff)
{
  lls_workspace *w;
  int info;
  double work[1];

  w = calloc(1, sizeof(lls_workspace));
  if (!w)
    return 0;

  /* A^T A matrix */
  w->ATA = gsl_matrix_alloc(ncoeff, ncoeff);
  w->work_A = gsl_matrix_alloc(ncoeff, ncoeff);
  w->AHA = gsl_matrix_complex_alloc(ncoeff, ncoeff);
  w->work_A_complex = gsl_matrix_complex_alloc(ncoeff, ncoeff);
  if (!w->ATA || !w->work_A || !w->AHA || !w->work_A_complex)
    {
      lls_free(w);
      return 0;
    }

  w->ATb = gsl_vector_alloc(ncoeff);
  w->work_b = gsl_vector_alloc(ncoeff);
  w->AHb = gsl_vector_complex_alloc(ncoeff);
  w->work_b_complex = gsl_vector_complex_alloc(ncoeff);
  w->lapack_ipiv = malloc(ncoeff * sizeof(int));
  if (!w->ATb || !w->AHb || !w->work_b || !w->lapack_ipiv)
    {
      lls_free(w);
      return 0;
    }

  w->lapack_nrhs = 1;
  w->lapack_trans = 'N';
  w->lapack_lwk = -1;
  dgels_(&(w->lapack_trans),
         (int *) &ncoeff,
         (int *) &ncoeff,
         &(w->lapack_nrhs),
         w->work_A->data,
         (int *) &ncoeff,
         w->work_b->data,
         (int *) &ncoeff,
         work,
         &(w->lapack_lwk),
         &info);

  if (info == 0)
    {
      w->lapack_lwk = (int) work[0];
      w->lapack_work = malloc(w->lapack_lwk * sizeof(double));
    }
  else
    {
      fprintf(stderr, "lls_alloc: dgels query failed: info = %d\n", info);
      lls_free(w);
      return 0;
    }

  w->ncoeff = ncoeff;

  w->residual = 0.0;
  w->cond = -1.0;
  w->eval_excluded = -1;

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

  if (w->AHA)
    gsl_matrix_complex_free(w->AHA);

  if (w->AHb)
    gsl_vector_complex_free(w->AHb);

  if (w->work_A_complex)
    gsl_matrix_complex_free(w->work_A_complex);

  if (w->work_b_complex)
    gsl_vector_complex_free(w->work_b_complex);

  if (w->lapack_ipiv)
    free(w->lapack_ipiv);

  if (w->lapack_work)
    free(w->lapack_work);

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
  gsl_matrix_complex_set_zero(w->AHA);
  gsl_vector_complex_set_zero(w->AHb);

  return 0;
}

/*
lls_calc_ATA()
  Add a set of nobs observations to the current A^T A matrix. This
routine is designed so it can be called multiple times, so the caller
can read in a set of observations, add them to the design matrix, and
then repeat

Inputs: A - matrix of basis functions (size nobs-by-ncoeff)
            A_{ij} = g_j(x_i)
        w - workspace

Notes: ATA <- ATA + A^T A
*/

int
lls_calc_ATA(const gsl_matrix *A, lls_workspace *w)
{
  if (A->size2 != w->ncoeff)
    {
      fprintf(stderr, "lls_calc_ATA: A has wrong size2\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;
      
      s = gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, A, 1.0, w->ATA);

      return s;
    }
} /* lls_calc_ATA() */

/*
lls_calc_ATb()
  Add a set of nobs observations to the current A^T b right hand side. This
routine is designed so it can be called multiple times, so the caller
can read in a set of observations, add them to the RHS vector, and
then repeat

Inputs: A - matrix of basis functions (size nobs-by-ncoeff)
            A_{ij} = g_j(x_i)
        b - right hand side vector
        w - workspace

Notes: ATb <- ATb + A^T b
*/

int
lls_calc_ATb(const gsl_matrix *A, const gsl_vector *b, lls_workspace *w)
{
  if (A->size2 != w->ncoeff)
    {
      fprintf(stderr, "lls_calc_ATb: A has wrong size2\n");
      return GSL_EBADLEN;
    }
  else if (A->size1 != b->size)
    {
      fprintf(stderr, "lls_calc_ATb: size1 of A does not match b\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;
      
      s = gsl_blas_dgemv(CblasTrans, 1.0, A, b, 1.0, w->ATb);

      return s;
    }
} /* lls_calc_ATb() */

/*
lls_apply_W()
  Apply weight matrix W = diag(wts) to least squares system

Inputs: wts - (input) weight vector
        b   - (output) rhs vector: b <- W b
        A   - (output) LS matrix: A <- W A

Return: success or error
*/

int
lls_apply_W(const gsl_vector *wts, gsl_vector *b, gsl_matrix *A)
{
  if (A->size1 != b->size)
    {
      fprintf(stderr, "lls_apply_W: size1 of A does not match b\n");
      return GSL_EBADLEN;
    }
  else if (b->size != wts->size)
    {
      fprintf(stderr, "lls_apply_W: size of b does not match wts\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;
      size_t i;

      for (i = 0; i < A->size1; ++i)
        {
          gsl_vector_view rv = gsl_matrix_row(A, i);
          double wi = gsl_vector_get(wts, i);
          double bi = gsl_vector_get(b, i);

          /* b <- W b */
          gsl_vector_set(b, i, wi * bi);

          /* A <- W A */
          gsl_vector_scale(&rv.vector, wi);
        }

      return s;
    }
} /* lls_apply_W() */

/*
lls_calc_ATWb()
  Add a set of nobs observations to the current A^T W b right hand side. This
routine is designed so it can be called multiple times, so the caller
can read in a set of observations, add them to the RHS vector, and
then repeat

Inputs: A   - matrix of basis functions (size nobs-by-ncoeff)
              A_{ij} = g_j(x_i)
        b   - right hand side vector
        wts - weight vector = diag(W)
        w   - workspace

Notes: ATb <- ATb + A^T W b
*/

int
lls_calc_ATWb(const gsl_matrix *A, gsl_vector *b, const gsl_vector *wts, lls_workspace *w)
{
  if (A->size2 != w->ncoeff)
    {
      fprintf(stderr, "lls_calc_ATWb: A has wrong size2\n");
      return GSL_EBADLEN;
    }
  else if (A->size1 != b->size)
    {
      fprintf(stderr, "lls_calc_ATWb: size1 of A does not match b\n");
      return GSL_EBADLEN;
    }
  else if (b->size != wts->size)
    {
      fprintf(stderr, "lls_calc_ATWb: size of b does not match wts\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;
 
      /* b <- W b */
      gsl_vector_mul(b, wts);

      /* ATb <- ATb + A^T W b */
      s = gsl_blas_dgemv(CblasTrans, 1.0, A, b, 1.0, w->ATb);

      return s;
    }
} /* lls_calc_ATWb() */

/*
lls_calc_AHA()
  Add a set of nobs observations to the current A^H A matrix. This
routine is designed so it can be called multiple times, so the caller
can read in a set of observations, add them to the design matrix, and
then repeat

Inputs: A - matrix of basis functions (size nobs-by-ncoeff)
            A_{ij} = g_j(x_i)
        w - workspace

Notes: AHA <- AHA + A^H A
*/

int
lls_calc_AHA(gsl_matrix_complex *A, lls_workspace *w)
{
  if (A->size2 != w->ncoeff)
    {
      fprintf(stderr, "lls_calc_AHA: A has wrong size2\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;
      gsl_complex one;

      GSL_SET_COMPLEX(&one, 1.0, 0.0);
 
      /* AHA += A^H A */
      s = gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, one, A, A, one, w->AHA);

      return s;
    }
} /* lls_calc_AHA() */

/*
lls_calc_AHb()
  Add a set of nobs observations to the current A^H b right hand side. This
routine is designed so it can be called multiple times, so the caller
can read in a set of observations, add them to the RHS vector, and
then repeat

Inputs: A - matrix of basis functions (size nobs-by-ncoeff)
            A_{ij} = g_j(x_i)
        b - right hand side vector
        w - workspace

Notes: AHb <- AHb + A^H b
*/

int
lls_calc_AHb(gsl_matrix_complex *A, gsl_vector_complex *b, lls_workspace *w)
{
  if (A->size2 != w->ncoeff)
    {
      fprintf(stderr, "lls_calc_AHb: A has wrong size2\n");
      return GSL_EBADLEN;
    }
  else if (A->size1 != b->size)
    {
      fprintf(stderr, "lls_calc_AHb: size1 of A does not match b\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;
      gsl_complex one;

      GSL_SET_COMPLEX(&one, 1.0, 0.0);
      
      s = gsl_blas_zgemv(CblasConjTrans, one, A, b, one, w->AHb);

      return s;
    }
} /* lls_calc_AHb() */

/*
lls_minmax()
  Return the minimum and maximum elements of A^T A and A^T b
*/

int
lls_minmax(double *min_ATA, double *max_ATA, double *min_rhs,
           double *max_rhs, lls_workspace *w)
{
  int s = 0;
  
  gsl_matrix_minmax(w->ATA, min_ATA, max_ATA);
  gsl_vector_minmax(w->ATb, min_rhs, max_rhs);

  return s;
} /* lls_maxmin() */

/*
lls_regularize()
  Perform Tikhonov regularization on the least-squares
problem by adding a positive constant to the diagonal
matrix terms

x = (A^T A + mu*I)^{-1) A^T b

Inputs: mu - positive constant
        w  - workspace
*/

int
lls_regularize(const double mu, lls_workspace *w)
{
  gsl_vector_view d = gsl_matrix_diagonal(w->ATA);
  gsl_vector_add_constant(&d.vector, mu);

  return GSL_SUCCESS;
} /* lls_regularize() */

/*
lls_solve()
  Solve the least squares system:

A^T A c = A^T b

where A^T A and A^T b have been previously computed with
lls_calc_ATA and lls_calc_ATb

Inputs: c - (output) coefficient vector
        w - workspace

Notes: on output, the residual || A^T A c - A^T b || is stored in
w->residual
*/

int
lls_solve(gsl_vector *c, lls_workspace *w)
{
  if (c->size != w->ncoeff)
    {
      fprintf(stderr, "lls_solve: coefficient vector has wrong size\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;

      if (w->ATA->size1 < 1000)
        {
          print_octave(w->ATA, "ATA");
          printv_octave(w->ATb, "ATb");
        }

#if 1
      s = lls_lapack_dgesv(c, w); /* matrix solve */
#else
      s = lls_lapack_dgels(c, w); /* least squares solve */
#endif

      /* compute residual || ATA c - ATb || */
      gsl_vector_memcpy(w->work_b, w->ATb);
      gsl_blas_dgemv(CblasNoTrans, 1.0, w->ATA, c, -1.0, w->work_b);
      w->residual = gsl_blas_dnrm2(w->work_b);

      printv_octave(c, "c");

      return s;
    }
} /* lls_solve() */

/*
lls_solve_complex()
  Solve the least squares system:

A^H A c = A^H b

where A^H A and A^H b have been previously computed with
lls_calc_AHA and lls_calc_AHb

Inputs: c - (output) coefficient vector
        w - workspace

Notes: on output, the residual || A^H A c - A^H b || is stored in
w->residual
*/

int
lls_solve_complex(gsl_vector_complex *c, lls_workspace *w)
{
  if (c->size != w->ncoeff)
    {
      fprintf(stderr, "lls_solve_complex: coefficient vector has wrong size\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;
      gsl_complex alpha, beta;

      if (w->AHA->size1 < 1000)
        {
          printc_octave(w->AHA, "AHA");
          printcv_octave(w->AHb, "AHb");
        }

      s = lls_lapack_zgesv(c, w); /* matrix solve */

      /* compute residual || AHA c - AHb || */
      GSL_SET_COMPLEX(&alpha, 1.0, 0.0);
      GSL_SET_COMPLEX(&beta, -1.0, 0.0);

      gsl_vector_complex_memcpy(w->work_b_complex, w->AHb);
      gsl_blas_zgemv(CblasNoTrans, alpha, w->AHA, c, beta,
                     w->work_b_complex);
      w->residual = gsl_blas_dznrm2(w->work_b_complex);

      printcv_octave(c, "c");

      return s;
    }
} /* lls_solve_complex() */

/*
lls_solve_spectral()
  Solve the least squares system:

A^T A c = A^T b

where A^T A and A^T b have been previously computed with
lls_calc_ATA and lls_calc_ATb. The method used is eigenvalue
decomposition of the A^T A matrix.

Inputs: c - (output) coefficient vector
        w - workspace

Notes:

1) on output, the residual || A^T A c - A^T b || is stored in
w->residual

2) on output, the matrix condition number is stored in w->cond
*/

int
lls_solve_spectral(gsl_vector *c, lls_workspace *w)
{
  const size_t N = w->ncoeff;

  if (c->size != N)
    {
      fprintf(stderr, "lls_solve_spectral: coefficient vector has wrong size\n");
      return GSL_EBADLEN;
    }
  else
    {
      int s = 0;
      size_t i;
      gsl_eigen_symmv_workspace *eigen_p = gsl_eigen_symmv_alloc(N);
      gsl_vector *lambda = gsl_vector_alloc(N);
      gsl_matrix *evec = gsl_matrix_alloc(N, N);
      double lambda_max, /* largest eigenvalue */
             lambda_min; /* smallest eigenvalue */
 
      if (w->ATA->size1 < 1000)
        {
          print_octave(w->ATA, "ATA");
          printv_octave(w->ATb, "ATb");
        }

      gsl_matrix_memcpy(w->work_A, w->ATA);
      s = gsl_eigen_symmv(w->work_A, lambda, evec, eigen_p);
      if (s)
        return s;

      /* sort eigenvalues in descending order */
      gsl_eigen_symmv_sort(lambda, evec, GSL_EIGEN_SORT_VAL_DESC);

      lambda_max = gsl_vector_get(lambda, 0);
      lambda_min = gsl_vector_get(lambda, N - 1);

      /* initialize coefficient vector to 0 */
      gsl_vector_set_zero(c);

      for (i = 0; i < N; ++i)
        {
          double li = gsl_vector_get(lambda, i);
          gsl_vector_view ei = gsl_matrix_column(evec, i);
          double bi;

#if 0
          /*
           * Once we get to eigenvalues sufficiently small, stop adding
           * them to the sum
           */
          if (li < lambda_max * 1.0e-3)
            break;
#endif

          /* compute bi = (A^T b) . ei */
          gsl_blas_ddot(w->ATb, &ei.vector, &bi);

          /* c += (b_i / lambda_i) * e_i */
          gsl_blas_daxpy(bi / li, &ei.vector, c);
        }

      /* count number of eigenvalues excluded from coefficient sum */
      w->eval_excluded = N - i;

      /* compute residual || ATA c - ATb || */
      gsl_vector_memcpy(w->work_b, w->ATb);
      gsl_blas_dgemv(CblasNoTrans, 1.0, w->ATA, c, -1.0, w->work_b);
      w->residual = gsl_blas_dnrm2(w->work_b);

      /* compute condition number of A^T A */
      w->cond = fabs(lambda_max / lambda_min);

      printv_octave(c, "c");
      printv_octave(lambda, "lambda");

      gsl_eigen_symmv_free(eigen_p);
      gsl_vector_free(lambda);
      gsl_matrix_free(evec);

      return s;
    }
} /* lls_solve_spectral() */

/*
lls_save()
  Save A^T A matrix and A^T b right hand side to a binary file
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

  s += gsl_matrix_fwrite(fp, w->ATA);
  s += gsl_vector_fwrite(fp, w->ATb);
  s += gsl_matrix_complex_fwrite(fp, w->AHA);
  s += gsl_vector_complex_fwrite(fp, w->AHb);

  fclose(fp);

  return s;
} /* lls_save() */

/*
lls_load()
  Load A^T A matrix and A^T b right hand side from a binary file
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

  s += gsl_matrix_fread(fp, w->ATA);
  s += gsl_vector_fread(fp, w->ATb);
  s += gsl_matrix_complex_fread(fp, w->AHA);
  s += gsl_vector_complex_fread(fp, w->AHb);

  fclose(fp);

  return s;
} /* lls_load() */

static int
lls_lapack_dgesv(gsl_vector *c, lls_workspace *w)
{
  int s = 0;
  int N = w->work_A->size1;
  int nrhs = 1;
  int lda = w->work_A->size1;
  int ldb = w->work_b->size;

  gsl_matrix_memcpy(w->work_A, w->ATA);
  gsl_vector_memcpy(w->work_b, w->ATb);

  dgesv_(&N,
         &nrhs,
         w->work_A->data,
         &lda,
         w->lapack_ipiv,
         w->work_b->data,
         &ldb,
         &s);

  /* save solution vector */
  gsl_vector_memcpy(c, w->work_b);

  return s;
} /* lls_lapack_dgesv() */

static int
lls_lapack_dgels(gsl_vector *c, lls_workspace *w)
{
  int s = 0;
  int N = w->work_A->size1;
  int nrhs = 1;
  int lda = w->work_A->size1;
  int ldb = w->work_b->size;
  char trans = 'N';

  gsl_matrix_memcpy(w->work_A, w->ATA);
  gsl_vector_memcpy(w->work_b, w->ATb);

  dgels_(&trans,
         &N,
         &N,
         &nrhs,
         w->work_A->data,
         &lda,
         w->work_b->data,
         &ldb,
         w->lapack_work,
         &(w->lapack_lwk),
         &s);

  /* save solution vector */
  gsl_vector_memcpy(c, w->work_b);

  return s;
} /* lls_lapack_dgels() */

static int
lls_lapack_zgesv(gsl_vector_complex *c, lls_workspace *w)
{
  int s = 0;
  int N = w->work_A_complex->size1;
  int nrhs = 1;
  int lda = w->work_A_complex->size1;
  int ldb = w->work_b_complex->size;

  gsl_matrix_complex_transpose_memcpy(w->work_A_complex, w->AHA);
  gsl_vector_complex_memcpy(w->work_b_complex, w->AHb);

  zgesv_(&N,
         &nrhs,
         w->work_A_complex->data,
         &lda,
         w->lapack_ipiv,
         w->work_b_complex->data,
         &ldb,
         &s);

  /* save solution vector */
  gsl_vector_complex_memcpy(c, w->work_b_complex);

  return s;
} /* lls_lapack_zgesv() */
