/*
 * superlu_mt.c
 * Patrick Alken
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <superlu_mt/slu_mt_ddefs.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>

#include "superlu.h"

/*
slu_alloc()
  Allocate a SuperLU workspace

Inputs: size1  - number of rows in sparse matrix
        size2  - number of columns in sparse matrix
        nprocs - number of processors to use

Return: pointer to workspace
*/

slu_workspace *
slu_alloc(int size1, int size2, int nprocs)
{
  slu_workspace *w;

  w = calloc(1, sizeof(slu_workspace));
  if (!w)
    {
      fprintf(stderr, "slu_alloc failed\n");
      return 0;
    }

  w->perm_r = malloc(size1 * sizeof(int));
  w->perm_c = malloc(size2 * sizeof(int));

  if (w->perm_r == 0 || w->perm_c == 0)
    {
      slu_free(w);
      return 0;
    }

  w->rhs_copy = gsl_vector_alloc(size1);
  w->cptr = malloc((size2 + 1) * sizeof(int));

  w->R = malloc(size1 * sizeof(double));
  w->C = malloc(size2 * sizeof(double));

  w->size1 = size1;
  w->size2 = size2;

  w->nprocs = nprocs;
  w->nrhs = 1;
  w->residual = 0.0;
  w->rcond = 0.0;

  w->options.nprocs = nprocs;
  w->options.fact = EQUILIBRATE;
  w->options.trans = NOTRANS;
  w->options.refact = NO;
  w->options.panel_size = sp_ienv(1);
  w->options.relax = sp_ienv(2);
  w->options.usepr = NO;
  w->options.drop_tol = 0.0;
  w->options.diag_pivot_thresh = 1.0;
  w->options.SymmetricMode = NO;
  w->options.PrintStat = NO;
  w->options.perm_c = w->perm_c;
  w->options.perm_r = w->perm_r;
  w->options.work = NULL;
  w->options.lwork = 0;

  w->options.etree = malloc(size2 * sizeof(int));
  w->options.colcnt_h = malloc(size2 * sizeof(int));
  w->options.part_super_h = malloc(size2 * sizeof(int));

  return (w);
} /* slu_alloc() */

void
slu_free(slu_workspace *w)
{
  if (w->perm_r)
    free(w->perm_r);

  if (w->perm_c)
    free(w->perm_c);

  if (w->rhs_copy)
    gsl_vector_free(w->rhs_copy);

  if (w->cptr)
    free(w->cptr);

  if (w->R)
    free(w->R);

  if (w->C)
    free(w->C);

  if (w->options.etree)
    free(w->options.etree);

  if (w->options.colcnt_h)
    free(w->options.colcnt_h);

  if (w->options.part_super_h)
    free(w->options.part_super_h);

  free(w);
} /* slu_free() */

double
slu_residual(slu_workspace *w)
{
  return w->residual;
}

#if 0

/*
slu_proc()
  Solve a sparse linear system A x = b

Inputs: A      - sparse matrix in CCS format
        rhs    - rhs vector (length w->size1)
        sol    - (output) solution vector (length w->size2)
        w      - workspace

Return: 0 on success, non-zero on error
*/

int
slu_proc(const gsl_spmatrix *A, const double *rhs, double *sol, slu_workspace *w)
{
  const size_t nnz = gsl_spmatrix_nnz(A);
  int info = 0;
  gsl_vector_const_view vrhs = gsl_vector_const_view_array(rhs, w->size1);
  gsl_vector_view vsol = gsl_vector_view_array(sol, w->size2);
  int *rind = malloc(nnz * sizeof(int));
  size_t i;

  /* make copy of rhs vector */
  gsl_vector_memcpy(w->rhs_copy, &vrhs.vector);

  /* have to copy arrays since sizeof(int) != sizeof(size_t) */
  for (i = 0; i < nnz; ++i)
    rind[i] = (int) A->i[i];
  for (i = 0; i < w->size2 + 1; ++i)
    w->cptr[i] = (int) A->p[i];

  dCreate_CompCol_Matrix(&(w->A),
                         (int) A->size1,
                         (int) A->size2,
                         nnz,
                         A->data,
                         rind,
                         w->cptr,
                         SLU_NC,
                         SLU_D,
                         SLU_GE);

  /* rhs matrix */
  dCreate_Dense_Matrix(&(w->B),
                       w->size1,
                       w->nrhs,
                       w->rhs_copy->data,
                       w->size1,
                       SLU_DN,
                       SLU_D,
                       SLU_GE);

  get_perm_c(1, &(w->A), w->perm_c);

  pdgssv(w->nprocs,
         &(w->A),
         w->perm_c,
         w->perm_r,
         &(w->L),
         &(w->U),
         &(w->B),
         &info);
  if (info != 0)
    {
      fprintf(stderr, "slu_proc: error in pdgssv: info = %d\n", info);
      return info; /* error */
    }

  /* now store solution in vector 'sol' */
  {
    size_t i;
    DNformat *Astore = (DNformat *) w->B.Store;
    double *dp = (double *) Astore->nzval;

    for (i = 0; i < A->size2; ++i)
      sol[i] = dp[i];
  }

  /* compute residual */
  gsl_vector_memcpy(w->rhs_copy, &vrhs.vector);
  gsl_spblas_dgemv(CblasNoTrans, 1.0, A, &vsol.vector, -1.0, w->rhs_copy);
  w->residual = gsl_blas_dnrm2(w->rhs_copy);

  Destroy_SuperMatrix_Store(&(w->A));
  Destroy_SuperMatrix_Store(&(w->B));
  Destroy_SuperNode_SCP(&(w->L));
  Destroy_CompCol_NCP(&(w->U));

  free(rind);

  return info;
} /* slu_proc() */

#endif

/*
slu_proc()
  Solve a sparse linear system A x = b

Inputs: A      - sparse matrix in CCS format
        rhs    - rhs vector (length w->size1)
        sol    - (output) solution vector (length w->size2)
        w      - workspace

Return: 0 on success, non-zero on error
*/

int
slu_proc(const gsl_spmatrix *A, const double *rhs, double *sol, slu_workspace *w)
{
  const size_t nnz = gsl_spmatrix_nnz(A);
  int info = 0;
  gsl_vector_const_view vrhs = gsl_vector_const_view_array(rhs, w->size1);
  gsl_vector_view vsol = gsl_vector_view_array(sol, w->size2);
  int *rind = malloc(nnz * sizeof(int));
  double *data = malloc(nnz * sizeof(double));
  size_t i;

  /* make copy of rhs vector */
  gsl_vector_memcpy(w->rhs_copy, &vrhs.vector);

  /* have to copy arrays since sizeof(int) != sizeof(size_t) */
  for (i = 0; i < nnz; ++i)
    {
      rind[i] = (int) A->i[i];
      data[i] = A->data[i];
    }

  for (i = 0; i < w->size2 + 1; ++i)
    w->cptr[i] = (int) A->p[i];

  dCreate_CompCol_Matrix(&(w->A),
                         (int) A->size1,
                         (int) A->size2,
                         nnz,
                         data,
                         rind,
                         w->cptr,
                         SLU_NC,
                         SLU_D,
                         SLU_GE);

  /* rhs matrix */
  dCreate_Dense_Matrix(&(w->B),
                       w->size1,
                       w->nrhs,
                       w->rhs_copy->data,
                       w->size1,
                       SLU_DN,
                       SLU_D,
                       SLU_GE);

  /* solution matrix */
  dCreate_Dense_Matrix(&(w->X),
                       w->size1,
                       w->nrhs,
                       sol,
                       w->size1,
                       SLU_DN,
                       SLU_D,
                       SLU_GE);

  get_perm_c(1, &(w->A), w->perm_c);

  {
    equed_t equed = NOEQUIL;
    double ferr, berr, rpg;
    superlu_memusage_t memusage;

    pdgssvx(w->nprocs, &(w->options), &(w->A), w->perm_c, w->perm_r,
            &equed, w->R, w->C, &(w->L), &(w->U), &(w->B), &(w->X), &rpg, &(w->rcond),
	          &ferr, &berr, &memusage, &info);
  }

  if (info != 0)
    {
      fprintf(stderr, "slu_proc: error in pdgssv: info = %d\n", info);
      return info; /* error */
    }

  /* compute residual */
  gsl_vector_memcpy(w->rhs_copy, &vrhs.vector);
  gsl_spblas_dgemv(CblasNoTrans, 1.0, A, &vsol.vector, -1.0, w->rhs_copy);
  w->residual = gsl_blas_dnrm2(w->rhs_copy);

  Destroy_SuperMatrix_Store(&(w->A));
  Destroy_SuperMatrix_Store(&(w->B));
  Destroy_SuperMatrix_Store(&(w->X));
  Destroy_SuperNode_SCP(&(w->L));
  Destroy_CompCol_NCP(&(w->U));

  free(rind);
  free(data);

  return GSL_SUCCESS;
} /* slu_proc() */
