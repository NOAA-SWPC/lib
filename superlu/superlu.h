/*
 * superlu.h
 * Patrick Alken
 */

#ifndef INCLUDED_superlu_h
#define INCLUDED_superlu_h

#include <gsl/gsl_spmatrix.h>

#include <superlu_mt/slu_mt_ddefs.h>

typedef struct
{
  size_t size1;  /* number of rows in A */
  size_t size2;  /* number of columns in A */

  SuperMatrix A;
  SuperMatrix B;
  SuperMatrix X;
  SuperMatrix L;
  SuperMatrix U;
  int nrhs;
  int *perm_r;
  int *perm_c;
  int *rind;
  int *cptr;
  double *R;
  double *C;
  int nprocs; /* number of processors */

  gsl_vector *rhs_copy;

  double rnorm; /* residual norm ||A*x - b|| */
  double rcond;

  superlumt_options_t options;
} slu_workspace;

/* Prototypes */

slu_workspace *slu_alloc(int size1, int size2, int nprocs);
void slu_free(slu_workspace *w);
double slu_residual(slu_workspace *w);
int slu_proc(const gsl_spmatrix *A, const double *rhs, double *sol, slu_workspace *w);

#endif /* INCLUDED_superlu_h */
