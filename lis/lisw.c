/*
 * lisw.c
 *
 * Wrapper functions for lis library
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>

#include "lisw.h"

lis_workspace *
lis_alloc(int size1, int size2)
{
  lis_workspace *w;

  w = calloc(1, sizeof(lis_workspace));
  if (!w)
    return 0;

  w->size1 = size1;
  w->size2 = size2;

  return w;
} /* lis_alloc() */

void
mylis_free(lis_workspace *w)
{
  free(w);
} /* lis_free() */

/*
lis_proc()
  Solve the system A x = b through QR reduction

Inputs: S   - sparse matrix in triplet format
        rhs - right hand side vector b
        tol - relative tolerance in solution
        sol - (output) where to store solution x
        w   - workspace
*/

int
lis_proc(const gsl_spmatrix *S, const double *rhs, const double tol,
         double *sol, lis_workspace *w)
{
  int s = 0;
  size_t i;
  LIS_MATRIX A;
  LIS_VECTOR b, x;
  LIS_SOLVER solver;
  LIS_INT size1 = w->size1;
  LIS_INT size2 = w->size2;
  LIS_REAL rrnorm; /* || b - Ax || / ||b|| */
  int argc = 0;
  char **argv = NULL;
  char str[2048];

  lis_initialize(&argc, &argv);

  lis_matrix_create(0, &A);
  lis_matrix_set_size(A, 0, size1);

  lis_vector_create(0, &b);
  lis_vector_create(0, &x);
  lis_vector_set_size(b, 0, size1);
  lis_vector_set_size(x, 0, size2);

  lis_solver_create(&solver);

  /* set solver parameters */
  lis_solver_set_option("-i fgmres -p ilut -f double", solver);
  lis_solver_set_option("-maxiter 2000", solver);
  lis_solver_set_option("-print 1", solver);
  sprintf(str, "-tol %e\n", tol);
  lis_solver_set_option(str, solver);
  /*lis_solver_set_optionC(solver);*/

  /* construct LIS_MATRIX type from A */
  for (i = 0; i < S->nz; ++i)
    lis_matrix_set_value(LIS_INS_VALUE, S->i[i], S->p[i], S->data[i], A);

  /* construct RHS */
  for (i = 0; i < w->size1; ++i)
    lis_vector_set_value(LIS_INS_VALUE, i, rhs[i], b);

  lis_matrix_set_type(A, LIS_MATRIX_CSR);
  lis_matrix_assemble(A);

  s = lis_solve(A, b, x, solver);
  s = solver->retcode; /*XXX bug in lis_solve */

  lis_solver_get_status(solver, &s);
  if (s != 0)
    fprintf(stderr, "lis_proc: error: status = %d\n", s);

  for (i = 0; i < w->size1; ++i)
    lis_vector_get_value(x, i, &sol[i]);

  /* compute residual norm */
  lis_solver_get_residualnorm(solver, &rrnorm);
  w->rrnorm = rrnorm;

  {
    gsl_vector_const_view bv = gsl_vector_const_view_array(rhs, w->size1);
    gsl_vector_view xv = gsl_vector_view_array(sol, w->size2);
    gsl_vector *r = gsl_vector_alloc(w->size1);

    gsl_vector_memcpy(r, &bv.vector);
    gsl_spblas_dgemv(CblasNoTrans, 1.0, S, &xv.vector, -1.0, r);
    w->rnorm = gsl_blas_dnrm2(r);

    gsl_vector_free(r);
  }

  lis_matrix_destroy(A);
  lis_vector_destroy(b);
  lis_vector_destroy(x);
  lis_solver_destroy(solver);

  lis_finalize();

  return s;
} /* lis_proc() */
