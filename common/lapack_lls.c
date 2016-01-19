/* lapack_lls.c
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * This module is a wrapper for the LAPACK linear least squares problem
 * (LLS).
 */

#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "lapack_lls.h"

void
dgelsd_(int *m, int *n, int *nrhs, double *a, int *lda, double *b,
        int *ldb, double *s, double *rcond, int *rank, double *work,
        int *lwork, int *iwork, int *info);

/*
lls()
  This routine solves

min_x || A x - c ||_2
*/

int
lls(const gsl_matrix *A, const gsl_vector *c, gsl_vector *x)
{
  int m = (int) A->size1;
  int n = (int) A->size2;
  int nrhs = 1;
  int info;
  int lwork;
  gsl_matrix *aa, *bb;
  gsl_vector *s;
  gsl_vector *work;
  double q[1];
  int ldb = GSL_MAX(m, n);
  int lda = m;
  double rcond = 1.0e-12;
  int rank;
  int *iwork = 0;
  gsl_vector_view v;
  gsl_vector *rhs;

  rhs = gsl_vector_alloc(c->size);
  aa = gsl_matrix_alloc(A->size2, A->size1);
  bb = gsl_matrix_alloc(nrhs, GSL_MAX(m, n));
  s = gsl_vector_alloc(GSL_MIN(m, n));

  gsl_matrix_transpose_memcpy(aa, A);
  gsl_vector_memcpy(rhs, c);

  v = gsl_matrix_subrow(bb, 0, 0, m);
  gsl_vector_memcpy(&v.vector, rhs);

  lwork = -1;
  dgelsd_(&m,
          &n,
          &nrhs,
          aa->data,
          &lda,
          bb->data,
          &ldb,
          s->data,
          &rcond,
          &rank,
          q,
          &lwork,
          iwork,
          &info);

  lwork = (int) q[0];
  work = gsl_vector_alloc((size_t) lwork);
  iwork = malloc(sizeof(int) * m);

  dgelsd_(&m,
          &n,
          &nrhs,
          aa->data,
          &lda,
          bb->data,
          &ldb,
          s->data,
          &rcond,
          &rank,
          work->data,
          &lwork,
          iwork,
          &info);

  v = gsl_matrix_subrow(bb, 0, 0, n);
  gsl_vector_memcpy(x, &v.vector);

  gsl_matrix_free(aa);
  gsl_matrix_free(bb);
  gsl_vector_free(s);
  gsl_vector_free(rhs);
  gsl_vector_free(work);
  free(iwork);

  if (info)
    fprintf(stderr, "ERROR: lls: info = %d\n", info);

  return (info);
} /* lls() */

/*
wlls()
  This routine solves the weighted least squares problem

min_x || A x - c ||_2

with weights w
*/

int
wlls(const gsl_matrix *A, const gsl_vector *w, const gsl_vector *c,
     gsl_vector *x)
{
  size_t i;
  int m = (int) A->size1;
  int n = (int) A->size2;
  int nrhs = 1;
  int info;
  int lwork;
  gsl_matrix *aa, *bb;
  gsl_vector *s;
  gsl_vector *work;
  double q[1];
  int ldb = GSL_MAX(m, n);
  int lda = m;
  double rcond = 1.0e-12;
  int rank;
  int *iwork = 0;
  gsl_vector_view v;
  gsl_vector *rhs;

  rhs = gsl_vector_alloc(c->size);
  aa = gsl_matrix_alloc(A->size2, A->size1);

  gsl_matrix_transpose_memcpy(aa, A);

  for (i = 0; i < A->size1; ++i)
    {
      double wi = gsl_vector_get(w, i);

      v = gsl_matrix_column(aa, i);
      gsl_vector_scale(&v.vector, sqrt(wi));

      gsl_vector_set(rhs, i, gsl_vector_get(c, i) * sqrt(wi));
    }

  bb = gsl_matrix_alloc(nrhs, GSL_MAX(m, n));
  s = gsl_vector_alloc(GSL_MIN(m, n));

  v = gsl_matrix_subrow(bb, 0, 0, m);
  gsl_vector_memcpy(&v.vector, rhs);

  lwork = -1;
  dgelsd_(&m,
          &n,
          &nrhs,
          aa->data,
          &lda,
          bb->data,
          &ldb,
          s->data,
          &rcond,
          &rank,
          q,
          &lwork,
          iwork,
          &info);

  lwork = (int) q[0];
  work = gsl_vector_alloc((size_t) lwork);
  iwork = malloc(sizeof(int) * m);

  dgelsd_(&m,
          &n,
          &nrhs,
          aa->data,
          &lda,
          bb->data,
          &ldb,
          s->data,
          &rcond,
          &rank,
          work->data,
          &lwork,
          iwork,
          &info);

  v = gsl_matrix_subrow(bb, 0, 0, n);
  gsl_vector_memcpy(x, &v.vector);

  gsl_matrix_free(aa);
  gsl_matrix_free(bb);
  gsl_vector_free(s);
  gsl_vector_free(rhs);
  gsl_vector_free(work);
  free(iwork);

  if (info)
    fprintf(stderr, "ERROR: wlls: info = %d\n", info);

  return (info);
} /* wlls() */
