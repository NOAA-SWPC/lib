/*
 * tsqr.h
 */

#ifndef INCLUDED_tsqr_h
#define INCLUDED_tsqr_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct
{
  size_t nmax;     /* maximum number of rows to accumulate at once */
  size_t p;        /* number of columns in matrix */
  int init;

  gsl_vector *tau; /* householder scalars */
  gsl_matrix *R;   /* [ R ; A_i ], size (nmax + p)-by-p */
  gsl_vector *QTb; /* [ Q^T b ; b_i ], size (nmax + p)-by-1 */
} tsqr_workspace;

/*
 * Prototypes
 */

tsqr_workspace *tsqr_alloc(const size_t nmax, const size_t p);
void tsqr_free(tsqr_workspace *w);
int tsqr_accumulate(const gsl_matrix * A, const gsl_vector * b, tsqr_workspace * w);
int tsqr_solve(const double lambda, gsl_vector * c, tsqr_workspace * w);

#endif /* INCLUDED_tsqr_h */
