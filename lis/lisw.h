/*
 * lisw.h
 * Patrick Alken
 */

#ifndef INCLUDED_lisw_h
#define INCLUDED_lisw_h

#include <gsl/gsl_spmatrix.h>

#include "lis.h"

typedef struct
{
  int size1;     /* number of rows */
  int size2;     /* number of columns */
  double rnorm;  /* || b - A*x || */
  double rrnorm; /* || b - A*x || / || b || */
} lis_workspace;

/*
 * Prototypes
 */

lis_workspace *lis_alloc(int size1, int size2);
void mylis_free(lis_workspace *w);
int lis_proc(const gsl_spmatrix *A, const double *rhs, const double tol,
             double *sol, lis_workspace *w);

#endif /* INCLUDED_lisw_h */
