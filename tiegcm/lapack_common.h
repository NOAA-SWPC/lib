/*
 * lapack_common.h
 */

#ifndef INCLUDED_lapack_common_h
#define INCLUDED_lapack_common_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

/* Prototypes */

int lapack_lls(const gsl_matrix * A, const gsl_matrix * B, gsl_matrix * X,
               int *rank);
int lapack_complex_lls(const gsl_matrix_complex * A, const gsl_matrix_complex * B,
                       gsl_matrix_complex * X, int *rank);
int lapack_eigen_herm(const gsl_matrix_complex * m, gsl_vector *eval, gsl_matrix_complex *evec,
                      int *eval_found);

#endif /* INCLUDED_lapack_common_h */
