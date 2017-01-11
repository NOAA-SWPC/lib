/*
 * lapack_common.h
 */

#ifndef INCLUDED_lapack_common_h
#define INCLUDED_lapack_common_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/* Prototypes */

int lapack_lls(const gsl_matrix * A, const gsl_matrix * B, gsl_matrix * X,
               int *rank, gsl_vector * rnorm);
int lapack_lls2(const gsl_matrix * A, const gsl_vector * b, gsl_vector * x,
                int *rank);
int lapack_complex_lls(const gsl_matrix_complex * A, const gsl_matrix_complex * B,
                       gsl_matrix_complex * X, int *rank);
int lapack_eigen_symm(const gsl_matrix * m, gsl_vector *eval, int *eval_found);
int lapack_eigen_symmv(const gsl_matrix * m, gsl_vector *eval, gsl_matrix *evec,
                       int *eval_found);
int lapack_eigen_herm(const gsl_matrix_complex * m, gsl_vector *eval, gsl_matrix_complex *evec,
                      int *eval_found);
int lapack_svd(const gsl_matrix * A, gsl_vector * S, gsl_matrix * U, gsl_matrix * V);
int lapack_complex_svd(const gsl_matrix_complex * A, gsl_vector * S,
                       gsl_matrix_complex * U, gsl_matrix_complex * V);

#endif /* INCLUDED_lapack_common_h */
