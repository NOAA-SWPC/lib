/*
 * io.h
 */

#ifndef INCLUDED_io_h
#define INCLUDED_io_h

#include <gsl/gsl_matrix.h>

/*
 * Prototypes
 */

int write_vector(const char *filename, const gsl_vector *v);
gsl_vector *read_vector(const char *filename);
int write_matrix(const char *filename, const gsl_matrix *m);
gsl_matrix *read_matrix(const char *filename);
int write_matrix_complex(const char *filename, const gsl_matrix_complex *m);
gsl_matrix_complex *read_matrix_complex(const char *filename);

#endif /* INCLUDED_io_h */
