/*
 * io.h
 */

#ifndef INCLUDED_io_h
#define INCLUDED_io_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>

/*
 * Prototypes
 */

int pca_write_data(const char *filename, const size_t nmax, const size_t mmax);
int pca_read_data(const char *filename, size_t *nmax, size_t *mmax);
int pca_write_vector(const char *filename, const gsl_vector *v);
gsl_vector *pca_read_vector(const char *filename);
int pca_write_matrix(const char *filename, const gsl_matrix *m);
gsl_matrix *pca_read_matrix(const char *filename);
int pca_write_matrix_complex(const char *filename, const gsl_matrix_complex *m);
gsl_matrix_complex *pca_read_matrix_complex(const char *filename);
int pca_write_S(const char *filename, const size_t nmax, const size_t mmax,
                const double freq_cpd, const double window_size,
                const double window_shift, const gsl_vector *S);
int pca_write_complex_U(const char *filename, const size_t nmax, const size_t mmax,
                        const double freq_cpd, const double window_size,
                        const double window_shift, const size_t nmodes, const gsl_matrix_complex *U);
int pca_write_complex_V(const char *filename, const size_t nmax, const size_t mmax,
                        const double freq_cpd, const double window_size,
                        const double window_shift, const size_t nmodes, const gsl_matrix_complex *V);

#endif /* INCLUDED_io_h */
