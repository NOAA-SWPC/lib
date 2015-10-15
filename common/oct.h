/*
 * oct.h
 */

#ifndef INCLUDED_oct_h
#define INCLUDED_oct_h

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void print_octave(const gsl_matrix *m, const char *str);
void printv_octave(const gsl_vector *v, const char *str);
void printc_octave(const gsl_matrix_complex *m, const char *str);
void printcv_octave(const gsl_vector_complex *v, const char *str);
void printherm_octave(const gsl_matrix_complex *m, const char *str);

#endif /* INCLUDED_oct_h */
