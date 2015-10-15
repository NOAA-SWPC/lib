#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

int fit2_proc(gsl_matrix *A, gsl_vector *b, gsl_vector *x);
int fit2_wproc(const gsl_matrix *X, const gsl_vector *y,
               const gsl_vector *w, gsl_vector *c);
