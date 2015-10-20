/*
 * lls.h
 * Patrick Alken
 */

#ifndef INCLUDED_lls_h
#define INCLUDED_lls_h

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>

typedef struct
{
  size_t ncoeff;
  gsl_matrix *ATA; /* A^T A */
  gsl_vector *ATb; /* A^T b */
  double residual;
  double cond;     /* condition number */
  int eval_excluded; /* number of eigenvalues excluded */

  gsl_matrix_complex *AHA; /* A^H A */
  gsl_vector_complex *AHb; /* A^H b */

  gsl_matrix *work_A;
  gsl_vector *work_b;
  int lapack_nrhs;

  gsl_matrix_complex *work_A_complex;
  gsl_vector_complex *work_b_complex;

  /* dgesv() */
  int *lapack_ipiv;

  /* dgels() */
  char lapack_trans;
  double *lapack_work;
  int lapack_lwk;
} lls_workspace;

/*
 * Prototypes
 */

lls_workspace *lls_alloc(size_t ncoeff);
void lls_free(lls_workspace *w);
int lls_reset(lls_workspace *w);
int lls_calc_ATA(const gsl_matrix *A, lls_workspace *w);
int lls_calc_ATb(const gsl_matrix *A, const gsl_vector *b, lls_workspace *w);
int lls_apply_W(const gsl_vector *wts, gsl_vector *b, gsl_matrix *A);
int lls_calc_AHA(gsl_matrix_complex *A, lls_workspace *w);
int lls_calc_AHb(gsl_matrix_complex *A, gsl_vector_complex *b,
                 lls_workspace *w);
int lls_minmax(double *min_ATA, double *max_ATA, double *min_rhs,
               double *max_rhs, lls_workspace *w);
int lls_regularize(const double mu, lls_workspace *w);
int lls_solve(gsl_vector *c, lls_workspace *w);
int lls_solve_complex(gsl_vector_complex *c, lls_workspace *w);
int lls_solve_spectral(gsl_vector *c, lls_workspace *w);
int lls_save(const char *filename, lls_workspace *w);
int lls_load(const char *filename, lls_workspace *w);

#endif /* INCLUDED_lls_h */
