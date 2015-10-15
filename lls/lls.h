/*
 * lls.h
 * Patrick Alken
 */

#ifndef INCLUDED_lls_h
#define INCLUDED_lls_h

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_eigen.h>

typedef struct
{
  size_t p;                /* number of coefficients */
  size_t max_block;        /* maximum observations at a time */
  gsl_matrix *ATA;         /* A^T W A */
  gsl_vector *ATb;         /* A^T W b */
  double bTb;              /* b^T W b */
  gsl_vector *c;           /* coefficient vector */
  double residual;
  double chisq;            /* residual ssq chi^2 */
  double cond;             /* condition number */

  gsl_vector *r;           /* residuals for robust regression */
  gsl_vector *w_robust;    /* robust weights */
  size_t niter;            /* number of robust iterations */

  /* extra workspace for LAPACK */
  gsl_matrix *work_A;
  gsl_vector *work_b;
  gsl_matrix *AF;
  gsl_vector *S;

  /* for computing L-curve */
  gsl_eigen_symm_workspace *eigen_p;
  gsl_vector *eval;

  gsl_multifit_robust_workspace *robust_workspace_p;
} lls_workspace;

typedef struct
{
  size_t p;                      /* number of coefficients */
  size_t max_block;              /* maximum observations at a time */
  gsl_vector *c;                 /* coefficient vector */
  double residual;
  double chisq;                  /* residual ssq chi^2 */
  double cond;                   /* condition number */

  gsl_matrix_complex *AHA;       /* A^H W A */
  gsl_vector_complex *AHb;       /* A^H W b */
  double bHb;                    /* rhs scalar product: b^H W b */

  gsl_vector_complex *c_complex; /* coefficient vector for robust iterations */
  gsl_vector_complex *r_complex; /* complex residuals for robust regression */
  gsl_vector *r;                 /* residuals for robust regression */
  gsl_vector *w_robust;          /* robust weights */
  size_t niter;                  /* number of robust iterations */

  /* extra workspace for LAPACK */
  gsl_matrix_complex *work_A_complex;
  gsl_vector_complex *work_b_complex;
  gsl_matrix_complex *AF;
  gsl_vector *S;

  gsl_multifit_robust_workspace *robust_workspace_p;
} lls_complex_workspace;

/*
 * Prototypes
 */

lls_workspace *lls_alloc(const size_t max_block, const size_t ncoeff);
void lls_free(lls_workspace *w);
int lls_reset(lls_workspace *w);
int lls_fold(gsl_matrix *A, gsl_vector *b,
             gsl_vector *wts, lls_workspace *w);
int lls_regularize(const double lambda, gsl_matrix *ATA);
int lls_regularize2(const gsl_vector *diag, lls_workspace *w);
int lls_solve(const double lambda, gsl_vector *c, lls_workspace *w);
int lls_lcurve(gsl_vector *reg_param, gsl_vector *rho, gsl_vector *eta,
               lls_workspace *w);
int lls_solve_spectral(gsl_vector *c, lls_workspace *w);
int lls_save(const char *filename, lls_workspace *w);
int lls_load(const char *filename, lls_workspace *w);

/* lls_complex.c */
lls_complex_workspace *lls_complex_alloc(const size_t max_block, const size_t p);
void lls_complex_free(lls_complex_workspace *w);
int lls_complex_reset(lls_complex_workspace *w);
int lls_complex_fold(gsl_matrix_complex *A, gsl_vector_complex *b,
                     gsl_vector *wts, lls_complex_workspace *w);
int lls_complex_solve(gsl_vector_complex *c, lls_complex_workspace *w);
int lls_complex_regularize2(const gsl_vector *diag, lls_complex_workspace *w);
int lls_complex_invert(gsl_matrix_complex *B, const lls_complex_workspace *w);
int lls_complex_correlation(gsl_matrix_complex *B, const lls_complex_workspace *w);
int lls_complex_save(const char *filename, lls_complex_workspace *w);
int lls_complex_load(const char *filename, lls_complex_workspace *w);

/* lls_lapack.c */
int lls_lapack_dposv(const double lambda, gsl_vector *c, lls_workspace *w);
int lls_lapack_zposv(gsl_vector_complex *c, lls_complex_workspace *w);
int lls_lapack_zinvert(gsl_matrix_complex *B, const lls_complex_workspace *w);

#endif /* INCLUDED_lls_h */
