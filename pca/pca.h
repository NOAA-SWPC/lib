/*
 * pca.h
 */

#ifndef INCLUDED_pca_h
#define INCLUDED_pca_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

#include "green.h"

/* Stage1: knm matrix file, nnm-by-nt */
#define PCA_STAGE1_KNM         "/data/palken/lib/pca/data/stage1_knm.dat"
#define PCA_STAGE1_DATA        "/data/palken/lib/pca/data/stage1_data.dat"

/* Stage2b: singular values, left and right singular vectors data files
 * These are file prefixes, the string %02zuUT.dat must be added to the end, to
 * specify UT
 */
#define PCA_STAGE2B_SVAL_TXT   "/data/palken/lib/pca/data/stage2b_sval_time"
#define PCA_STAGE2B_SVAL       "/data/palken/lib/pca/data/stage2b_sval"
#define PCA_STAGE2B_U          "/data/palken/lib/pca/data/stage2b_U"
#define PCA_STAGE2B_V          "/data/palken/lib/pca/data/stage2b_V"

typedef struct
{
  size_t nmax;       /* maximum spherical harmonic degree */
  size_t mmax;       /* maximum spherical harmonic order */
  size_t nnm;        /* number of spherical harmonic coefficients */
  double R;          /* reference radius (km) */
  double b;          /* current shell radius (km) */
  gsl_matrix *K;     /* knm matrix, nnm-by-nt */

  /* these arrays contain PC components for each UT hour [0,23] */
  gsl_vector *S[24]; /* singular values */
  gsl_matrix *U[24]; /* left singular vectors in terms of knm, nnm-by-nnm */
  gsl_matrix *G[24]; /* left singular vectors in terms of gnm, nnm-by-nnm */

  size_t ut;         /* current UT hour in [0,23] */

  double *X;         /* temporary workspace, size nnm */
  double *Y;
  double *Z;
  gsl_vector *work;  /* size nnm */

  green_workspace *green_workspace_p;
} pca_workspace;

/*
 * Prototypes
 */

pca_workspace *pca_alloc();
void pca_free(pca_workspace *w);
int pca_set_UT(const size_t ut, pca_workspace *w);
int pca_sval(gsl_vector * T, pca_workspace *w);
int pca_variance(const char *filename, const double thresh, size_t * nsing,
                 pca_workspace *w);
int pca_print_map(FILE *fp, const double r, const gsl_vector * alpha,
                  pca_workspace *w);
int pca_pc_B(const size_t pcidx, const double r, const double theta, const double phi,
             double B[3], pca_workspace *w);
int pca_B(const gsl_vector *alpha, const double r, const double theta, const double phi,
          double B[3], pca_workspace *w);
double pca_chi(const gsl_vector *alpha, const double theta, const double phi, pca_workspace *w);
int pca_K(const gsl_vector *alpha, const double theta, const double phi, double K[3], pca_workspace *w);
int pca_print_pc_map(const char *filename, const double r, const size_t pcidx,
                     pca_workspace *w);

#endif /* INCLUDED_pca_h */
