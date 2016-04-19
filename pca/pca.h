/*
 * pca.h
 */

#ifndef INCLUDED_pca_h
#define INCLUDED_pca_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

#include "green.h"

/* Stage1: knm matrix file, nnm-by-nt */
#define PCA_STAGE1_KNM     "/data/palken/lib/pca/data/stage1_knm.dat"
#define PCA_STAGE1_DATA    "/data/palken/lib/pca/data/stage1_data.dat"

/* Stage2b: singular values, left and right singular vectors data files */
#define PCA_STAGE2B_SVAL   "/data/palken/lib/pca/data/stage2b_sval.dat"
#define PCA_STAGE2B_U      "/data/palken/lib/pca/data/stage2b_U.dat"
#define PCA_STAGE2B_V      "/data/palken/lib/pca/data/stage2b_V.dat"

typedef struct
{
  size_t nmax;     /* maximum spherical harmonic degree */
  size_t mmax;     /* maximum spherical harmonic order */
  size_t nnm;      /* number of spherical harmonic coefficients */
  double R;        /* reference radius (km) */
  double b;        /* current shell radius (km) */
  gsl_matrix *K;   /* knm matrix, nnm-by-nt */
  gsl_vector *S;   /* singular values */
  gsl_matrix *U;   /* left singular vectors in terms of knm, nnm-by-nnm */
  gsl_matrix *V;   /* right singular vectors in terms of knm, nt-by-nt */
  gsl_matrix *G;   /* left singular vectors in terms of gnm, nnm-by-nnm */

  double *X;       /* temporary workspace, size nnm */
  double *Y;
  double *Z;

  green_workspace *green_workspace_p;
} pca_workspace;

/*
 * Prototypes
 */

pca_workspace *pca_alloc();
void pca_free(pca_workspace *w);
int pca_variance(const char *filename, const double thresh, size_t * nsing,
                 pca_workspace *w);
int pca_print_map(const char *filename, const double r, const size_t pcidx,
                  pca_workspace *w);

#endif /* INCLUDED_pca_h */
