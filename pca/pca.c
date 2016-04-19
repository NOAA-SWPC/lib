/*
 * pca.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_legendre.h>

#include "green.h"

#include "io.h"
#include "pca.h"

static int pca_calc_G(const double b, pca_workspace *w);
static double pca_chi(const double b, const double theta, const double phi,
                      const gsl_vector *k, pca_workspace *w);

pca_workspace *
pca_alloc()
{
  pca_workspace *w;

  w = calloc(1, sizeof(pca_workspace));
  if (!w)
    return 0;

  fprintf(stderr, "pca_alloc: reading %s...", PCA_STAGE1_DATA);
  pca_read_data(PCA_STAGE1_DATA, &(w->nmax), &(w->mmax));
  fprintf(stderr, "done (nmax = %zu mmax = %zu)\n", w->nmax, w->mmax);

  fprintf(stderr, "pca_alloc: reading %s...", PCA_STAGE1_KNM);
  w->K = pca_read_matrix(PCA_STAGE1_KNM);
  fprintf(stderr, "done (%zu-by-%zu matrix)\n", w->K->size1, w->K->size2);

  fprintf(stderr, "pca_alloc: reading singular values from %s...", PCA_STAGE2B_SVAL);
  w->S = pca_read_vector(PCA_STAGE2B_SVAL);
  fprintf(stderr, "done (%zu singular values read)\n", w->S->size);

  fprintf(stderr, "pca_alloc: reading left singular vectors from %s...", PCA_STAGE2B_U);
  w->U = pca_read_matrix(PCA_STAGE2B_U);
  fprintf(stderr, "done (%zu-by-%zu matrix read)\n", w->U->size1, w->U->size2);

  w->green_workspace_p = green_alloc(w->nmax, w->mmax);

  w->nnm = green_nnm(w->green_workspace_p);
  w->R = w->green_workspace_p->R;
  w->b = w->R + 110.0;

  w->X = malloc(w->nnm * sizeof(double));
  w->Y = malloc(w->nnm * sizeof(double));
  w->Z = malloc(w->nnm * sizeof(double));

  w->G = gsl_matrix_alloc(w->nnm, w->nnm);

  pca_calc_G(w->b, w);

  return w;
}

void
pca_free(pca_workspace *w)
{
  if (w->K)
    gsl_matrix_free(w->K);

  if (w->G)
    gsl_matrix_free(w->G);

  if (w->S)
    gsl_vector_free(w->S);

  if (w->U)
    gsl_matrix_free(w->U);

  if (w->V)
    gsl_matrix_free(w->V);

  if (w->green_workspace_p)
    green_free(w->green_workspace_p);

  if (w->X)
    free(w->X);

  if (w->Y)
    free(w->Y);

  if (w->Z)
    free(w->Z);

  free(w);
}

/*
pca_variance()
  Print cumulative variance based on singular values to a file.
Also determine number of singular vectors required to explain
a desired threshold of variance in the data.

Inputs: filename - output file
        S        - singular values
        thresh   - desired variance threshold in [0,1]
        nsing    - (output) number of singular values required to
                   explain a variance of 'thresh'
*/

int
pca_variance(const char *filename, const double thresh, size_t * nsing,
             pca_workspace *w)
{
  int s = 0;
  const gsl_vector *S = w->S;
  const size_t n = S->size;
  FILE *fp = fopen(filename, "w");
  size_t i;
  double sum_all = 0.0;
  double cumsum = 0.0;
  size_t nvec = 0;

  i = 1;
  fprintf(fp, "# Field %zu: eigenvalue number i\n", i++);
  fprintf(fp, "# Field %zu: ( sum_{j=1}^i lambda_j ) / ( sum_j lambda_j )\n", i++);

  /* compute sum_i lambda_i, with lambda_i = sigma_i^2 */
  gsl_blas_ddot(S, S, &sum_all);

  for (i = 0; i < n; ++i)
    {
      double sigma = gsl_vector_get(S, i);
      double cumvar;

      cumsum += sigma * sigma;
      cumvar = cumsum / sum_all;

      /* compute number of singular vectors / eigenvectors needed
       * to explain desired variance threshold */
      if (nvec == 0 && cumvar > thresh)
        nvec = i + 1;

      fprintf(fp, "%zu %.12e\n", i + 1, cumvar);
    }

  fclose(fp);

  *nsing = nvec;

  return s;
}

/*
pca_print_map()
  Print lat/lon map of one PC

Inputs: filename - output file
        r        - radius of shell for magnetic field maps
        pcidx    - which principal component, in [0,nnm-1]
        w        - workspace
*/

int
pca_print_map(const char *filename, const double r, const size_t pcidx,
              pca_workspace *w)
{
  int s = 0;
  const size_t nnm = w->nnm;
  const double b = w->b; /* radius of current shell */
  FILE *fp;
  double lat, lon;
  green_workspace *green_p = w->green_workspace_p;
  gsl_vector_view Xv = gsl_vector_view_array(w->X, nnm);
  gsl_vector_view Yv = gsl_vector_view_array(w->Y, nnm);
  gsl_vector_view Zv = gsl_vector_view_array(w->Z, nnm);
  gsl_vector_const_view pc_knm = gsl_matrix_const_column(w->U, pcidx);
  gsl_vector_const_view pc_gnm = gsl_matrix_const_column(w->G, pcidx);
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "pca_print_map: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Principle component %zu/%zu\n", pcidx + 1, w->nnm);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: chi (kA / nT)\n", i++);
  fprintf(fp, "# Field %zu: B_x (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: B_y (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: B_z (dimensionless)\n", i++);

  for (lon = -180.0; lon <= 180.0; lon += 1.0)
    {
      double phi = lon * M_PI / 180.0;

      for (lat = -89.9; lat <= 89.9; lat += 1.0)
        {
          double theta = M_PI / 2.0 - lat * M_PI / 180.0;
          double B_pc[3];
          double chi;

          chi = pca_chi(b, theta, phi, &pc_knm.vector, w);

          /*
           * If r < b, the current shell is an external source so
           * we can directly use the knm coefficients in the U matrix.
           *
           * If r > b, the current shell is an internal source, and
           * we must first compute the gnm coefficients from knm (done
           * in pca_alloc via pca_calc_G, and then use internal Green's functions
           * for the dot product
           */
          if (r <= b)
            {
              green_calc_ext(r, theta, phi, w->X, w->Y, w->Z, green_p);

              gsl_blas_ddot(&pc_knm.vector, &Xv.vector, &B_pc[0]);
              gsl_blas_ddot(&pc_knm.vector, &Yv.vector, &B_pc[1]);
              gsl_blas_ddot(&pc_knm.vector, &Zv.vector, &B_pc[2]);
            }
          else
            {
              green_calc_int(r, theta, phi, w->X, w->Y, w->Z, green_p);

              gsl_blas_ddot(&pc_gnm.vector, &Xv.vector, &B_pc[0]);
              gsl_blas_ddot(&pc_gnm.vector, &Yv.vector, &B_pc[1]);
              gsl_blas_ddot(&pc_gnm.vector, &Zv.vector, &B_pc[2]);
            }

          fprintf(fp, "%f %f %f %f %f %f\n",
                  lon,
                  lat,
                  chi,
                  B_pc[0],
                  B_pc[1],
                  B_pc[2]);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

/*
pca_calc_G()
  Convert PCs in a basis of external Gauss coefficients representing
a current shell at radius b into internal Gauss coefficients, suitable
for computing the magnetic field above the current shell.

gnm = - (n / (n+1)) (b/R)^{2n + 1} knm

Inputs: b - radius of current shell (km)
        w - workspace
*/

static int
pca_calc_G(const double b, pca_workspace *w)
{
  size_t n;
  green_workspace *green_p = w->green_workspace_p;

  for (n = 1; n <= w->nmax; ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax);
      int m;
      double nfac = -(double)n / (n + 1.0);
      double rfac = pow(b / w->R, 2.0*n + 1.0);

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = green_nmidx(n, m, green_p);
          gsl_vector_view k = gsl_matrix_row(w->U, cidx);
          gsl_vector_view g = gsl_matrix_row(w->G, cidx);

          gsl_vector_memcpy(&g.vector, &k.vector);
          gsl_vector_scale(&g.vector, nfac * rfac);
        }
    }

  return 0;
}

/*
pca_chi()
  Compute current stream function
*/

static double
pca_chi(const double b, const double theta, const double phi,
        const gsl_vector *k, pca_workspace *w)
{
  const double mu0 = 400.0 * M_PI; /* units of nT / (kA km^{-1}) */
  const size_t nmax = w->nmax;
  const double ratio = b / w->R;
  green_workspace *green_p = w->green_workspace_p;
  double *Plm = green_p->Plm;
  double *dPlm = green_p->dPlm;
  size_t n;
  double chi = 0.0;
  double rfac = 1.0;

  /* compute associated legendres */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax, cos(theta), Plm, dPlm);

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax);
      int m;
      double nfac = (2.0 * n + 1.0) / (n + 1.0);

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t cidx = green_nmidx(n, m, green_p);
          double knm = gsl_vector_get(k, cidx);
          double qnm = nfac * rfac * knm;
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          double Snm = Plm[pidx];

          if (m < 0)
            Snm *= sin(mabs * phi);
          else
            Snm *= cos(mabs * phi);

          chi += qnm * Snm;
        }

      /* (b/R)^{n-2} */
      rfac *= ratio;
    }

  /* units of kA */
  chi *= -(b / mu0);

  return chi;
}
