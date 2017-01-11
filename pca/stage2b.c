/*
 * stage2b.c
 *
 * 1. Read in spherical harmonic time series knm from stage1 matrix file (nnm-by-nt)
 * 2. Subtract mean from each knm(t) time series
 * 2. Compute SVD of centered knm matrix
 * 3. Write singular values and singular vectors to disk
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <getopt.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>

#include <fftw3.h>
#include <lapacke/lapacke.h>

#include "common.h"
#include "green.h"
#include "lapack_wrapper.h"
#include "oct.h"

#include "io.h"
#include "pca.h"

#include "pca_common.c"

/* define to subtract mean from knm(t) time series prior to SVD */
#define SUBTRACT_MEAN           1

/* define to taper the high degree SH coefficients with a cosine taper (to correct
 * TIEGCM ringing issue) */
#define TAPER_COEFFS            0

/* define to construct and output covariance matrix
 * 2 Jan 2017: after plotting covariance and correlation matrix,
 * there isn't much structure to see */
#define PLOT_COVAR              0

int
subtract_mean(gsl_matrix * K)
{
  const size_t N = K->size1;
  const size_t M = K->size2;
  size_t i;

  for (i = 0; i < N; ++i)
    {
      gsl_vector_view v = gsl_matrix_row(K, i);
      double mean = gsl_stats_mean(v.vector.data, v.vector.stride, M);
      gsl_vector_add_constant(&v.vector, -mean);
    }

  return 0;
}

int
correlation_matrix(gsl_matrix * A)
{
  const size_t N = A->size1;
  size_t i, j;
  gsl_vector *d = gsl_vector_alloc(N);

  /* store inverse std deviations D = diag(A)^{-1/2} and set diag(A) = 1 */
  for (i = 0; i < N; ++i)
    {
      double *Aii = gsl_matrix_ptr(A, i, i);

      gsl_vector_set(d, i, 1.0 / sqrt(*Aii));
      *Aii = 1.0;
    }

  /* compute correlation matrix corr = D A D, using lower triangle */
  for (i = 0; i < N; ++i)
    {
      double di = gsl_vector_get(d, i);

      for (j = 0; j < i; ++j)
        {
          double dj = gsl_vector_get(d, j);
          double *Aij = gsl_matrix_ptr(A, i, j);

          *Aij *= di * dj;
        }
    }

  gsl_vector_free(d);

  return 0;
}

int
main(int argc, char *argv[])
{
  size_t nnm;                             /* number of spherical harmonic time series */
  size_t nt;                              /* number of time steps */
  char *infile = PCA_STAGE1_KNM;
  gsl_matrix *K;
  struct timeval tv0, tv1;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          default:
            fprintf(stderr, "Usage: %s <-i stage1_matrix_file>\n", argv[0]);
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i stage1_matrix_file>\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

  fprintf(stderr, "main: reading %s...", infile);
  K = pca_read_matrix(infile);
  fprintf(stderr, "done (%zu-by-%zu matrix)\n", K->size1, K->size2);

  nnm = K->size1;
  nt = K->size2;

#if TAPER_COEFFS

  /* taper spherical harmonic coefficients to eliminate ringing effect from TIEGCM */
  {
    green_workspace *green_p = green_alloc(60, 30, R_EARTH_KM);
    char *spectrum_file = "spectrum_taper.s";
    gsl_vector_view x = gsl_matrix_column(K, 0);

    fprintf(stderr, "main: tapering knm coefficients...");
    taper_knm(20, K, green_p);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing tapered spectrum to %s...", spectrum_file);
    green_print_spectrum(spectrum_file, &x.vector, green_p);
    fprintf(stderr, "done\n");

    green_free(green_p);
  }

#endif

#if SUBTRACT_MEAN
  fprintf(stderr, "main: subtracting mean from each knm time series...");
  subtract_mean(K);
  fprintf(stderr, "done\n");
#endif

  /* compute 1/sqrt(nt) K for SVD computation */
  gsl_matrix_scale(K, 1.0 / sqrt((double) nt));

#if PLOT_COVAR
  /* build and output covariance matrix */
  {
    const char *covar_file = "covar.dat";
    const char *eval_file = "eval_time.txt";
    gsl_matrix *covar = gsl_matrix_alloc(nnm, nnm);
    gsl_vector *eval = gsl_vector_alloc(nnm);
    int eval_found;
    FILE *fp;

    fprintf(stderr, "main: constructing covariance matrix...");
    gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, K, 0.0, covar);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: computing eigenvalues of covariance matrix...");
    lapack_eigen_symm(covar, eval, &eval_found);
    fprintf(stderr, "done (%d eigenvalues found)\n", eval_found);

    fprintf(stderr, "main: writing eigenvalues format to %s...", eval_file);
    fp = fopen(eval_file, "w");
    gsl_vector_fprintf(fp, eval, "%.12e");
    fclose(fp);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: converting to correlation matrix...");
    correlation_matrix(covar);
    fprintf(stderr, "done\n");

    /* copy lower triangle to upper */
    gsl_matrix_transpose_tricpy('L', 0, covar, covar);

    fprintf(stderr, "main: printing correlation matrix to %s...", covar_file);
    octave_plot(covar, covar_file);
    fprintf(stderr, "done\n");

    gsl_matrix_free(covar);
    gsl_vector_free(eval);
  }
#endif /* PLOT_COVAR */

  {
    const char *sval_txt_file = "sval_time.txt";
    const char *sval_file = PCA_STAGE2B_SVAL;
    const char *U_file = PCA_STAGE2B_U;
    const char *V_file = PCA_STAGE2B_V;
    gsl_vector *S = gsl_vector_alloc(GSL_MIN(nnm, nt));
    gsl_matrix *U = gsl_matrix_alloc(nnm, nnm);
    gsl_matrix *V = gsl_matrix_alloc(nt, nt);
    int status;
    FILE *fp;

    fprintf(stderr, "main: performing SVD of K matrix (%zu-by-%zu)...", nnm, nt);
    gettimeofday(&tv0, NULL);
    status = lapack_svd(K, S, U, V);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds, status = %d)\n",
            time_diff(tv0, tv1), status);

    fprintf(stderr, "main: writing singular values in text format to %s...",
            sval_txt_file);
    fp = fopen(sval_txt_file, "w");
    gsl_vector_fprintf(fp, S, "%.12e");
    fclose(fp);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing singular values in binary format to %s...",
            sval_file);
    pca_write_vector(sval_file, S);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing left singular vectors in binary format to %s...",
            U_file);
    pca_write_matrix(U_file, U);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing right singular vectors in binary format to %s...",
            V_file);
    pca_write_matrix(V_file, V);
    fprintf(stderr, "done\n");

    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_vector_free(S);
  }

  gsl_matrix_free(K);

  return 0;
}
