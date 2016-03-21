/*
 * stage3.c
 *
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
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

#include "common.h"
#include "green.h"

#include "io.h"
#include "lapack_common.h"

int
main(int argc, char *argv[])
{
  size_t nnm;
  char *qnm_file = "data/stage1_qnm.dat";
  char *eval_file = "data/stage2_eval.dat";
  char *evec_file = "data/stage2_evec.dat";
  char *Q_file = "data/stage2_Q.dat";
  struct timeval tv0, tv1;

  /* nnm-by-T matrix storing power at frequency omega for each time segment and each
   * (n,m) channel */
  gsl_matrix_complex *Q;
  gsl_matrix_complex *Qt;    /* Q~ = U*alpha */

  gsl_matrix_complex *evec;  /* eigenvectors of SDM */
  gsl_vector *eval;          /* eigenvalues of SDM */

  gsl_matrix_complex *alpha; /* alpha matrix, P-by-T */
  gsl_matrix *qnm;           /* qnm(t) matrix */

  size_t T;                  /* number of time segments */
  size_t P;                  /* number of principal eigenvectors to use (<= T) */
  gsl_matrix_complex_view U; /* largest P eigenvectors of SDM */
  int status, rank;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          default:
            fprintf(stderr, "Usage: %s <-i stage1_matrix_file> [-t time_segment_length (days)]\n", argv[0]);
            break;
        }
    }

  fprintf(stderr, "main: reading qnm matrix from %s...", qnm_file);
  qnm = read_matrix(qnm_file);
  fprintf(stderr, "done (%zu-by-%zu matrix read)\n", qnm->size1, qnm->size2);

  fprintf(stderr, "main: reading eigenvalues from %s...", eval_file);
  eval = read_vector(eval_file);
  fprintf(stderr, "done (%zu eigenvalues read)\n", eval->size);

  fprintf(stderr, "main: reading eigenvectors from %s...", evec_file);
  evec = read_matrix_complex(evec_file);
  fprintf(stderr, "done (%zu-by-%zu matrix read)\n", evec->size1, evec->size2);

  fprintf(stderr, "main: reading Q matrix from %s...", Q_file);
  Q = read_matrix_complex(Q_file);
  fprintf(stderr, "done (%zu-by-%zu matrix read)\n", Q->size1, Q->size2);

  T = Q->size2;
  P = GSL_MIN(20, T);
  nnm = evec->size1;

  fprintf(stderr, "main: using %zu largest eigenvectors\n", P);

  alpha = gsl_matrix_complex_alloc(P, T);
  Qt = gsl_matrix_complex_alloc(nnm, T);

  /* select the P largest eigenvectors of SDM */
  U = gsl_matrix_complex_submatrix(evec, 0, nnm - P, nnm, P);

  /* solve: U*alpha = Q */
  fprintf(stderr, "main: solving PCA problem for alpha...");
  gettimeofday(&tv0, NULL);
  status = lapack_complex_lls(&U.matrix, Q, alpha, &rank);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, status = %d, rank = %d)\n",
          time_diff(tv0, tv1), status, rank);

  /* compute: Q~ = U*alpha */
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, &U.matrix, alpha,
                 GSL_COMPLEX_ZERO, Qt);

  {
    size_t k;
    size_t cidx = green_nmidx(2, 0);

    for (k = 0; k < T; ++k)
      {
        gsl_complex Qnm = gsl_matrix_complex_get(Q, cidx, k);
        gsl_complex Qtnm = gsl_matrix_complex_get(Qt, cidx, k);

        printf("%f %f %f %f\n",
               GSL_REAL(Qnm),
               GSL_IMAG(Qnm),
               GSL_REAL(Qtnm),
               GSL_IMAG(Qtnm));
      }
  }

  gsl_matrix_complex_free(Q);
  gsl_matrix_complex_free(evec);
  gsl_matrix_complex_free(alpha);
  gsl_matrix_complex_free(Qt);
  gsl_matrix_free(qnm);
  gsl_vector_free(eval);

  return 0;
}
