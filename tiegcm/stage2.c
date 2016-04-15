/*
 * stage2.c
 *
 * 1. Read in spherical harmonic time series from stage1 matrix file (nnm-by-nt)
 * 2. Divide each time series into T smaller segments and perform
 *    windowed Fourier transform of each segment
 * 3. Build Q(omega) matrix, nnm-by-T
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
#include <gsl/gsl_eigen.h>

#include <fftw3.h>
#include <lapacke/lapacke.h>

#include "common.h"
#include "green.h"

#include "io.h"
#include "lapack_wrapper.h"

static double
hamming_window(const double alpha, const double beta,
               const size_t n, const size_t N)
{
  double ratio = 2.0 * M_PI * n / ((double)N - 1.0);
  double w = alpha - beta*cos(ratio);

  return w;
}

static int
apply_window(gsl_vector *v)
{
  const size_t N = v->size;
  const double alpha = 0.53836;
  const double beta = 0.46164;
  size_t n;

  for (n = 0; n < N; ++n)
    {
      double wn = hamming_window(alpha, beta, n, N);
      double *ptr = gsl_vector_ptr(v, n);

      *ptr *= wn;
    }

  return 0;
}

/*
do_transform()
  Divide input time series into smaller segments and computed
windowed FFT of each segment. Store the complex coefficients
of the FFT for each time segment, with a fixed frequency omega,
into an output array

Inputs: omega - fixed frequency (per hour)
        qnm   - input time series, size nt
        Qnm   - (output) array of FFT values for each time
                segment at frequency omega, size T where T
                is the number of time segments
*/

static int
do_transform(FILE *fp, const double omega, const gsl_vector *qnm, gsl_vector_complex *Qnm)
{
  const size_t nsamples = qnm->size;                              /* number of time samples */
  const size_t T = Qnm->size;                                     /* number of time segments */
  const size_t nwindow = (size_t) ((double)nsamples / (double)T); /* size of window */
  const double fs = 1.0;                                          /* sampling rate = 1/hr */
  const double fny = 0.5 * fs;                                    /* Nyquist frequency (per hour) */
  size_t start_idx = 0;                                           /* index of starting sample for current time segment */
  size_t k;
  gsl_vector *work = gsl_vector_alloc(nwindow);

  /* loop over time segments */
  for (k = 0; k < T; ++k)
    {
      size_t end_idx = GSL_MIN(start_idx + nwindow - 1, nsamples - 1);
      size_t n = end_idx - start_idx + 1;
      double sqrtn = sqrt((double) n);
      gsl_vector_const_view vqnm = gsl_vector_const_subvector(qnm, start_idx, n);
      gsl_vector_view vwork = gsl_vector_subvector(work, 0, n);
      fftw_complex *fft_out = fftw_malloc(sizeof(fftw_complex) * (n/2 + 1));
      fftw_plan plan;
      size_t fft_idx = (n / 2) * (omega / fny); /* index of desired frequency in fft_out */
      gsl_complex z;

      gsl_vector_memcpy(&vwork.vector, &vqnm.vector);

      /* apply window to time segment */
      apply_window(&vwork.vector);

      plan = fftw_plan_dft_r2c_1d(n, vwork.vector.data, fft_out, FFTW_ESTIMATE);
      fftw_execute(plan);

      /* normalize by dividing by sqrt(n) */
      GSL_SET_COMPLEX(&z, fft_out[fft_idx][0] / sqrtn, fft_out[fft_idx][1] / sqrtn);
      gsl_vector_complex_set(Qnm, k, z);

      if (k == 0)
        {
          size_t i;

          for (i = 0; i < n; ++i)
            {
              double freqi, powi, periodi;
              size_t idx;

              if (i <= n/2 + 1)
                idx = i;
              else
                idx = n - i;

              freqi = 2.0 * M_PI * idx * (fs / n);
              periodi = 2.0 * M_PI / (freqi * 24.0);
              powi = fft_out[idx][0]*fft_out[idx][0] +
                     fft_out[idx][1]*fft_out[idx][1];

              fprintf(fp, "%f %f %f %f\n",
                      gsl_vector_get(&vqnm.vector, i),
                      gsl_vector_get(&vwork.vector, i),
                      periodi,
                      powi);
            }

          fprintf(fp, "\n\n");
        }

      fftw_destroy_plan(plan);
      fftw_free(fft_out);

      start_idx = end_idx + 1;
    }

  gsl_vector_free(work);

  return 0;
}

static int
build_sdm(const gsl_matrix_complex *Q, gsl_matrix_complex *S)
{
  const size_t T = Q->size2;
  size_t k;
  gsl_complex z;

  gsl_matrix_complex_set_zero(S);

  for (k = 0; k < T; ++k)
    {
      gsl_vector_complex_const_view Xk = gsl_matrix_complex_const_column(Q, k);

      /* S += X_k X_k^H */
      gsl_blas_zher(CblasLower, 1.0, &Xk.vector, S);
    }

  GSL_SET_COMPLEX(&z, 1.0 / (double)T, 0.0);
  gsl_matrix_complex_scale(S, z);

  return 0;
}

int
main(int argc, char *argv[])
{
  const size_t nmax = 60;                 /* maximum spherical harmonic degree */
  const size_t mmax = GSL_MIN(nmax, 30);  /* maximum spherical harmonic order */
  green_workspace *green_p = green_alloc(nmax, mmax);
  char *infile = "data/stage1_qnm.dat";
  gsl_matrix *A;
  double time_segment_length = 5.0;       /* number of days in each time segment */
  size_t nt, T, nnm;
  struct timeval tv0, tv1;

  /* nnm-by-T matrix storing power at frequency omega for each time segment and each
   * (n,m) channel */
  gsl_matrix_complex *Q;

  gsl_matrix_complex *S; /* spectral density matrix */

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:t:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 't':
            time_segment_length = atof(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s <-i stage1_matrix_file> [-t time_segment_length (days)]\n", argv[0]);
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i stage1_matrix_file> [-t time_segment_length (days)]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

  fprintf(stderr, "main: reading %s...", infile);
  A = read_matrix(infile);
  fprintf(stderr, "done (%zu-by-%zu matrix)\n", A->size1, A->size2);

  /* compute number of time segments */
  nt = A->size2;
  T = (nt / 24.0) / time_segment_length;

  fprintf(stderr, "main: time segment length: %g [days]\n", time_segment_length);
  fprintf(stderr, "main: number of time segments: %zu\n", T);

  /*nnm = A->size1;*/
  nnm = nmax*(nmax+2);

  Q = gsl_matrix_complex_alloc(nnm, T);
  S = gsl_matrix_complex_alloc(nnm, nnm);

  {
    const char *fft_file = "fft_data.txt";
    FILE *fp = fopen(fft_file, "w");
    const double omega = 1.0 / 24.0; /* 1 cpd */
    size_t n;

    n = 1;
    fprintf(fp, "# Field %zu: qnm(t) for first time segment\n", n++);
    fprintf(fp, "# Field %zu: Hamming-windowed qnm(t) for first time segment\n", n++);
    fprintf(fp, "# Field %zu: Period (days)\n", n++);
    fprintf(fp, "# Field %zu: Power (nT^2)\n", n++);

    fprintf(stderr, "main: building the Q matrix by performing windowed FFTs...");

    for (n = 1; n <= nmax; ++n)
      {
        int m, ni = (int) n;

        for (m = -ni; m <= ni; ++m)
          {
            size_t cidx = green_nmidx(n, m, green_p);
            gsl_vector_view qnm = gsl_matrix_row(A, cidx);
            gsl_vector_complex_view Qnm = gsl_matrix_complex_row(Q, cidx);

            fprintf(fp, "# q(%zu,%d)\n", n, m);
            do_transform(fp, omega, &qnm.vector, &Qnm.vector);
          }
      }

    fprintf(stderr, "done (data written to %s)\n", fft_file);

    fclose(fp);
  }

  fprintf(stderr, "main: building spectral density matrix...");
  build_sdm(Q, S);
  fprintf(stderr, "done\n");

  {
    const size_t N = S->size1;
    const char *eval_txt_file = "eval.txt";
    const char *eval_file = "data/stage2_eval.dat";
    const char *evec_file = "data/stage2_evec.dat";
    const char *Q_file = "data/stage2_Q.dat";
    gsl_vector *eval = gsl_vector_alloc(N);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(N, N);
    int status, eval_found;
    FILE *fp;

    fprintf(stderr, "main: performing eigendecomposition of SDM...");
    gettimeofday(&tv0, NULL);
    status = lapack_eigen_herm(S, eval, evec, &eval_found);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds, status = %d, %d eigenvalues found)\n",
            time_diff(tv0, tv1), status, eval_found);

    fprintf(stderr, "main: writing eigenvalues in text format to %s...",
            eval_txt_file);
    fp = fopen(eval_txt_file, "w");
    gsl_vector_fprintf(fp, eval, "%.12e");
    fclose(fp);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing eigenvalues in binary format to %s...",
            eval_file);
    write_vector(eval_file, eval);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing eigenvectors in binary format to %s...",
            evec_file);
    write_matrix_complex(evec_file, evec);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing Q matrix in binary format to %s...",
            Q_file);
    write_matrix_complex(Q_file, Q);
    fprintf(stderr, "done\n");

    gsl_vector_free(eval);
    gsl_matrix_complex_free(evec);
  }

  gsl_matrix_free(A);
  gsl_matrix_complex_free(Q);
  green_free(green_p);

  return 0;
}
