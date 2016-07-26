/*
 * stage2.c
 *
 * 1. Read in spherical harmonic time series k_nm(t) from stage1 matrix file (nnm-by-nt)
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
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

#include <fftw3.h>
#include <lapacke/lapacke.h>

#include "common.h"
#include "green.h"
#include "lapack_wrapper.h"

#include "io.h"
#include "pca.h"
#include "window.h"

#define MAX_FREQ          100

/*
do_transform()
  Divide input time series into smaller segments and computed
windowed FFT of each segment. Store the complex coefficients
of the FFT for each time segment and each frequency into an
output matrix

Inputs:
        knm          - input time series, size nt
        rowidx       - row of Q matrices corresponding to (n,m)
        fs           - sampling frequency in 1/days
        window_size  - number of days per window
        window_shift - number of days to advance/slide forward
        Q            - (output) for each frequency omega_i,
                       Q[i](rowidx,:) is set to the FFT value at frequency
                       omega_i for each time segment
*/

static int
do_transform(FILE *fp, const gsl_vector *knm, const size_t rowidx,
             const double fs, const double window_size, const double window_shift,
             gsl_matrix_complex *Q[MAX_FREQ])
{
  const size_t nsamples = knm->size;              /* number of time samples */
  size_t start_idx = 0;
  size_t nwindow = (size_t) (window_size * fs);   /* optimal number of samples per window */
  size_t nforward = (size_t) (window_shift * fs); /* number of samples to slide forward */
  int done = 0;
  size_t k = 0;                                   /* current window index */
  size_t nfreq = nwindow / 2 + 1;
  gsl_vector *work = gsl_vector_alloc(nwindow);
  fftw_complex *fft_out = fftw_malloc(sizeof(fftw_complex) * nfreq);

  while (!done)
    {
      size_t end_idx = GSL_MIN(start_idx + nwindow - 1, nsamples - 1);
      size_t n = end_idx - start_idx + 1;         /* size of actual window */
      double sqrtn = sqrt((double) n);
      gsl_vector_const_view vknm = gsl_vector_const_subvector(knm, start_idx, n);
      gsl_vector_view vwork = gsl_vector_subvector(work, 0, n);
      fftw_plan plan;
      size_t i;

      /* apply window to current time segment */
      gsl_vector_memcpy(&vwork.vector, &vknm.vector);
      apply_ps1(&vwork.vector);

      /* compute windowed FFT */
      plan = fftw_plan_dft_r2c_1d(n, vwork.vector.data, fft_out, FFTW_ESTIMATE);
      fftw_execute(plan);

      /* loop over frequencies omega_i and store FFT output in Q[i] matrix */
      for (i = 0; i < nfreq; ++i)
        {
          gsl_vector_complex_view Qinm = gsl_matrix_complex_row(Q[i], rowidx);
          gsl_complex z = gsl_complex_rect(fft_out[i][0] / sqrtn, fft_out[i][1] / sqrtn);

          gsl_vector_complex_set(&Qinm.vector, k, z);
        }

      if (k == 0)
        {
          for (i = 0; i < n; ++i)
            {
              double freqi, powi, periodi;
              size_t idx;

              if (i <= n/2)
                idx = i;
              else
                idx = n - i;

              freqi = 2.0 * M_PI * idx * (fs / n);
              periodi = 2.0 * M_PI / freqi;
              powi = fft_out[idx][0]*fft_out[idx][0] +
                     fft_out[idx][1]*fft_out[idx][1];

              fprintf(fp, "%f %f %f %f\n",
                      gsl_vector_get(&vknm.vector, i),
                      gsl_vector_get(&vwork.vector, i),
                      periodi,
                      powi);
            }

          fprintf(fp, "\n\n");
        }

      fftw_destroy_plan(plan);

      ++k;

      start_idx += nforward;
      if (start_idx >= nsamples)
        done = 1;
    }

  gsl_vector_free(work);
  fftw_free(fft_out);

  return 0;
}

/*
count_windows()
  Count number of time segment windows in FFT analysis. So
for a sliding 2 day window with a 1 day overlap, set

window_size = 2
window_shift = 1

Inputs: nsamples     - total number of samples in time series
        fs           - sampling frequency in 1/days
        window_size  - number of days per window
        window_shift - number of days to advance/slide forward

Return: total number of windows
*/

size_t
count_windows(const size_t nsamples, const double fs,
              const double window_size, const double window_shift)
{
  size_t T = 0;                                   /* number of windows */
  size_t start_idx = 0;
  size_t nforward = (size_t) (window_shift * fs); /* number of samples to slide forward */
  int done = 0;

  while (!done)
    {
      ++T;

      start_idx += nforward;
      if (start_idx >= nsamples)
        done = 1;
    }

  return T;
}

int
main(int argc, char *argv[])
{
  const double R = R_EARTH_KM;
  const double fs = 24.0;         /* sample frequency in 1/days */
  size_t nmax, mmax;
  green_workspace *green_p;
  char *infile = PCA_STAGE1_KNM;
  gsl_matrix *A;
  double window_size = 2.0;       /* number of days in each time segment */
  double window_shift = 1.0;      /* number of days to shift forward in time */
  size_t nwindow = (size_t) (window_size * fs); /* number of samples per time window segment */
  size_t nfreq = nwindow / 2 + 1; /* number of frequencies returned from FFT */
  size_t nt, T, nnm;
  struct timeval tv0, tv1;

  /* nnm-by-T matrix storing power at frequency omega for each time segment and each
   * (n,m) channel, Q = [ X_1 X_2 ... X_T ] */
  gsl_matrix_complex *Q[MAX_FREQ];

  gsl_vector *S[MAX_FREQ];          /* singular values of Q */
  gsl_matrix_complex *U[MAX_FREQ];  /* left singular vectors of Q */
  gsl_matrix_complex *V[MAX_FREQ];  /* right singular vectors of Q */

  size_t i;

  if (nfreq > MAX_FREQ)
    {
      fprintf(stderr, "main: error: MAX_FREQ not large enough (%zu)\n", nfreq);
      exit(1);
    }

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
            window_size = atof(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s <-i stage1_matrix_file> [-t window_size (days)]\n", argv[0]);
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i stage1_matrix_file> [-t window_size (days)]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: reading %s...", PCA_STAGE1_DATA);
  pca_read_data(PCA_STAGE1_DATA, &nmax, &mmax);
  fprintf(stderr, "done (nmax = %zu mmax = %zu)\n", nmax, mmax);

  green_p = green_alloc(nmax, mmax, R);
  nnm = green_nnm(green_p);

  fprintf(stderr, "main: reading %s...", infile);
  A = pca_read_matrix(infile);
  fprintf(stderr, "done (%zu-by-%zu matrix)\n", A->size1, A->size2);

  /* compute number of time segments */
  nt = A->size2;
  T = count_windows(nt, fs, window_size, window_shift);

  fprintf(stderr, "main: time segment length: %g [days]\n", window_size);
  fprintf(stderr, "main: time segment slide:  %g [days]\n", window_shift);
  fprintf(stderr, "main: number of time segments: %zu\n", T);
  fprintf(stderr, "main: number of SH coefficients: %zu\n", nnm);
  fprintf(stderr, "main: number of frequencies: %zu\n", nfreq);

  /* allocate a matrix for each frequency */
  for (i = 0; i < nfreq; ++i)
    {
      Q[i] = gsl_matrix_complex_alloc(nnm, T);
      S[i] = gsl_vector_alloc(nnm);
      U[i] = gsl_matrix_complex_alloc(nnm, nnm);
      V[i] = gsl_matrix_complex_alloc(T, T);
    }

  {
    const char *fft_file = "fft_data.txt";
    FILE *fp = fopen(fft_file, "w");
    gsl_complex z = gsl_complex_rect(1.0 / sqrt((double) T), 0.0);
    size_t n;

    n = 1;
    fprintf(fp, "# Field %zu: knm(t) for first time segment\n", n++);
    fprintf(fp, "# Field %zu: Hamming-windowed knm(t) for first time segment\n", n++);
    fprintf(fp, "# Field %zu: Period (days)\n", n++);
    fprintf(fp, "# Field %zu: Power (nT^2)\n", n++);

    fprintf(stderr, "main: building the %zu-by-%zu Q matrices by performing windowed FFTs...",
            nnm, T);

    for (n = 1; n <= nmax; ++n)
      {
        int M = (int) GSL_MIN(n, mmax);
        int m;

        for (m = -M; m <= M; ++m)
          {
            size_t cidx = green_nmidx(n, m, green_p);
            gsl_vector_view knm = gsl_matrix_row(A, cidx);

            fprintf(fp, "# k(%zu,%d)\n", n, m);
            do_transform(fp, &knm.vector, cidx, fs, window_size, window_shift, Q);
          }
      }

    /* scale by 1/sqrt(T) - this acts as a weight factor */
    for (i = 0; i < nfreq; ++i)
      {
        gsl_matrix_complex_scale(Q[i], z);
      }

    fprintf(stderr, "done (data written to %s)\n", fft_file);

    fclose(fp);
  }

  {
    const size_t nmodes = 25; /* number of left singular vectors to output */
    int status;
    char buf[2048];

    for (i = 0; i < nfreq; ++i)
      {
        double freq = (fs / nwindow) * i;

        fprintf(stderr, "main: performing SVD of Q for frequency %g [cpd]...", freq);
        gettimeofday(&tv0, NULL);
        status = lapack_complex_svd(Q[i], S[i], U[i], V[i]);
        gettimeofday(&tv1, NULL);
        fprintf(stderr, "done (%g seconds, status = %d)\n",
                time_diff(tv0, tv1), status);

        sprintf(buf, "modes/S_%zu", i);
        fprintf(stderr, "main: writing singular values for frequency %g [cpd] in text format to %s...",
                freq, buf);
        pca_write_S(buf, nmax, mmax, freq, window_size, window_shift, S[i]);
        fprintf(stderr, "done\n");

        sprintf(buf, "modes/U_%zu", i);
        fprintf(stderr, "main: writing left singular vectors for frequency %g [cpd] in text format to %s...",
                freq, buf);
        pca_write_complex_U(buf, nmax, mmax, freq, window_size, window_shift, nmodes, U[i]);
        fprintf(stderr, "done\n");

        sprintf(buf, "modes/V_%zu", i);
        fprintf(stderr, "main: writing right singular vectors for frequency %g [cpd] in text format to %s...",
                freq, buf);
        pca_write_complex_V(buf, nmax, mmax, freq, window_size, window_shift, nmodes, V[i]);
        fprintf(stderr, "done\n");
      }

#if 0
    {
    const char *sval_file = PCA_STAGE2_SVAL;
    const char *U_file = PCA_STAGE2_U;
    const char *V_file = PCA_STAGE2_V;
    const char *Q_file = "data/stage2_Q.dat";

    fprintf(stderr, "main: writing singular values in binary format to %s...",
            sval_file);
    pca_write_vector(sval_file, S);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing left singular vectors in binary format to %s...",
            U_file);
    pca_write_matrix_complex(U_file, U);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing right singular vectors in binary format to %s...",
            V_file);
    pca_write_matrix_complex(V_file, V);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing Q matrix in binary format to %s...",
            Q_file);
    pca_write_matrix_complex(Q_file, Q);
    fprintf(stderr, "done\n");
    }
#endif
  }

  gsl_matrix_free(A);
  green_free(green_p);

  return 0;
}
