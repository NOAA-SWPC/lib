/*
 * stage3a.c
 *
 * Use results of FFT analysis (stage2a) on J grids to build
 * spectral density matrix and compute SVD
 *
 * ./stage3a [-i fft_data_file]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <assert.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <common/common.h>
#include <common/bsearch.h>

#include "lapack_wrapper.h"

#include "io.h"
#include "pca3d.h"
#include "tiegcm3d.h"

/*
build_X()

  Build matrix X (3N-by-T, with grid size N = nr * ntheta * nphi)
for a given frequency

X(w_p) = [ X_1(w_p) X_2(w_p) ... X_T(w_p) ]

with length 3N column vectors

X_t(w_p) = [ Qr^{t}_{ir,ilat,ilon}(w_p) ]
           [ Qt^{t}_{ir,ilat,ilon}(w_p) ]
           [ Qp^{t}_{ir,ilat,ilon}(w_p) ]

Inputs: ifreq - index of desired frequency
        data  - fft data
        X     - (output) matrix, 3N-by-T

Return: success/error
*/

int
build_X(const size_t ifreq, const tiegcm3d_fft_data * data, gsl_matrix_complex * X)
{
  int s = 0;
  const size_t T = data->T;   /* number of time window segments */
  const double fs = data->fs;
  const size_t nwindow = (size_t) (data->window_size * fs);   /* optimal number of samples per window */
  const size_t nfreq = nwindow / 2 + 1;                       /* FFT output buffer size (number of frequencies) */
  size_t t;

  /*
   * The time window index t is the column index of X; the spatial indices
   * (ir,ilat,ilon) give the row index of X
   */
  for (t = 0; t < T; ++t)
    {
      size_t ir, ilat, ilon;

      for (ir = 0; ir < data->nr; ++ir)
        {
          for (ilat = 0; ilat < data->nlat; ++ilat)
            {
              for (ilon = 0; ilon < data->nlon; ++ilon)
                {
                  size_t row_idx = CIDX3(ir, data->nr, ilat, data->nlat, ilon, data->nlon);
                  size_t fft_idx = TIEGCM3D_FREQIDX(t, ifreq, ir, ilat, ilon, data, T, nfreq);
                  gsl_complex Qr = data->Qr[fft_idx];
                  gsl_complex Qt = data->Qt[fft_idx];
                  gsl_complex Qp = data->Qp[fft_idx];

                  gsl_matrix_complex_set(X, row_idx, t, Qr);
                  gsl_matrix_complex_set(X, 2*row_idx, t, Qt);
                  gsl_matrix_complex_set(X, 3*row_idx, t, Qp);
                }
            }
        }
    }

  return s;
}

int
main(int argc, char *argv[])
{
  tiegcm3d_fft_data data;
  struct timeval tv0, tv1;
  char *infile = PCA3D_STAGE2A_FFT_DATA;

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
            fprintf(stderr, "Usage: %s [-i fft_data_file]\n", argv[0]);
            exit(1);
            break;
        }
    }

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);
  data = pca3d_read_fft_data(infile);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  {
    int status;
    const double window_size = data.window_size;
    const double window_shift = data.window_shift;
    const double freq = 1.0;                                  /* desired frequency in cpd */
    const size_t ifreq = (size_t) (freq * data.window_size);  /* index of desired frequency */
    const size_t N = data.nr * data.nlat * data.nlon;         /* spatial grid size */
    const size_t T = data.T;                                  /* number of time window segments */
    gsl_matrix_complex *X = gsl_matrix_complex_alloc(3 * N, T);
    gsl_vector *S = gsl_vector_alloc(T);
    gsl_matrix_complex *U = gsl_matrix_complex_alloc(3 * N, T);
    gsl_matrix_complex *V = gsl_matrix_complex_alloc(T, T);
    char buf[2048];

    fprintf(stderr, "main: building matrix X (%zu-by-%zu) for frequency %.2f [cpd]...",
            X->size1, X->size2, freq);
    gettimeofday(&tv0, NULL);
    build_X(ifreq, &data, X);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

    fprintf(stderr, "main: performing SVD of Q for frequency %g [cpd]...", freq);
    gettimeofday(&tv0, NULL);
    status = lapack_complex_svd_thin(X, S, U, V);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds, status = %d)\n", time_diff(tv0, tv1), status);

    sprintf(buf, "%s_%zu", PCA3D_STAGE3A_SVAL_TXT, ifreq);
    fprintf(stderr, "main: writing singular values for frequency %g [cpd] in text format to %s...", freq, buf);
    pca3d_write_S(buf, 0, 0, freq, window_size, window_shift, S);
    fprintf(stderr, "done\n");

    sprintf(buf, "%s_%zu", PCA3D_STAGE3A_U, ifreq);
    fprintf(stderr, "main: writing U matrix for frequency %g [cpd] in binary format to %s...", freq, buf);
    pca3d_write_matrix_complex(buf, U);
    fprintf(stderr, "done\n");

    gsl_matrix_complex_free(X);
    gsl_matrix_complex_free(U);
    gsl_matrix_complex_free(V);
    gsl_vector_free(S);
  }

  free(data.t);
  free(data.r);
  free(data.glat);
  free(data.glon);
  free(data.Jr);
  free(data.Jt);
  free(data.Jp);
  free(data.Qp);
  gsl_vector_free(data.window);

  return 0;
}
