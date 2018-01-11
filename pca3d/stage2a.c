/*
 * stage2.c
 *
 * Perform PCA (frequency domain) directly on the 3D current grids from TIEGCM
 *
 * 1. Read 3D current grids for some time interval
 * 2. Divide total time interval into T smaller segments; perform windowed FFT
 *    of each time series segment for all grid points
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include <assert.h>

#include <fftw3.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <common/common.h>

#include "pca3d.h"
#include "tiegcm3d.h"
#include "window.h"

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

static size_t
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

/*
do_transforms()
  Perform FFTs on each time window segment for all grid points in TIEGCM grid

Inputs: data_file    - output data file for FFT grids
        fs           - sampling frequency (samples/day)
        window_size  - number of days in each time window
        window_shift - how many days to advance window
        data         - TIEGCM data
*/

static int
do_transforms(const char *data_file, const double fs, const double window_size, const double window_shift,
              const tiegcm3d_data * data)
{
  int s = 0;
  const size_t nt = data->nt;
  size_t ir, ilat, ilon;
  size_t nwindow = (size_t) (window_size * fs);   /* optimal number of samples per window */
  size_t nforward = (size_t) (window_shift * fs); /* number of samples to slide forward */
  size_t nfreq = nwindow / 2 + 1;
  size_t T;                                       /* number of time windows */
  gsl_vector *work = gsl_vector_alloc(nwindow);
  gsl_vector *window = gsl_vector_alloc(nwindow);
  fftw_complex *fft_out = fftw_malloc(sizeof(fftw_complex) * nfreq);
  fftw_plan plan;
  struct timeval tv0, tv1;
  gsl_complex *Qp; /* FT'd Jp grid, nfreq-by-nr-by-nlat-by-nlon */
  tiegcm3d_fft_data fft_data;

  Qp = malloc(nfreq * data->nr * data->nlat * data->nlon * sizeof(gsl_complex));

  T = count_windows(nt, fs, window_size, window_shift);

  fprintf(stderr, "do_transforms: samples per window   = %zu\n", nwindow);
  fprintf(stderr, "do_transforms: sample slide forward = %zu\n", nforward);
  fprintf(stderr, "do_transforms: number of freqs      = %zu\n", nfreq);
  fprintf(stderr, "do_transforms: time samples         = %zu [hourly]\n", nt);
  fprintf(stderr, "do_transforms: window segments (T)  = %zu\n", T);

  /* compute window function */
  /*apply_ps1(NULL, window);*/
  apply_hamming(NULL, window);

  fft_data.nt = data->nt;
  fft_data.nfreq = nfreq;
  fft_data.nr = data->nr;
  fft_data.nlat = data->nlat;
  fft_data.nlon = data->nlon;
  fft_data.T = T;
  fft_data.fs = fs;
  fft_data.window_size = window_size;
  fft_data.window_shift = window_shift;
  fft_data.nwindow = nwindow;
  fft_data.t = data->t;
  fft_data.r = data->r;
  fft_data.glat = data->glat;
  fft_data.glon = data->glon;
  fft_data.window = window;
  fft_data.Jr = data->Jr;
  fft_data.Jt = data->Jt;
  fft_data.Jp = data->Jp;
  fft_data.Qp = Qp;

  plan = fftw_plan_dft_r2c_1d(nwindow, work->data, fft_out, FFTW_ESTIMATE);

  fprintf(stderr, "do_transforms: computing FFTs of windowed data...");
  gettimeofday(&tv0, NULL);

  for (ir = 0; ir < data->nr; ++ir)
    {
      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          for (ilon = 0; ilon < data->nlon; ++ilon)
            {
              size_t start_idx = 0; /* starting time index */
              size_t t;

              for (t = 0; t < T; ++t)
                {
                  size_t end_idx = GSL_MIN(start_idx + nwindow - 1, nt - 1);
                  size_t n = end_idx - start_idx + 1; /* size of actual window */
                  double sqrtn = sqrt((double) n);
                  size_t it, ifreq;

                  assert(start_idx < end_idx);

                  if (n < nwindow)
                    {
                      /* could happen at the end of the time series; zero pad input buffer */
                      gsl_vector_set_zero(work);
                    }

                  /* copy current time window into work array */
                  for (it = start_idx; it <= end_idx; ++it)
                    {
                      size_t idx = TIEGCM3D_IDX(it, ir, ilat, ilon, data);
                      double Jp = data->Jp[idx];

                      gsl_vector_set(work, it - start_idx, Jp);
                    }

                  /* apply window function */
                  gsl_vector_mul(work, window);

#if 0
                  for (it = 0; it < n; ++it)
                    {
                      size_t idx = TIEGCM3D_IDX(it + start_idx, ir, ilat, ilon, data);
                      printf("%ld %.12e %.12e\n",
                             data->t[it + start_idx],
                             data->Jp[idx],
                             gsl_vector_get(work, it));
                    }
                  printf("\n\n");
#endif

                  /* compute FFT of this windowed data */
                  fftw_execute(plan);

                  /* store FFT result in Q grids */
                  for (ifreq = 0; ifreq < nfreq; ++ifreq)
                    {
                      size_t idx = TIEGCM3D_FREQIDX(ifreq, ir, ilat, ilon, data, nfreq);
                      gsl_complex z = gsl_complex_rect(fft_out[ifreq][0] / sqrtn, fft_out[ifreq][1] / sqrtn);
                      double freq = ifreq * (fs / n);
                      double period = 1.0 / freq;
                      double power = gsl_complex_abs2(z);
                      Qp[idx] = z;

#if 1
                      printf("%f %f %.12e\n", freq, period, power);
#endif
                    }
#if 1
                  printf("\n\n");
#endif

                  start_idx += nforward;
                }
#if 1
              exit(1);
#endif
            }
        }
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "do_transforms: writing FFT grids to %s...", data_file);
  pca3d_write_fft_data(data_file, &fft_data);
  fprintf(stderr, "done\n");

  free(Qp);
  fftw_free(fft_out);
  gsl_vector_free(work);
  gsl_vector_free(window);
  fftw_destroy_plan(plan);

  fftw_cleanup();

  return s;
}

int
main(int argc, char *argv[])
{
  const char *fft_data_file = PCA3D_STAGE2A_FFT_DATA;
  const double fs = 24.0;    /* sample frequency (samples/day) */
  char *infile = NULL;
  double window_size = 2.0;  /* number of days in each time segment */
  double window_shift = 1.0; /* number of days to shift forward in time */
  struct timeval tv0, tv1;
  tiegcm3d_data *data;

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
            fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file> [-t window_size (days)]\n", argv[0]);
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file> [-t window_size (days)]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: input file          = %s\n", infile);
  fprintf(stderr, "main: sample frequency    = %g [samples/day]\n", fs);
  fprintf(stderr, "main: window size         = %g [days]\n", window_size);
  fprintf(stderr, "main: window shift        = %g [days]\n", window_shift);

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = tiegcm3d_read(infile, NULL);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->nt,
          time_diff(tv0, tv1));

  do_transforms(fft_data_file, fs, window_size, window_shift, data);

  return 0;
}
