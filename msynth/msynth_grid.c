/*
 * msynth_grid.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_legendre.h>

#include "msynth.h"
#include "msynth_grid.h"

/*
msynth_grid_alloc()
  Allocate msynth grid workspace

Inputs: nphi     - number of equally spaced longitude nodes from [0,2pi]
        msynth_p - msynth workspace containing model coefficients

Return: pointer to new workspace

Notes:
1) phi_k = 2*pi*k/nphi, k \in [0,1,...,nphi-1]
*/

msynth_grid_workspace *
msynth_grid_alloc(const size_t nphi, const msynth_workspace *msynth_p)
{
  msynth_grid_workspace *w;

  w = calloc(1, sizeof(msynth_grid_workspace));
  if (!w)
    return 0;

  w->nphi = nphi;
  w->nmin = msynth_p->eval_nmin;
  w->nmax = msynth_p->eval_nmax;
  w->msynth_workspace_p = msynth_p;

  w->fft_x = fftw_malloc(2 * (nphi / 2 + 1) * sizeof(double));
  w->fft_y = fftw_malloc(2 * (nphi / 2 + 1) * sizeof(double));
  w->fft_z = fftw_malloc(2 * (nphi / 2 + 1) * sizeof(double));

  {
    complex double *v = (complex double *) w->fft_x;
    w->fftw_p = fftw_plan_dft_c2r_1d((int) nphi, v, w->fft_x, FFTW_ESTIMATE);
  }

  return w;
} /* msynth_grid_alloc() */

void
msynth_grid_free(msynth_grid_workspace *w)
{
  if (w->fft_x)
    fftw_free(w->fft_x);

  if (w->fft_y)
    fftw_free(w->fft_y);

  if (w->fft_z)
    fftw_free(w->fft_z);

  if (w->fftw_p)
    fftw_destroy_plan(w->fftw_p);

  fftw_cleanup();

  free(w);
}

/*
msynth_grid_calc()
  Synthesize magnetic field values at all longitude nodes for
a given r,theta

Inputs: r     - radius (km)
        theta - colatitude (radians)
        B     - (output) nphi-by-3 matrix
                B(:,1) = B_x (nT)
                B(:,2) = B_y (nT)
                B(:,3) = B_z (nT)
        w     - workspace

Return: success or error

Notes:
1) w->Plm and w->dPlm must be already initialized with ALFs for cos(theta),
   to enable fast loop calculations on a grid
*/

int
msynth_grid_calc(const double r, const double theta, gsl_matrix *B,
                 msynth_grid_workspace *w)
{
  int s = 0;
  const msynth_workspace *msynth_p = w->msynth_workspace_p;
  const size_t nphi = w->nphi;
  const size_t nmin = w->nmin;
  const size_t nmax = w->nmax;
  const double sint = sin(theta);
  const double ratio = msynth_p->R / r; /* (a/r) */
  const double *g = msynth_p->c;        /* Gauss coefficients */
  fftw_complex *x = (fftw_complex *) w->fft_x;
  fftw_complex *y = (fftw_complex *) w->fft_y;
  fftw_complex *z = (fftw_complex *) w->fft_z;
  size_t mmax = w->nmax;
  int m;
  size_t k;

  assert(nphi > 2 * nmax);

  /* initialize fft arrays to 0 */
  for (k = 0; k < nphi / 2 + 1; ++k)
    {
      x[k] = 0.0 + I * 0.0;
      y[k] = 0.0 + I * 0.0;
      z[k] = 0.0 + I * 0.0;
    }

  for (m = 0; m <= (int) mmax; ++m)
    {
      size_t n;
      size_t n0 = GSL_MAX(nmin, (size_t) m);
      double rterm = gsl_pow_int(ratio, n0 + 1); /* initialize to (a/r)^{n+1} */
      complex double sx = 0.0;
      complex double sy = 0.0;
      complex double sz = 0.0;

      for (n = n0; n <= nmax; ++n)
        {
          size_t cidx = msynth_nmidx(n, m, msynth_p);
          size_t aidx = gsl_sf_legendre_array_index(n, m);
          complex double fnm;

          /* (a/r)^{n+2} */
          rterm *= ratio;

          if (m == 0)
            {
              fnm = g[cidx];
            }
          else /* m > 0 */
            {
              fnm = 0.5 * g[cidx];
              cidx = msynth_nmidx(n, -m, msynth_p);
              fnm -= 0.5 * I * g[cidx];
            }

          sx += rterm * fnm * msynth_p->dPlm[aidx];
          sy -= rterm * fnm * msynth_p->Plm[aidx];
          sz -= (n + 1.0) * rterm * fnm * msynth_p->Plm[aidx];
        }

      x[m] = sx;
      y[m] = I * m / sint * sy;
      z[m] = sz;
    }

  /* perform inverse FFTs */
  fftw_execute_dft_c2r(w->fftw_p, x, w->fft_x);
  fftw_execute_dft_c2r(w->fftw_p, y, w->fft_y);
  fftw_execute_dft_c2r(w->fftw_p, z, w->fft_z);

  for (k = 0; k < nphi; ++k)
    {
      gsl_matrix_set(B, k, 0, w->fft_x[k]);
      gsl_matrix_set(B, k, 1, w->fft_y[k]);
      gsl_matrix_set(B, k, 2, w->fft_z[k]);
    }

  return s;
}

/*
msynth_grid_calc2()
  Synthesize magnetic field values at all longitude nodes for
a given r, and the point conjugate to theta across the equator:

B(r, pi - theta, phi_k)

If the Legendre polynomials have already been computed for theta,
then we can avoid computing them again for pi - theta using
the simple identity:

P_n^m(cos(pi-theta)) = (-1)^{n+m} P_n^m(cos(theta))

Inputs: r     - radius (km)
        theta - colatitude (radians)
        B     - (output) nphi-by-3 matrix
                B(:,1) = B_x (nT)
                B(:,2) = B_y (nT)
                B(:,3) = B_z (nT)
        w     - workspace

Return: success or error

Notes:
1) w->Plm and w->dPlm must be already initialized with ALFs for cos(theta);
   they will automatically be adjusted for cos(pi - theta)
*/

int
msynth_grid_calc2(const double r, const double theta, gsl_matrix *B,
                  msynth_grid_workspace *w)
{
  int s = 0;
  const msynth_workspace *msynth_p = w->msynth_workspace_p;
  const size_t nphi = w->nphi;
  const size_t nmin = w->nmin;
  const size_t nmax = w->nmax;
  const double sint = sin(theta); /* sin(pi - theta) = sin(theta) */
  const double ratio = msynth_p->R / r; /* (a/r) */
  const double *g = msynth_p->c;        /* Gauss coefficients */
  fftw_complex *x = (fftw_complex *) w->fft_x;
  fftw_complex *y = (fftw_complex *) w->fft_y;
  fftw_complex *z = (fftw_complex *) w->fft_z;
  size_t mmax = w->nmax;
  int m;
  size_t k;

  assert(nphi > 2 * nmax);

  /* initialize fft arrays to 0 */
  for (k = 0; k < nphi / 2 + 1; ++k)
    {
      x[k] = 0.0 + I * 0.0;
      y[k] = 0.0 + I * 0.0;
      z[k] = 0.0 + I * 0.0;
    }

  for (m = 0; m <= (int) mmax; ++m)
    {
      size_t n;
      size_t n0 = GSL_MAX(nmin, (size_t) m);
      double rterm = gsl_pow_int(ratio, n0 + 1); /* initialize to (a/r)^{n+1} */
      complex double sx = 0.0;
      complex double sy = 0.0;
      complex double sz = 0.0;
      double fac; /* (-1)^{n + m} */

      if ((n0 + m) % 2 == 0)
        fac = 1.0;
      else
        fac = -1.0;

      for (n = n0; n <= nmax; ++n)
        {
          size_t cidx = msynth_nmidx(n, m, msynth_p);
          size_t aidx = gsl_sf_legendre_array_index(n, m);
          complex double fnm;
#if 0
          double fac = gsl_pow_int(-1.0, (int) n + m); /*XXX: FIX!!! */
#endif

          /* (a/r)^{n+2} */
          rterm *= ratio;

          if (m == 0)
            {
              fnm = g[cidx];
            }
          else /* m > 0 */
            {
              fnm = 0.5 * g[cidx];
              cidx = msynth_nmidx(n, -m, msynth_p);
              fnm -= 0.5 * I * g[cidx];
            }

          sx -= rterm * fnm * fac * msynth_p->dPlm[aidx];
          sy -= rterm * fnm * fac * msynth_p->Plm[aidx];
          sz -= (n + 1.0) * rterm * fnm * fac * msynth_p->Plm[aidx];

          fac *= -1.0;
        }

      x[m] = sx;
      y[m] = I * m / sint * sy;
      z[m] = sz;
    }

  /* perform inverse FFTs */
  fftw_execute_dft_c2r(w->fftw_p, x, w->fft_x);
  fftw_execute_dft_c2r(w->fftw_p, y, w->fft_y);
  fftw_execute_dft_c2r(w->fftw_p, z, w->fft_z);

  for (k = 0; k < nphi; ++k)
    {
      gsl_matrix_set(B, k, 0, w->fft_x[k]);
      gsl_matrix_set(B, k, 1, w->fft_y[k]);
      gsl_matrix_set(B, k, 2, w->fft_z[k]);
    }

  return s;
}
