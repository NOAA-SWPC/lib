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

static int calc_fnm(const size_t nmax, const double * g, msynth_grid_workspace * w);

/*
msynth_grid_alloc()
  Allocate msynth grid workspace

Inputs: nphi - number of equally spaced longitude nodes from [0,2pi]
        nmin - minimum spherical harmonic degree for evaluation
        nmax - maximum spherical harmonic degree for evaluation
        g    - internal Gauss coefficients (static, no SV or SA)

Return: pointer to new workspace

Notes:
1) phi_k = 2*pi*k/nphi, k \in [0,1,...,nphi-1]
*/

msynth_grid_workspace *
msynth_grid_alloc(const size_t nphi, const size_t nmin, const size_t nmax,
                  const double * g)
{
  msynth_grid_workspace *w;

  w = calloc(1, sizeof(msynth_grid_workspace));
  if (!w)
    return 0;

  w->nphi = nphi;
  w->nmin = nmin;
  w->nmax = nmax;
  w->msynth_workspace_p = msynth_alloc2(nmax, 1, 1, NULL);

  w->fnm = malloc(w->msynth_workspace_p->nnm * sizeof(complex double));
  w->fnm_mirror = malloc(w->msynth_workspace_p->nnm * sizeof(complex double));
  w->rterm = malloc((w->nmax + 1) * sizeof(double));

  w->fft_r = fftw_malloc(2 * (nphi / 2 + 1) * sizeof(double));
  w->fft_t = fftw_malloc(2 * (nphi / 2 + 1) * sizeof(double));
  w->fft_p = fftw_malloc(2 * (nphi / 2 + 1) * sizeof(double));

  {
    complex double *v = (complex double *) w->fft_r;
    w->fftw_p = fftw_plan_dft_c2r_1d((int) nphi, v, w->fft_r, FFTW_ESTIMATE);
  }

  /* calculate complex Gauss coefficients */
  calc_fnm(w->nmax, g, w);

  return w;
}

void
msynth_grid_free(msynth_grid_workspace *w)
{
  if (w->rterm)
    free(w->rterm);

  if (w->fnm)
    free(w->fnm);

  if (w->fnm_mirror)
    free(w->fnm_mirror);

  if (w->fft_r)
    fftw_free(w->fft_r);

  if (w->fft_t)
    fftw_free(w->fft_t);

  if (w->fft_p)
    fftw_free(w->fft_p);

  if (w->fftw_p)
    fftw_destroy_plan(w->fftw_p);

  if (w->msynth_workspace_p)
    msynth_free(w->msynth_workspace_p);

  fftw_cleanup();

  free(w);
}

/*
msynth_grid_init_r()
  This function should be called each time a new radius is
used, prior to calling msynth_grid_calc(). This function will
initialize the array

w->rterm[n] = (a/r)^{n+2}

needed for the grid calculation inner loop

Inputs: r - radius (km)
        w - workspace
*/

int
msynth_grid_init_r(const double r, msynth_grid_workspace * w)
{
  const msynth_workspace *msynth_p = w->msynth_workspace_p;
  const double ratio = msynth_p->R / r; /* (a/r) */
  double rterm = ratio; /* init: (a/r) */
  size_t n;

  for (n = 0; n <= w->nmax; ++n)
    {
      /* (a/r)^{n+2} */
      rterm *= ratio;
      w->rterm[n] = rterm;
    }

  return GSL_SUCCESS;
}

/*
msynth_grid_init_theta()
  This function should be called each time a new colatitude is
used, prior to calling msynth_grid_calc(). This function will
calculate the associated Legendre functions needed for the grid
calculations

Inputs: theta - colatitude (radians)
        w     - workspace
*/

int
msynth_grid_init_theta(const double theta, msynth_grid_workspace * w)
{
  int s;
  msynth_workspace *msynth_p = w->msynth_workspace_p;

  /* compute associated Legendre functions */
  s = gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, w->nmax, cos(theta),
                                      msynth_p->Plm, msynth_p->dPlm);

  return s;
}

/*
msynth_grid_calc()
  Synthesize magnetic field values at all longitude nodes for
a given r, theta. This function can be called in two modes:

(1) Compute points (r, theta, phi_k), k = 0,...,nphi - 1      (mirror = 0)

(2) Compute points (r, pi - theta, phi_k), k = 0,...,nphi - 1 (mirror = 1)

This allows the calling function to only compute the ALF functions once
for each theta in [0,pi/2], and then automatically compute the terms for
the point (pi - theta) using the identity:

P_{nm}(cos(pi - theta)) = (-1)^{n + m} P_{nm}(cos(theta))

Inputs: mirror - 0 if we are computing the point theta
                 1 if we are computing the mirror point, pi - theta
                 Since the Plm[] and dPlm[] are built for theta,
                 if we compute the mirror point we have to multiply
                 by the phase factor (-1)^{n+m}
        r      - radius (km)
        theta  - colatitude (radians)
        B      - (output) nphi-by-3 matrix
                 B(k,1) = B_x(r, theta, phi_k) (nT)
                 B(k,2) = B_y(r, theta, phi_k) (nT)
                 B(k,3) = B_z(r, theta, phi_k) (nT)
                 where: phi_k = 2*pi*k/nphi
        w      - workspace

Return: success or error

Notes:
1) The functions:

msynth_grid_init_r()
msynth_grid_init_theta()

must be called prior to this function to initialize various
arrays needed for the grid calculation
*/

int
msynth_grid_calc(const int mirror, const double r, const double theta, gsl_matrix *B,
                 msynth_grid_workspace *w)
{
  int s = 0;
  const msynth_workspace *msynth_p = w->msynth_workspace_p;
  const size_t nphi = w->nphi;
  const size_t nmin = w->nmin;
  const size_t nmax = w->nmax;
  const double sint = sin(theta);
  const complex double *fnm = mirror ? w->fnm_mirror : w->fnm; /* f_{nm} or (-1)^{n+m} f_{nm} */
  const double mirror_phase = mirror ? -1.0 : 1.0;
  fftw_complex *r_coeff = (fftw_complex *) w->fft_r;
  fftw_complex *t_coeff = (fftw_complex *) w->fft_t;
  fftw_complex *p_coeff = (fftw_complex *) w->fft_p;
  size_t mmax = w->nmax;
  int m;
  size_t k;

  (void) r;

  assert(nphi > 2 * nmax);

  /* initialize fft arrays to 0 */
  for (k = 0; k < nphi / 2 + 1; ++k)
    {
      r_coeff[k] = 0.0 + I * 0.0;
      t_coeff[k] = 0.0 + I * 0.0;
      p_coeff[k] = 0.0 + I * 0.0;
    }

  for (m = 0; m <= (int) mmax; ++m)
    {
      size_t n;
      size_t n0 = GSL_MAX(nmin, (size_t) m);
      complex double sr = 0.0;
      complex double st = 0.0;
      complex double sp = 0.0;

      for (n = n0; n <= nmax; ++n)
        {
          size_t cidx = msynth_nmidx(n, m, msynth_p);
          size_t aidx = gsl_sf_legendre_array_index(n, m);

          sr += (n + 1.0) * w->rterm[n] * fnm[cidx] * msynth_p->Plm[aidx];
          st -= w->rterm[n] * fnm[cidx] * msynth_p->dPlm[aidx];
          sp -= w->rterm[n] * fnm[cidx] * msynth_p->Plm[aidx];
        }

      /*
       * We need an additional factor of -1 (mirror_phase) in 'st' for the mirror
       * case to account for the d/dtheta P_{nm}(cos(pi-theta)) introducing
       * an extra factor of -1
       */
      r_coeff[m] = sr;
      t_coeff[m] = st * mirror_phase;
      p_coeff[m] = I * m / sint * sp;
    }

  /* perform inverse FFTs */
  fftw_execute_dft_c2r(w->fftw_p, r_coeff, w->fft_r);
  fftw_execute_dft_c2r(w->fftw_p, t_coeff, w->fft_t);
  fftw_execute_dft_c2r(w->fftw_p, p_coeff, w->fft_p);

  for (k = 0; k < nphi; ++k)
    {
      double Br = w->fft_r[k];
      double Bt = w->fft_t[k];
      double Bp = w->fft_p[k];

      gsl_matrix_set(B, k, 0, -Bt);
      gsl_matrix_set(B, k, 1,  Bp);
      gsl_matrix_set(B, k, 2, -Br);
    }

  return s;
}

/*
calc_fnm()
  Compute complex Gauss coefficients

Inputs: nmax - maximum spherical harmonic degree
        g    - real Gauss coefficients, size nnm
        w    - msynth grid workspace

Notes:
1) fnm coefficients are stored in w->fnm on output
2) (-1)^{n+m} * fnm coefficients are stored in w->fnm_mirror on output
*/

static int
calc_fnm(const size_t nmax, const double * g, msynth_grid_workspace * w)
{
  const msynth_workspace * msynth_p = w->msynth_workspace_p;
  size_t n;
  int phasen, /* (-1)^n */
      phasem; /* (-1)^m */

  phasen = -1;
  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;
      int m;
      size_t gidx, hidx;

      /* m = 0 case */
      gidx = msynth_nmidx(n, 0, msynth_p);
      w->fnm[gidx] = g[gidx];
      w->fnm_mirror[gidx] = phasen * g[gidx];

      /* m > 0 case */
      phasem = -1;
      for (m = 1; m <= ni; ++m)
        {
          gidx = msynth_nmidx(n, m, msynth_p);
          hidx = msynth_nmidx(n, -m, msynth_p);
          w->fnm[gidx] = 0.5 * (g[gidx] - I * g[hidx]);
          w->fnm_mirror[gidx] = phasen * phasem * w->fnm[gidx];
          phasem = -phasem;
        }

      /* m < 0 case */
      phasem = -1;
      for (m = -1; m >= -ni; --m)
        {
          int mabs = abs(m);
          gidx = msynth_nmidx(n, mabs, msynth_p);
          hidx = msynth_nmidx(n, -mabs, msynth_p);
          w->fnm[hidx] = 0.5 * (g[gidx] - I * g[hidx]);
          w->fnm_mirror[hidx] = phasen * phasem * w->fnm[hidx];
          phasem = -phasem;
        }

      phasen = -phasen;
    }

  return GSL_SUCCESS;
}
