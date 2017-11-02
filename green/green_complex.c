/*
 * green_complex.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <complex.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_legendre.h>

#include "green.h"
#include "green_complex.h"

static int green_complex_Ynm(const double theta, const double phi, green_complex_workspace *w);

/*
green_complex_alloc()
  Allocate Green's function workspace

Inputs: nmax - maximum spherical harmonic degree
        mmax - maximum spherical harmonic order
        R    - reference radius (km)
*/

green_complex_workspace *
green_complex_alloc(const size_t nmax, const size_t mmax, const double R)
{
  green_complex_workspace *w;
  size_t plm_array_size = gsl_sf_legendre_array_n(nmax);

  if (mmax > nmax)
    {
      fprintf(stderr, "green_complex_alloc: error: mmax > nmax\n");
      return 0;
    }

  w = calloc(1, sizeof(green_complex_workspace));
  if (!w)
    return 0;

  /* total number of SH coefficients */
  w->nnm = mmax * (mmax + 2) + (nmax - mmax) * (2*mmax + 1);

  w->nmax = nmax;
  w->mmax = mmax;
  w->R = R;

  w->Pnm = malloc(plm_array_size * sizeof(double));
  w->dPnm = malloc(plm_array_size * sizeof(double));
  if (!w->Pnm || !w->dPnm)
    {
      green_complex_free(w);
      return 0;
    }

  w->Ynm = malloc(plm_array_size * sizeof(complex double));
  w->dYnm = malloc(plm_array_size * sizeof(complex double));
  if (!w->Ynm || !w->dYnm)
    {
      green_complex_free(w);
      return 0;
    }

  return w;
}

void
green_complex_free(green_complex_workspace *w)
{
  if (w->Pnm)
    free(w->Pnm);

  if (w->dPnm)
    free(w->dPnm);

  if (w->Ynm)
    free(w->Ynm);

  if (w->dYnm)
    free(w->dYnm);

  free(w);
}

/*
green_complex_calc_int()
  Compute Green's functions for X,Y,Z spherical harmonic expansion. These
are simply the basis functions multiplying the g_{nm} coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        X     - (output) array of X Green's functions, size nnm
        Y     - (output) array of Y Green's functions, size nnm
        Z     - (output) array of Z Green's functions, size nnm
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Pnm
w->dPnm
w->Ynm
w->dYnm
*/

int
green_complex_calc_int(const double r, const double theta, const double phi,
                       complex double *X, complex double *Y, complex double *Z,
                       green_complex_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax;
  const size_t mmax = w->mmax;
  const complex double invisint = I / sin(theta);
  const double ratio = w->R / r;
  double rterm = ratio * ratio; /* (R/r)^{n+2} */
  size_t n;

  /* compute Y_{nm} and d/dtheta Y_{nm} */
  green_complex_Ynm(theta, phi, w);

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;

      /* (R/r)^(n+2) */
      rterm *= ratio;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = green_idx(n, m, mmax);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          complex double Ynm, dYnm;

          if (m >= 0)
            {
              Ynm = w->Ynm[pidx];
              dYnm = w->dYnm[pidx];
            }
          else
            {
              Ynm = conj(w->Ynm[pidx]);
              dYnm = conj(w->dYnm[pidx]);
            }

          X[nmidx] = rterm * dYnm;
          Y[nmidx] = -rterm * m * invisint * Ynm;
          Z[nmidx] = -rterm * (n + 1.0) * Ynm;
        }
    }

  return s;
}

/*
green_complex_calc_ext()
  Compute Green's functions for X,Y,Z spherical harmonic expansion due to
external current source. These are simply the basis functions multiplying
the k_{nm} coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        X     - (output) array of X Green's functions, size nnm
        Y     - (output) array of Y Green's functions, size nnm
        Z     - (output) array of Z Green's functions, size nnm
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Pnm
w->dPnm
w->Ynm
w->dYnm
*/

int
green_complex_calc_ext(const double r, const double theta, const double phi,
                       complex double *X, complex double *Y, complex double *Z,
                       green_complex_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax;
  const size_t mmax = w->mmax;
  const complex double invisint = I / sin(theta);
  const double ratio = r / w->R;
  double rterm = 1.0; /* (r/R)^{n-1} */
  size_t n;

  /* compute Y_{nm} and d/dtheta Y_{nm} */
  green_complex_Ynm(theta, phi, w);

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = green_idx(n, m, mmax);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          complex double Ynm, dYnm;

          if (m >= 0)
            {
              Ynm = w->Ynm[pidx];
              dYnm = w->dYnm[pidx];
            }
          else
            {
              Ynm = conj(w->Ynm[pidx]);
              dYnm = conj(w->dYnm[pidx]);
            }

          X[nmidx] = rterm * dYnm;
          Y[nmidx] = -rterm * m * invisint * Ynm;
          Z[nmidx] = rterm * n * Ynm;
        }

      /* (r/R)^{n-1} */
      rterm *= ratio;
    }

  return s;
}

/*
green_nmidx()
  This function returns a unique index in [0,p-1] corresponding
to a given (n,m) pair. The array will look like:

[(1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) ...]

(the (0,0) coefficient is not solved for)

Things get a little more tricky when mmax != nmax, so the
base indices of each 'n' block are precomputed in _alloc and
stored for easy reference here. The offset of (n,m) in a given
'n' block is:

offset = m + min(n,mmax)

which defaults to the standard m + n for the case mmax = nmax

Inputs: n - SH degree in [1,nmax]
        m - SH order (-n <= m <= n)
        w - workspace

Return: index in [0,nnm-1]
*/

inline size_t
green_complex_nmidx(const size_t n, const int m, const green_complex_workspace *w)
{
  return green_idx(n, m, w->mmax);
}

size_t
green_complex_nnm(const green_complex_workspace *w)
{
  return w->nnm;
}

/*
green_complex_Ynm()
  Compute Y_{nm}(theta,phi) = P_{nm}(cos(theta)) * exp(i m phi) as
well as d/dtheta Y_{nm}(theta,phi)

Inputs: theta - colatitude (radians)
        phi   - longitude (radians)
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Plm
w->dPlm
w->Ylm
w->dYlm
*/

static int
green_complex_Ynm(const double theta, const double phi, green_complex_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax;
  const size_t mmax = w->mmax;
  size_t n;
  int m;

  /* compute associated Legendres */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT,
                                  nmax, cos(theta), w->Pnm, w->dPnm);

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= (int) mmax; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);

      for (n = GSL_MAX(m, 1); n <= nmax; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;
        }
    }

  return s;
}
