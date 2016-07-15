/*
 * green.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_legendre.h>

#include "green.h"

/*
green_alloc()
  Allocate Green's function workspace

Inputs: nmax - maximum spherical harmonic degree
        mmax - maximum spherical harmonic order
        R    - reference radius (km)
*/

green_workspace *
green_alloc(const size_t nmax, const size_t mmax, const double R)
{
  green_workspace *w;
  size_t plm_array_size = gsl_sf_legendre_array_n(nmax);

  if (mmax > nmax)
    {
      fprintf(stderr, "green_alloc: error: mmax > nmax\n");
      return 0;
    }

  w = calloc(1, sizeof(green_workspace));
  if (!w)
    return 0;

  /* total number of SH coefficients */
  w->nnm = mmax * (mmax + 2) + (nmax - mmax) * (2*mmax + 1);

  w->nmax = nmax;
  w->mmax = mmax;
  w->R = R;

  w->cosmphi = malloc((mmax + 1) * sizeof(double));
  w->sinmphi = malloc((mmax + 1) * sizeof(double));

  w->Plm = malloc(plm_array_size * sizeof(double));
  w->dPlm = malloc(plm_array_size * sizeof(double));
  if (!w->Plm || !w->dPlm)
    {
      green_free(w);
      return 0;
    }

  return w;
}

void
green_free(green_workspace *w)
{
  if (w->cosmphi)
    free(w->cosmphi);

  if (w->sinmphi)
    free(w->sinmphi);

  if (w->Plm)
    free(w->Plm);

  if (w->dPlm)
    free(w->dPlm);

  free(w);
}

/*
green_calc_int()
  Compute Green's functions for X,Y,Z spherical harmonic expansion. These
are simply the basis functions multiplying the g_{nm} and h_{nm} coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        X     - (output) array of X Green's functions, size nnm
        Y     - (output) array of Y Green's functions, size nnm
        Z     - (output) array of Z Green's functions, size nnm
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Plm
w->dPlm
w->sinmphi
w->cosmphi
*/

int
green_calc_int(const double r, const double theta, const double phi,
               double *X, double *Y, double *Z, green_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax;
  const size_t mmax = w->mmax;
  size_t n;
  int m;
  const double sint = sin(theta);
  const double cost = cos(theta);
  double ratio = w->R / r;
  double term = ratio * ratio;     /* (a/r)^{n+2} */

  /* precompute cos(m phi) and sin(m phi) */
  for (m = 0; m <= (int) mmax; ++m)
    {
      w->cosmphi[m] = cos(m * phi);
      w->sinmphi[m] = sin(m * phi);
    }

  /* compute associated legendres */
  gsl_sf_legendre_deriv_array(GSL_SF_LEGENDRE_SCHMIDT,
                              nmax, cost, w->Plm, w->dPlm);

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(mmax, n);

      /* (a/r)^{n+2} */
      term *= ratio;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t cidx = green_nmidx(n, m, w);
          size_t aidx = gsl_sf_legendre_array_index(n, mabs);

          if (m < 0)
            {
              /* h_{nm} */
              X[cidx] = term * w->sinmphi[mabs] * w->dPlm[aidx] * (-sint);
              Y[cidx] = -term / sint * mabs * w->cosmphi[mabs] * w->Plm[aidx];
              Z[cidx] = -(n + 1.0) * term * w->sinmphi[mabs] * w->Plm[aidx];
            }
          else
            {
              /* g_{nm} */
              X[cidx] = term * w->cosmphi[mabs] * w->dPlm[aidx] * (-sint);
              Y[cidx] = term / sint * mabs * w->sinmphi[mabs] * w->Plm[aidx];
              Z[cidx] = -(n + 1.0) * term * w->cosmphi[mabs] * w->Plm[aidx];
            }
        }
    }

  return s;
} /* green_calc_int() */

/*
green_calc_ext()
  Compute Green's functions for X,Y,Z spherical harmonic expansion due to
external current source. These are simply the basis functions multiplying
the k_{nm} and q_{nm} coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        X     - (output) array of X Green's functions, size nnm
        Y     - (output) array of Y Green's functions, size nnm
        Z     - (output) array of Z Green's functions, size nnm
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Plm
w->dPlm
w->sinmphi
w->cosmphi
*/

int
green_calc_ext(const double r, const double theta, const double phi,
               double *X, double *Y, double *Z, green_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax;
  const size_t mmax = w->mmax;
  size_t n;
  int m;
  const double sint = sin(theta);
  const double cost = cos(theta);
  double ratio = r / w->R;
  double term = 1.0;     /* (r/a)^{n-1} */

  /* precompute cos(m phi) and sin(m phi) */
  for (m = 0; m <= (int) mmax; ++m)
    {
      w->cosmphi[m] = cos(m * phi);
      w->sinmphi[m] = sin(m * phi);
    }

  /* compute associated legendres */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax, cost,
                                  w->Plm, w->dPlm);

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(mmax, n);

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t cidx = green_nmidx(n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);

          if (m < 0)
            {
              /* q_{nm} */
              X[cidx] = term * w->sinmphi[mabs] * w->dPlm[pidx];
              Y[cidx] = -term / sint * mabs * w->cosmphi[mabs] * w->Plm[pidx];
              Z[cidx] = (double) n * term * w->sinmphi[mabs] * w->Plm[pidx];
            }
          else
            {
              /* k_{nm} */
              X[cidx] = term * w->cosmphi[mabs] * w->dPlm[pidx];
              Y[cidx] = term / sint * mabs * w->sinmphi[mabs] * w->Plm[pidx];
              Z[cidx] = (double) n * term * w->cosmphi[mabs] * w->Plm[pidx];
            }
        }

      /* (r/a)^{n-1} */
      term *= ratio;
    }

  return s;
} /* green_calc_ext() */

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

size_t
green_nmidx(const size_t n, const int m, const green_workspace *w)
{
  const size_t mmax = w->mmax;
  size_t nmidx;

  if (n == 0)
    {
      fprintf(stderr, "green_nmidx: error: n = 0\n");
      return 0;
    }
  else if (abs(m) > (int) mmax)
    {
      fprintf(stderr, "green_nmidx: error: m = %d [mmax = %zu]\n", m, mmax);
      return 0;
    }

  if (n <= mmax)
    {
      size_t base = n * n; /* index of block for this n */
      int offset = m + n;  /* offset within block for this m */

      /* subtract 1 to exclude (0,0) coefficient */
      nmidx = base + offset - 1;
    }
  else
    {
      size_t base1 = (mmax + 1) * (mmax + 1);
      size_t base2 = (n - mmax - 1) * (2 * mmax + 1);
      int offset = m + (int)mmax;

      /* subtract 1 to exclude (0,0) coefficient */
      nmidx = base1 + base2 + offset - 1;
    }

  return nmidx;
}

size_t
green_nnm(const green_workspace *w)
{
  return w->nnm;
}

int
green_print_spectrum(const char *filename, const gsl_vector *c,
                     const green_workspace *w)
{
  const size_t nmax = w->nmax;
  size_t n;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_spectrum: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  n = 1;
  fprintf(fp, "# Field %zu: spherical harmonic degree n\n", n++);
  fprintf(fp, "# Field %zu: power (nT^2)\n", n++);

  for (n = 1; n <= nmax; ++n)
    {
      double sum = 0.0;
      int M = (int) GSL_MIN(w->mmax, n);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = green_nmidx(n, m, w);
          double gnm = gsl_vector_get(c, cidx);

          sum += gnm * gnm;
        }

      sum *= (n + 1.0);

      fprintf(fp, "%zu %.12e\n", n, sum);
    }

  fclose(fp);

  return 0;
}

int
green_print_spectrum_azim(const char *filename, const gsl_vector * c,
                          const green_workspace *w)
{
  size_t n;
  FILE *fp;
  
  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "green_print_spectrum_azim: cannot open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  n = 1;
  fprintf(fp, "# Field %zu: spherical harmonic degree n\n", n++);
  fprintf(fp, "# Field %zu: spherical harmonic order m\n", n++);
  fprintf(fp, "# Field %zu: m/n\n", n++);
  fprintf(fp, "# Field %zu: n + 0.25*(m/n)\n", n++);
  fprintf(fp, "# Field %zu: gnm^2 (nT^2)\n", n++);

  for (n = 1; n <= w->nmax; ++n)
    {
      int M = (int) GSL_MIN(w->mmax, n);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = green_nmidx(n, m, w);
          double gnm = gsl_vector_get(c, cidx);
          double ratio = (double)m / (double)n;

          fprintf(fp, "%4zu %4d %8.4f %8.4f %.12e\n",
                  n,
                  m,
                  ratio,
                  n + 0.25 * ratio,
                  gnm * gnm);
        }
    }

  fclose(fp);

  return 0;
}
