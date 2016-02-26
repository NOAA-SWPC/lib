/*
 * green.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>

#include "green.h"

green_workspace *
green_alloc(const size_t nmax)
{
  green_workspace *w;
  size_t plm_array_size = gsl_sf_legendre_array_n(nmax);
  size_t nnm_tot;

  w = calloc(1, sizeof(green_workspace));
  if (!w)
    return 0;

  nnm_tot = (nmax + 1) * (nmax + 1);

  w->nmax = nmax;
  w->nnm = nnm_tot - 1; /* exclude (0,0) coefficient */
  w->R = 6371.2;

  w->cosmphi = malloc((nmax + 1) * sizeof(double));
  w->sinmphi = malloc((nmax + 1) * sizeof(double));

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
  size_t n;
  int m;
  const double sint = sin(theta);
  const double cost = cos(theta);
  double ratio = w->R / r;
  double term = ratio * ratio;     /* (a/r)^{n+2} */

  /* precompute cos(m phi) and sin(m phi) */
  for (n = 0; n <= nmax; ++n)
    {
      w->cosmphi[n] = cos(n * phi);
      w->sinmphi[n] = sin(n * phi);
    }

  /* compute associated legendres */
  gsl_sf_legendre_deriv_array(GSL_SF_LEGENDRE_SCHMIDT,
                              nmax, cost, w->Plm, w->dPlm);

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;

      /* (a/r)^{n+2} */
      term *= ratio;

      for (m = -ni; m <= ni; ++m)
        {
          int mabs = abs(m);
          size_t cidx = green_nmidx(n, m);
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
green_calc_ext(const double r, const double theta, const double phi,
               double *X, double *Y, double *Z, green_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax;
  size_t n;
  int m;
  const double sint = sin(theta);
  const double cost = cos(theta);
  double ratio = r / w->R;
  double term = 1.0;     /* (r/a)^{n-1} */

  /* precompute cos(m phi) and sin(m phi) */
  for (n = 0; n <= nmax; ++n)
    {
      w->cosmphi[n] = cos(n * phi);
      w->sinmphi[n] = sin(n * phi);
    }

  /* compute associated legendres */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax, cost,
                                  w->Plm, w->dPlm);

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          int mabs = abs(m);
          size_t cidx = green_nmidx(n, m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);

          if (m < 0)
            {
              /* h_{nm} */
              X[cidx] = term * w->sinmphi[mabs] * w->dPlm[pidx];
              Y[cidx] = -term / sint * mabs * w->cosmphi[mabs] * w->Plm[pidx];
              Z[cidx] = (double) n * term * w->sinmphi[mabs] * w->Plm[pidx];
            }
          else
            {
              /* g_{nm} */
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
  This function returns a unique index in [0,w->p-1] corresponding
to a given (l,m) pair. The array will look like:

[(1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) ...]

(the (0,0) coefficient is not solved for)

Inputs: n - SH degree (> 0)
        m - SH order (-l <= m <= l)

Return: index in [0,nnm-1]
*/

size_t
green_nmidx(const size_t n, const int m)
{
  size_t base = n * n; /* index of block for this n */
  int offset = m + n;  /* offset within block for this m */
  size_t nmidx;

  if (n == 0)
    {
      fprintf(stderr, "green_nmidx: error: n = 0\n");
      return 0;
    }

  nmidx = base + offset;

  /* subtract 1 to exclude (0,0) coefficient */
  return nmidx - 1;
}
