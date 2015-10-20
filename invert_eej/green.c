/*
 * green.c
 *
 * This module contains routines for computing the Green's (basis)
 * functions of spherical harmonic expansions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_test.h>

#include "green.h"

green_workspace *
green_alloc(const size_t nmax)
{
  green_workspace *w;
  size_t nnm_tot, array_size;

  w = calloc(1, sizeof(green_workspace));
  if (!w)
    return 0;

  w->nmax = nmax;

  nnm_tot = (nmax + 1) * (nmax + 1);

  /* exclude the (0,0) coefficient */
  w->nnm = nnm_tot - 1;

  w->dX = malloc(w->nnm * sizeof(double));
  w->dY = malloc(w->nnm * sizeof(double));
  w->dZ = malloc(w->nnm * sizeof(double));
  if (!w->dX || !w->dY || !w->dZ)
    {
      green_free(w);
      fprintf(stderr, "green_alloc: failed to alloc X,Y,Z arrays\n");
      return 0;
    }

  w->dX_ext = malloc(w->nnm * sizeof(double));
  w->dY_ext = malloc(w->nnm * sizeof(double));
  w->dZ_ext = malloc(w->nnm * sizeof(double));
  if (!w->dX_ext || !w->dY_ext || !w->dZ_ext)
    {
      green_free(w);
      fprintf(stderr, "green_alloc: failed to alloc X_ext,Y_ext,Z_ext arrays\n");
      return 0;
    }

  w->cosmphi = malloc((nmax + 1) * sizeof(double));
  w->sinmphi = malloc((nmax + 1) * sizeof(double));
  if (!w->cosmphi || !w->sinmphi)
    {
      green_free(w);
      fprintf(stderr, "green_alloc: failed to alloc sin/cos arrays\n");
      return 0;
    }

  array_size = gsl_sf_legendre_array_n(nmax);
  w->Plm = malloc(array_size * sizeof(double));
  w->dPlm = malloc(array_size * sizeof(double));
  if (!w->Plm || !w->dPlm)
    {
      green_free(w);
      fprintf(stderr, "green_alloc: failed to alloc legendre arrays\n");
      return 0;
    }

  return w;
} /* green_alloc() */

void
green_free(green_workspace *w)
{
  if (w->dX)
    free(w->dX);

  if (w->dY)
    free(w->dY);

  if (w->dZ)
    free(w->dZ);

  if (w->dX_ext)
    free(w->dX_ext);

  if (w->dY_ext)
    free(w->dY_ext);

  if (w->dZ_ext)
    free(w->dZ_ext);

  if (w->cosmphi)
    free(w->cosmphi);

  if (w->sinmphi)
    free(w->sinmphi);

  if (w->Plm)
    free(w->Plm);

  if (w->dPlm)
    free(w->dPlm);

  free(w);
} /* green_free() */

size_t
green_nnm(const green_workspace *w)
{
  return w->nnm;
}

size_t
green_nmax(const green_workspace *w)
{
  return w->nmax;
}

/*
green_calc()
  Compute Green's functions for internal X,Y,Z spherical harmonic
expansion. These are simply the basis functions multiplying the g_{nm}
and h_{nm} coefficients

Inputs: r     - radius (any units)
        theta - colatitude (radians)
        phi   - longitude (radians)
        R     - reference radius (same units as r)
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Plm
w->dPlm
w->sinmphi
w->cosmphi

2) The output Green's functions are stored in w->dX, w->dY, w->dZ

3) The radii r and R can have any units as long as they match
*/

int
green_calc(const double r, const double theta, const double phi,
           const double R, green_workspace *w)
{
  int s = 0;
  size_t n;
  int m;
  const double sint = sin(theta);
  const double cost = cos(theta);
  double ratio = R / r;
  double term = ratio * ratio;     /* (R/r)^{n+2} */

  /* precompute cos(m phi) and sin(m phi) */
  for (n = 0; n <= w->nmax; ++n)
    {
      w->cosmphi[n] = cos(n * phi);
      w->sinmphi[n] = sin(n * phi);
    }

  /* compute associated legendres */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, w->nmax, cost,
                                  w->Plm, w->dPlm);

  for (n = 1; n <= w->nmax; ++n)
    {
      int ni = (int) n;

      /* (a/r)^{n+2} */
      term *= ratio;

      for (m = -ni; m <= ni; ++m)
        {
          int mabs = abs(m);
          size_t cidx = green_nmidx(n, m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);

          if (m < 0)
            {
              /* h_{nm} */
              w->dX[cidx] = term * w->sinmphi[mabs] * w->dPlm[pidx];
              w->dY[cidx] = -term / sint * mabs * w->cosmphi[mabs] * w->Plm[pidx];
              w->dZ[cidx] = -(n + 1.0) * term * w->sinmphi[mabs] * w->Plm[pidx];
            }
          else
            {
              /* g_{nm} */
              w->dX[cidx] = term * w->cosmphi[mabs] * w->dPlm[pidx];
              w->dY[cidx] = term / sint * mabs * w->sinmphi[mabs] * w->Plm[pidx];
              w->dZ[cidx] = -(n + 1.0) * term * w->cosmphi[mabs] * w->Plm[pidx];
            }
        }
    }

  return s;
} /* green_calc() */

/*
green_eval()
  Evaluate spherical harmonic expansion at a given point with
given coefficients

Inputs: r      - radius (any units)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        R      - reference radius (same units as r)
        coeffs - spherical harmonic coefficients
                 must be indexed with green_nmidx()
        B      - (output)
                 B[0] = X component of expansion
                 B[1] = Y component of expansion
                 B[2] = Z component of expansion
        w      - workspace
*/

int
green_eval(const double r, const double theta, const double phi,
           const double R, const double *coeffs,
           double B[3], green_workspace *w)
{
  int s = 0;
  size_t n;

  /* calculate basis functions */
  s = green_calc(r, theta, phi, R, w);

  /* initialize to 0 */
  B[0] = B[1] = B[2] = 0.0;

  for (n = 1; n <= w->nmax; ++n)
    {
      int ni = (int) n;
      int m;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = green_nmidx(n, m);

          B[0] += coeffs[cidx] * w->dX[cidx];
          B[1] += coeffs[cidx] * w->dY[cidx];
          B[2] += coeffs[cidx] * w->dZ[cidx];
        }
    }

  return s;
} /* green_eval() */

/*
green_calc_ext()
  Compute Green's functions for external X,Y,Z spherical harmonic
expansion. These are simply the basis functions multiplying the q_{nm}
and k_{nm} coefficients. Specifically,

V_ext = R sum_{nm} (r/R)^n (q_{nm} cos(m phi) + k_{nm} sin(m phi)) P_{nm}(cos(theta))

and

B_ext = -grad(V_ext)

Inputs: r     - radius (any units)
        theta - colatitude (radians)
        phi   - longitude (radians)
        R     - reference radius (same units as r)
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Plm
w->dPlm
w->sinmphi
w->cosmphi

2) The output Green's functions are stored in w->dX_ext, w->dY_ext,
w->dZ_ext

3) The radii r and R can have any units as long as they match
*/

int
green_calc_ext(const double r, const double theta, const double phi,
               const double R, green_workspace *w)
{
  int s = 0;
  size_t n;
  int m;
  const double sint = sin(theta);
  const double cost = cos(theta);
  double ratio = r / R;
  double term = 1.0;     /* (r/R)^{n-1} */

  /* precompute cos(m phi) and sin(m phi) */
  for (n = 0; n <= w->nmax; ++n)
    {
      w->cosmphi[n] = cos(n * phi);
      w->sinmphi[n] = sin(n * phi);
    }

  /* compute associated legendres */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, w->nmax, cost,
                                  w->Plm, w->dPlm);

  for (n = 1; n <= w->nmax; ++n)
    {
      int ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          int mabs = abs(m);
          size_t cidx = green_nmidx(n, m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);

          if (m < 0)
            {
              /* k_{nm} */
              w->dX_ext[cidx] = term * w->sinmphi[mabs] * w->dPlm[pidx];
              w->dY_ext[cidx] = -term / sint * mabs * w->cosmphi[mabs] * w->Plm[pidx];
              w->dZ_ext[cidx] = n * term * w->sinmphi[mabs] * w->Plm[pidx];
            }
          else
            {
              /* q_{nm} */
              w->dX_ext[cidx] = term * w->cosmphi[mabs] * w->dPlm[pidx];
              w->dY_ext[cidx] = term / sint * mabs * w->sinmphi[mabs] * w->Plm[pidx];
              w->dZ_ext[cidx] = n * term * w->cosmphi[mabs] * w->Plm[pidx];
            }
        }

      /* (r/R)^{n-1} */
      term *= ratio;
    }

  return s;
} /* green_calc_ext() */

/*
green_nmidx()
  This function returns a unique index in [0,w->nnm-1] corresponding
to a given (l,m) pair. The array will look like:

[(1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) ...]

(the (0,0) coefficient is not stored)

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
} /* green_nmidx() */
