/*
 * msynth_green.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>

#include "msynth.h"

/*
msynth_green()
  Compute Green's functions for X,Y,Z spherical harmonic expansion. These
are simply the basis functions multiplying the g_{nm} and h_{nm} coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Plm
w->dPlm
w->sinmphi
w->cosmphi

2) The output Green's functions are stored in w->dX, w->dY, w->dZ

3) Green's functions are computed only up to degree w->eval_nmax
*/

int
msynth_green(const double r, const double theta, const double phi,
             msynth_workspace *w)
{
  int s = 0;
  const size_t nmax = w->eval_nmax;
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
          size_t cidx = msynth_nmidx(n, m, w);
          size_t aidx = gsl_sf_legendre_array_index(n, mabs);

          if (m < 0)
            {
              /* h_{nm} */
              w->dX[cidx] = term * w->sinmphi[mabs] * w->dPlm[aidx] * (-sint);
              w->dY[cidx] = -term / sint * mabs * w->cosmphi[mabs] * w->Plm[aidx];
              w->dZ[cidx] = -(n + 1.0) * term * w->sinmphi[mabs] * w->Plm[aidx];
            }
          else
            {
              /* g_{nm} */
              w->dX[cidx] = term * w->cosmphi[mabs] * w->dPlm[aidx] * (-sint);
              w->dY[cidx] = term / sint * mabs * w->sinmphi[mabs] * w->Plm[aidx];
              w->dZ[cidx] = -(n + 1.0) * term * w->cosmphi[mabs] * w->Plm[aidx];
            }
        }
    }

  return s;
} /* msynth_green() */

/*
msynth_green_ext()
  Compute Green's functions for X,Y,Z spherical harmonic expansion. These
are simply the basis functions multiplying the g_{nm} and h_{nm} coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Plm
w->dPlm
w->sinmphi
w->cosmphi

2) The output Green's functions are stored in w->dX_ext, w->dY_ext, w->dZ_ext
*/

int
msynth_green_ext(const double r, const double theta, const double phi,
                 msynth_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax_ext;
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
          size_t cidx = msynth_nmidx(n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);

          if (m < 0)
            {
              /* h_{nm} */
              w->dX_ext[cidx] = term * w->sinmphi[mabs] * w->dPlm[pidx];
              w->dY_ext[cidx] = -term / sint * mabs * w->cosmphi[mabs] * w->Plm[pidx];
              w->dZ_ext[cidx] = (double) n * term * w->sinmphi[mabs] * w->Plm[pidx];
            }
          else
            {
              /* g_{nm} */
              w->dX_ext[cidx] = term * w->cosmphi[mabs] * w->dPlm[pidx];
              w->dY_ext[cidx] = term / sint * mabs * w->sinmphi[mabs] * w->Plm[pidx];
              w->dZ_ext[cidx] = (double) n * term * w->cosmphi[mabs] * w->Plm[pidx];
            }
        }

      /* (r/a)^{n-1} */
      term *= ratio;
    }

  return s;
} /* msynth_green_ext() */
