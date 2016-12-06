/*
 * current.c
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
green_eval_sheet_int()
  Compute equivalent sheet current density at a given radius using
internal Gauss coefficients gnm

Inputs: b     - radius of spherical shell (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        gnm   - gnm coefficients
        K     - (output) sheet current density in A/km in (X,Y,Z)
        w     - workspace

Return: success/error
*/

int
green_eval_sheet_int(const double b, const double theta, const double phi,
                     const double *gnm, double K[3], green_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax;
  const size_t mmax = w->mmax;
  const double ratio = w->R / b;
  const double sint = sin(theta);
  size_t n;
  int m;
  double rterm = pow(ratio, 2.0);

  K[0] = 0.0;
  K[1] = 0.0;
  K[2] = 0.0;

  /* precompute cos(m phi) and sin(m phi) */
  for (m = 0; m <= (int) mmax; ++m)
    {
      w->cosmphi[m] = cos(m * phi);
      w->sinmphi[m] = sin(m * phi);
    }

  /* compute legendre functions */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                  cos(theta), w->Plm, w->dPlm);

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;
      double nfac = -(2.0 * n + 1.0) / (double) n;

      /* (R/b)^(n+3) */
      rterm *= ratio;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = green_nmidx(n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          double qnm;
          double Ynm, dYnm;

          qnm = nfac * rterm * gnm[nmidx];

          if (m >= 0)
            {
              Ynm = w->Plm[pidx] * w->cosmphi[mabs];
              dYnm = w->dPlm[pidx] * w->cosmphi[mabs];
            }
          else
            {
              Ynm = w->Plm[pidx] * w->sinmphi[mabs];
              dYnm = w->dPlm[pidx] * w->sinmphi[mabs];
            }

          K[0] += (double)m / sint * qnm * Ynm;
          K[1] += qnm * dYnm;
        }
    }

  /* convert to A/km */
  K[0] *= 1.0e3 / (GREEN_MU_0 * ratio);
  K[1] *= 1.0e3 / (GREEN_MU_0 * ratio);

  return s;
}

/*
green_eval_chi_int()
  Compute equivalent current stream function at a given point using
internal Gauss coefficients gnm

Inputs: b     - radius of spherical current shell (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        gnm   - gnm coefficients
        w     - workspace

Return: current stream function in units of kA
*/

double
green_eval_chi_int(const double b, const double theta, const double phi,
                   const double *gnm, green_workspace *w)
{
  const size_t nmax = w->nmax;
  const size_t mmax = w->mmax;
  const double ratio = b / w->R;
  const double sint = sin(theta);
  size_t n;
  int m;
  double rterm = pow(ratio, 2.0);
  double chi = 0.0;

  /* precompute cos(m phi) and sin(m phi) */
  for (m = 0; m <= (int) mmax; ++m)
    {
      w->cosmphi[m] = cos(m * phi);
      w->sinmphi[m] = sin(m * phi);
    }

  /* compute legendre functions */
  gsl_sf_legendre_array(GSL_SF_LEGENDRE_SCHMIDT, nmax, cos(theta), w->Plm);

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;
      double nfac = -(2.0 * n + 1.0) / (double) n;

      /* (R/b)^(n+3) */
      rterm *= ratio;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = green_nmidx(n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          double qnm;
          double Ynm;

          qnm = nfac * rterm * gnm[nmidx];

          if (m >= 0)
            Ynm = w->Plm[pidx] * w->cosmphi[mabs];
          else
            Ynm = w->Plm[pidx] * w->sinmphi[mabs];

          chi += qnm * Ynm;
        }
    }

  /* convert to A/km */
  chi *= -ratio * (b / GREEN_MU_0);

  return chi;
}
