#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>

#include "mfield_green.h"

mfield_green_workspace *
mfield_green_alloc(const size_t nmax, const double R)
{
  mfield_green_workspace *w;
  const size_t plm_size = gsl_sf_legendre_array_n(nmax);
  const size_t nnm_tot = (nmax + 1) * (nmax + 1);

  w = calloc(1, sizeof(mfield_green_workspace));
  if (!w)
    return 0;

  w->nmax = nmax;
  w->nnm = nnm_tot - 1; /* exclude (0,0) coefficient */
  w->R = R;

  w->cosmphi = malloc((nmax + 1) * sizeof(double));
  w->sinmphi = malloc((nmax + 1) * sizeof(double));
  if (!w->cosmphi || !w->sinmphi)
    {
      mfield_green_free(w);
      fprintf(stderr, "mfield_green_alloc: cannot allocate space for cos/sin arrays: %s\n",
              strerror(errno));
      return 0;
    }

  w->Plm = malloc(plm_size * sizeof(double));
  w->dPlm = malloc(plm_size * sizeof(double));
  if (!w->Plm || !w->dPlm)
    {
      mfield_green_free(w);
      fprintf(stderr, "mfield_green_alloc: cannot allocate space for legendre arrays: %s\n",
              strerror(errno));
      return 0;
    }

  w->dX = malloc(w->nnm * sizeof(double));
  w->dY = malloc(w->nnm * sizeof(double));
  w->dZ = malloc(w->nnm * sizeof(double));
  w->dX_ext = malloc(w->nnm * sizeof(double));
  w->dY_ext = malloc(w->nnm * sizeof(double));
  w->dZ_ext = malloc(w->nnm * sizeof(double));
  if (!w->dX || !w->dY || !w->dZ || !w->dX_ext || !w->dY_ext || !w->dZ_ext)
    {
      mfield_green_free(w);
      fprintf(stderr, "mfield_green_alloc: cannot allocate space for result arrays: %s\n",
              strerror(errno));
      return 0;
    }

  return w;
} /* mfield_green_alloc() */

void
mfield_green_free(mfield_green_workspace *w)
{
  if (w->cosmphi)
    free(w->cosmphi);

  if (w->sinmphi)
    free(w->sinmphi);

  if (w->Plm)
    free(w->Plm);

  if (w->dPlm)
    free(w->dPlm);

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

  free(w);
}

/*
mfield_green_calc()
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
*/

int
mfield_green_calc(const double r, const double theta, const double phi,
                  mfield_green_workspace *w)
{
  int s = 0;
  size_t n;
  int m;
  const double sint = sin(theta);
  const double cost = cos(theta);
  double ratio = w->R / r;
  double term = ratio * ratio;     /* (a/r)^{n+2} */

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
          size_t cidx = mfield_green_nmidx(n, m);
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
} /* mfield_green() */

/*
mfield_green_ext()
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
mfield_green_ext(const double r, const double theta, const double phi,
                 mfield_green_workspace *w)
{
  int s = 0;
  size_t n;
  int m;
  const double sint = sin(theta);
  const double cost = cos(theta);
  double ratio = r / w->R;
  double term = 1.0;     /* (r/a)^{n-1} */

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
          size_t cidx = mfield_green_nmidx(n, m);
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
} /* mfield_green_ext() */

/*
mfield_green_nmidx()
  This function returns a unique index in [0,w->p-1] corresponding
to a given (l,m) pair. The array will look like:

[(1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) ...]

(the (0,0) coefficient is not solved for)

Inputs: l - SH degree (> 0)
        m - SH order (-l <= m <= l)

Return: index in [0,nnm-1]
*/

size_t
mfield_green_nmidx(const size_t n, const int m)
{
  size_t base = n * n; /* index of block for this n */
  int offset = m + n;  /* offset within block for this m */
  size_t nmidx;

  if (n == 0)
    {
      fprintf(stderr, "mfield_green_nmidx: error: n = 0\n");
      return 0;
    }

  nmidx = base + offset;

  /* subtract 1 to exclude (0,0) coefficient */
  return nmidx - 1;
} /* mfield_green_nmidx() */

double
mfield_green_X(const size_t idx, const mfield_green_workspace *w)
{
  return w->dX[idx];
}

double
mfield_green_Y(const size_t idx, const mfield_green_workspace *w)
{
  return w->dY[idx];
}

double
mfield_green_Z(const size_t idx, const mfield_green_workspace *w)
{
  return w->dZ[idx];
}
