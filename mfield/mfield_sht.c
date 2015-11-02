/*
 * mfield_sht.c
 *
 * This module is contains wrapper functions for the NFFT library,
 * to perform spherical harmonic expansions on scattered data points
 * quickly
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <complex.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>

#include <nfft3.h>

#include "mfield_eval.h"
#include "mfield_sht.h"

static int transform_Y(const gsl_vector *g, const double r,
                       const gsl_vector *theta, gsl_vector *Y,
                       mfield_sht_workspace *w);
static int transform_Z(const gsl_vector *g, const double r,
                       gsl_vector *Z, mfield_sht_workspace *w);

/*
mfield_sht_alloc()
  Allocate a mfield_sht workspace

Inputs: N - maximum spherical harmonic degree
        M - number of data points (r,theta,phi) to process
            at a time
*/

mfield_sht_workspace *
mfield_sht_alloc(const size_t N, const size_t M)
{
  mfield_sht_workspace *w;

  w = calloc(1, sizeof(mfield_sht_workspace));
  if (!w)
    {
      fprintf(stderr, "mfield_sht_alloc: calloc failed: %s\n",
              strerror(errno));
      return 0;
    }

  /* precomputation step - inits global data */
  nfsft_precompute(N, 1000, 0U, 0U);

  nfsft_init(&(w->plan), N, M);

  w->N = N;
  w->M = M;

  return w;
}

void
mfield_sht_free(mfield_sht_workspace *w)
{
  nfsft_finalize(&(w->plan));
  nfsft_forget();

  free(w);
}

/*
mfield_sht_calc()

  Calculate spherical harmonic expansions on scattered data

(r,theta_j,phi_j), j = 0, ..., M - 1

Currently, the radius is fixed due to lack of support in NFFT
for arbitrary radii

The magnetic field components X,Y,Z evaluated at the M points
are returned

Inputs: g     - model coefficients (semi-Schmidt normalized, nT)
        r     - radial point (dimensionless: r~ = r / a)
        theta - colatitude points (radians)
        phi   - longitude points (radians) in [0,2pi]
        Y     - (output) Y field values at node points (nT)
        Z     - (output) Z field values at node points (nT)
        w     - workspace

Notes:
1) Vector lengths of theta,phi must match previously specified w->m
*/

int
mfield_sht_calc(const gsl_vector *g, const double r,
                const gsl_vector *theta,
                const gsl_vector *phi, gsl_vector *Y,
                gsl_vector *Z,
                mfield_sht_workspace *w)
{
  const size_t M = w->M;

  if (theta->size != M)
    {
      GSL_ERROR("theta vector has wrong length", GSL_EBADLEN);
    }
  else if (phi->size != M)
    {
      GSL_ERROR("phi vector has wrong length", GSL_EBADLEN);
    }
  else if (Y->size != M)
    {
      GSL_ERROR("Y vector has wrong length", GSL_EBADLEN);
    }
  else if (Z->size != M)
    {
      GSL_ERROR("Z vector has wrong length", GSL_EBADLEN);
    }
  else
    {
      int s = GSL_SUCCESS;
      size_t i;

      /* scale phi to [-1/2,1/2] and theta to [0,1/2] for input to nfsft */
      for (i = 0; i < M; ++i)
        {
          double ti = gsl_vector_get(theta, i) / (2.0 * M_PI);
          double pi = gsl_vector_get(phi, i) / (2.0 * M_PI);

          if (pi >= 0.5)
            pi -= 1.0;

          w->plan.x[2*i] = pi;
          w->plan.x[2*i + 1] = ti;
        }

      /* node-dependent precomputation */
      nfsft_precompute_x(&(w->plan));

      /* perform vector transforms */
      transform_Y(g, r, theta, Y, w);
      transform_Z(g, r, Z, w);

      return s;
    }
}

/*
transform_Y()
  Perform spherical harmonic transform for Y component

Inputs: g     - model coefficients (nT)
        r     - radial point
        theta - colatitude node points (radians)
        Y     - (output) where to store Y vector (nT)
        w     - workspace

Notes:
1) The (theta,phi) values must be input to w->plan prior to
calling this function
*/

static int
transform_Y(const gsl_vector *g, const double r, const gsl_vector *theta,
            gsl_vector *Y, mfield_sht_workspace *w)
{
  const size_t N = w->N;
  const size_t M = w->M;
  size_t n, i;
  int m;
  double rterm = 1.0 / (r * r);

  /* fill in Y coefficients in w->plan */
  for (n = 1; n <= N; ++n)
    {
      int ni = (int) n;

      rterm /= r;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_eval_nmidx(n, m);
          double gnm = gsl_vector_get(g, cidx);

          /* radial term (a/r)^(n+2) */
          gnm *= rterm;

          /* Y prefactor */
          gnm *= (double) m;

          /*
           * account for differences in normalization of ALFs; nfft
           * uses ALFs which differ by sqrt(2) from semi-Schmidt
           * for m != 0
           */
          if (m != 0)
            gnm *= sqrt(2.0);

          if (m >= 0)
            w->plan.f_hat[NFSFT_INDEX(n, m, &(w->plan))] = gnm + 0.0 * I;
          else
            w->plan.f_hat[NFSFT_INDEX(n, m, &(w->plan))] = gnm * I;
        }
    }

  /* set the (0,0) coefficient to 0 */
  w->plan.f_hat[NFSFT_INDEX(0, 0, &(w->plan))] = 0.0;

  /* perform the Y transform */
  nfsft_trafo(&(w->plan));

  /* store solution in output vectors */
  for (i = 0; i < M; ++i)
    {
      double complex f = w->plan.f[i];
      double ti = gsl_vector_get(theta, i);

      gsl_vector_set(Y, i, cimag(f) / sin(ti));
    }

  return GSL_SUCCESS;
} /* transform_Y() */

/*
transform_Z()
  Perform spherical harmonic transform for Z component

Inputs: g - model coefficients (nT)
        Z - (output) where to store Z vector (nT)
        w - workspace

Notes:
1) The (theta,phi) values must be input to w->plan prior to
calling this function
*/

static int
transform_Z(const gsl_vector *g, const double r,
            gsl_vector *Z, mfield_sht_workspace *w)
{
  const size_t N = w->N;
  const size_t M = w->M;
  size_t n, i;
  int m;
  double rterm = 1.0 / (r * r);

  /* fill in Z coefficients in w->plan */
  for (n = 1; n <= N; ++n)
    {
      int ni = (int) n;

      rterm /= r;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_eval_nmidx(n, m);
          double gnm = gsl_vector_get(g, cidx);

          /* radial term (a/r)^(n+2) */
          gnm *= rterm;

          /* Z prefactor */
          gnm *= -(n + 1.0);

          /*
           * account for differences in normalization of ALFs; nfft
           * uses ALFs which differ by sqrt(2) from semi-Schmidt
           * for m != 0
           */
          if (m != 0)
            gnm *= sqrt(2.0);

          if (m >= 0)
            w->plan.f_hat[NFSFT_INDEX(n, m, &(w->plan))] = gnm + 0.0 * I;
          else
            w->plan.f_hat[NFSFT_INDEX(n, m, &(w->plan))] = gnm * I;
        }
    }

  /* set the (0,0) coefficient to 0 */
  w->plan.f_hat[NFSFT_INDEX(0, 0, &(w->plan))] = 0.0;

  /* perform the Z transform */
  nfsft_trafo(&(w->plan));

  /* store solution in output vectors */
  for (i = 0; i < M; ++i)
    {
      double complex f = w->plan.f[i];
      gsl_vector_set(Z, i, creal(f));
    }

  return GSL_SUCCESS;
} /* transform_Z() */
