/*
 * ellipsoid.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "ellipsoid.h"

static int ellipsoid_cart2ell(const double r[3], double *mu, double *nu,
                              double *phi, double *a);
static int ellipsoid_coord2basis(const double mu, const double nu, const double phi,
                                 double e_mu[3], double e_nu[3], double e_phi[3]);

/*
ellipsoid_cart2ell()
  Convert Cartesian coordinates (X, Y, Z) to WGS84 ellipsoidal
coordinates (mu,nu,phi) using the relations

X = a cosh(mu) cos(nu) cos(phi)
Y = a cosh(mu) cos(nu) sin(phi)
Z = a sinh(mu) sin(nu)

and the inverse transformation which is detailed on the
Wikipedia page for Oblate Spheroidal coordinates

Inputs: r   - Cartesian position vector in km (X, Y, Z)
        mu  - (output) mu distance parameter (dimensionless)
        nu  - (output) nu eccentric anomaly angle in radians
        phi - (output) phi azimuth angle in radians \in [-pi,pi]
        a   - (output) constant defining reference ellipsoid in km

Return: success or error
*/

static int
ellipsoid_cart2ell(const double r[3], double *mu, double *nu,
                   double *phi, double *a)
{
  int s = 0;
  double X = r[0];
  double Y = r[1];
  double Z = r[2];
  double ratio = WGS84_B / WGS84_A;
  double rho;
  double complex arg;

  /* calculate constant a for WGS84 system */
  *a = WGS84_A * sqrt(1.0 - ratio * ratio);

  /* calculate aximuth angle */
  *phi = atan2(Y, X);

  /*
   * calculate mu and eccentric anomaly angle nu; these
   * formulas were taken from the Mathematica package
   */

  /* compute cylindrical radius sqrt(X^2 + Y^2) */
  rho = sqrt(X * X + Y * Y);

  arg = cacosh((rho + Z * I) / *a);

  *mu = creal(arg);
  *nu = cimag(arg);

  return s;
} /* ellipsoid_cart2ell() */

/*
ellipsoid_coord2basis()
  Compute oblate spheroidal basis vectors e_mu,e_nu,e_phi,
given coordinates (mu,nu,phi)

Inputs: mu   - distance parameter (dimensionless)
        nu   - nu eccentric anomaly angle in radians
        phi  - longitude (radians)
        e_mu - (output)
*/

static int
ellipsoid_coord2basis(const double mu, const double nu, const double phi,
                      double e_mu[3], double e_nu[3], double e_phi[3])
{
  int s = 0;
  double C; /* normalization constant */

  C = sqrt(sinh(mu) * sinh(mu) + sin(nu) * sin(nu));

  e_mu[0] = sinh(mu) * cos(nu) * cos(phi) / C;
  e_mu[1] = sinh(mu) * cos(nu) * sin(phi) / C;
  e_mu[2] = cosh(mu) * sin(nu) / C;

  e_nu[0] = -cosh(mu) * sin(nu) * cos(phi) / C;
  e_nu[1] = -cosh(mu) * sin(nu) * sin(phi) / C;
  e_nu[2] = sinh(mu) * cos(nu) / C;

  e_phi[0] = -sin(phi);
  e_phi[1] = cos(phi);
  e_phi[2] = 0.0;

  return s;
} /* ellipsoid_coord2basis() */

/*
ellipsoid_basis()
  Compute basis vectors in ellipsoidal coordinates defined by
WGS84 (actually oblate spheroidal coordinates) at a given point.
Using the relations

X = a cosh(mu) cos(nu) cos(phi)
Y = a cosh(mu) cos(nu) sin(phi)
Z = a sinh(mu) sin(nu)

r = X i + Y j + Z k

e_mu = dr/dmu / |dr/dmu|
e_nu = dr/dnu / |dr/dnu|
e_phi = dr/dphi / |dr/dphi|

Inputs: r     - ECEF Cartesian position vector (X,Y,Z) in km
        e_mu  - (output) Cartesian components of e_mu
        e_nu  - (output) Cartesian components of e_nu
        e_phi - (output) Cartesian components of e_phi

Return: success or error
*/

int
ellipsoid_basis(const double r[3], double e_mu[3], double e_nu[3],
                double e_phi[3])
{
  int s = 0;
  double mu, nu, phi, a;

  /* calculate ellipsoid parameters from Cartesian */
  s += ellipsoid_cart2ell(r, &mu, &nu, &phi, &a);

  /* compute basis vectors */
  s += ellipsoid_coord2basis(mu, nu, phi, e_mu, e_nu, e_phi);

  return s;
} /* ellipsoid_basis() */

/*
ellipsoid_basis_mu()
  Compute basis vectors in ellipsoidal coordinates defined by
WGS84 (actually oblate spheroidal coordinates) at a given point,
except with a fixed distance parameter mu, using the relations:

X = a cosh(mu) cos(nu) cos(phi)
Y = a cosh(mu) cos(nu) sin(phi)
Z = a sinh(mu) sin(nu)

r = X i + Y j + Z k

e_mu = dr/dmu / |dr/dmu|
e_nu = dr/dnu / |dr/dnu|
e_phi = dr/dphi / |dr/dphi|

Specifying a value of 'mu' is useful, for example, to compute
basis vectors on the Earth's surface ellipsoid rather than
the local ellipsoid to the point 'r'

Inputs: r     - Cartesian position vector (X,Y,Z) in km
        mu    - distance parameter specifying ellipsoid on which
                to compute basis vectors
        e_mu  - (output) Cartesian components of e_mu
        e_nu  - (output) Cartesian components of e_nu
        e_phi - (output) Cartesian components of e_phi

Return: success or error
*/

int
ellipsoid_basis_mu(const double r[3], const double mu, double e_mu[3],
                   double e_nu[3], double e_phi[3])
{
  int s = 0;
  double dummy, nu, phi, a;

  /* calculate ellipsoid parameters from Cartesian */
  s += ellipsoid_cart2ell(r, &dummy, &nu, &phi, &a);

  /* compute basis vectors */
  s += ellipsoid_coord2basis(mu, nu, phi, e_mu, e_nu, e_phi);

  return s;
} /* ellipsoid_basis_mu() */
