/*
 * geo.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "geo.h"

#define WGS84_A 6378.1370       /* in km */
#define WGS84_B 6356.752314     /* in km */

/*
geo2geodetic()
  Convert geographic/geocentric (lat,lon,r) to geodetic (lat,lon,alt)
Longitude is the same in both coordinate systems

Inputs: lat    - geocentric latitude (radians)
        lon    - longitude (radians)
        radius - geocentric radius (km)
        latd   - (output) geodetic latitude (radians)
        altd   - (output) geodetic altitude (km)

Return: success or error
*/

int
geo2geodetic(double lat, double lon, double radius, double *latd, double *altd)
{
  int status = 0;
  double Re = 6378.137; /* equatorial radius km */
  double Rp = 6356.752; /* polar radius km */
  double e;             /* ellipse eccentricity */
  double x, y, z, r;
  double s, t0, dzeta1, xi1, rho1, c1, s1, b1;
  double u1, w1;

  /* calculate planet eccentricity */
  e = sqrt(Re*Re - Rp*Rp) / Re;

  x = radius * cos(lat) * cos(lon);
  y = radius * cos(lat) * sin(lon);
  z = radius * sin(lat);
  r = sqrt(x*x + y*y);

  s = sqrt(r*r + z*z) * (1.0 - Re * sqrt((1.0 - e*e) / ((1.0 - e*e) * r * r + z * z)));
  t0 = 1.0 + s * sqrt(1.0 - (e*z)*(e*z) / (r*r + z*z)) / Re;
  dzeta1 = z * t0;
  xi1 = r * (t0 - e*e);
  rho1 = sqrt(xi1*xi1 + dzeta1*dzeta1);
  c1 = xi1 / rho1;
  s1 = dzeta1 / rho1;
  b1 = Re / sqrt(1.0 - e*e*s1*s1);
  u1 = b1 * c1;
  w1 = b1 * s1 * (1.0 - e*e);

  *altd = sqrt(pow(r - u1, 2.0) + pow(z - w1, 2.0));
  *latd = atan2(s1, c1);

  return status;
} /* geo2geodetic() */

#if 0
/*
geodetic2geo()
  Convert geodetic height/latitude to geocentric height/latitude

Inputs: phi    - geodetic latitude (radians)
        h      - geodetic altitude (km)
        latrad - (output) geocentric latitude (radians)
        r      - (output) geocentric radius (km)
*/

void
geodetic2geo(double phi, double h, double *latrad, double *r)
{
  double cosphi2, sinphi2, A2, B2, c;

  A2 = WGS84_A * WGS84_A;
  B2 = WGS84_B * WGS84_B;
  cosphi2 = pow(cos(phi), 2.0);
  sinphi2 = pow(sin(phi), 2.0);

  c = h * sqrt(A2 * cosphi2 + B2 * sinphi2);

  *r = sqrt(h*h + 2*c + (A2*A2*cosphi2 + B2*B2*sinphi2)/(A2*cosphi2 + B2*sinphi2));

  if (fabs(phi-M_PI/2.0) < 0.00001)
    *latrad = phi;
  else
    {
      *latrad = atan((c + B2)/(c + A2) * tan(phi));
    }
} /* geodetic2geo() */
#endif /* 0 */

/*
geodetic2geo()
  Convert geodetic to spherical geocentric using WMM report equations
(7)-(8)

Inputs: latd   - geodetic latitude (radians)
        h      - geodetic height above ellipsoid (km)
        latrad - (output) geocentric latitude
        r      - (output) geocentric radius (km)
*/

void
geodetic2geo(double latd, double h, double *latrad, double *r)
{
  const double a = WGS84_A;
  const double b = WGS84_B;
  const double epssq = 0.0066943800042605908;
  double CosLat, SinLat, rc, xp, zp; /*all local variables */

  /*
   ** Convert geodetic coordinates, (defined by the WGS-84
   ** reference ellipsoid), to Earth Centered Earth Fixed Cartesian
   ** coordinates, and then to spherical coordinates.
   */

  CosLat = cos(latd);
  SinLat = sin(latd);

  /* compute the local radius of curvature on the WGS-84 reference ellipsoid */

  rc = a / sqrt(1.0 - epssq * SinLat * SinLat);

  /* compute ECEF Cartesian coordinates of specified point (for longitude=0) */

  xp = (rc + h) * CosLat;
  zp = (rc * (1.0 - epssq) + h) * SinLat;

  /* compute spherical radius and angle lambda and latd of specified point */
  *r = sqrt(xp * xp + zp * zp);
  *latrad = asin(zp / *r); /* geocentric latitude */
}
