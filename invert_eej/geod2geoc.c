#include <math.h>
#include <gsl/gsl_math.h>

#include "geod2geoc.h"

#define WGS84_A 6378.1370       /* in km */
#define WGS84_B 6356.752314     /* in km */

/*
geod2geoc()
  Convert geodetic height/latitude to geocentric height/latitude

Inputs: phi    - geodetic latitude (radians)
        h      - geodetic altitude (km)
        latrad - (output) geocentric latitude (radians)
        r      - (output) geocentric altitude (km)
*/

void
geod2geoc(double phi, double h, double *latrad, double *r)
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
} /* geod2geoc() */
