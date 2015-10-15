#include <math.h>

#define WGS84_A 6378.1370       /* in km */
#define WGS84_B 6356.752314     /* in km */
#define WGS84_E sqrt(WGS84_A*WGS84_A - WGS84_B*WGS84_B)
#define WGS84_E2 (WGS84_E  * WGS84_E)
#define WGS84_E4 (WGS84_E2 * WGS84_E2)

void geodetic2geocentric(double phi, double h,double *latrad, double *r) /* phi is geodetic latitude in radian*/

{
  double cosphi2, sinphi2, A2, B2, c;

  A2 = WGS84_A*WGS84_A;
  B2 = WGS84_B*WGS84_B;
  cosphi2 = pow(cos(phi),2);
  sinphi2 = pow(sin(phi),2);

  c = h*sqrt(A2*cosphi2+ B2*sinphi2);

  *r = sqrt(h*h + 2*c + (A2*A2*cosphi2 + B2*B2*sinphi2)/(A2*cosphi2 + B2*sinphi2));

  if (fabs(phi-M_PI/2.0) < 0.00001) *latrad = phi;
  else
    {
      *latrad = atan((c + B2)/(c + A2) * tan(phi));
    }
} /* geodetic2geocentric */


void geocentric2geodetic_vec(double theta, double delta, double bx, double by, double bz,
                             double *x, double *y, double *z)

{
  double psi, sp, cp;

  psi = delta-theta; /* compatible with WMM2005 report, sign reversed for colatitudes */
  sp = sin(psi);
  cp = cos(psi);
  *x = bx*cp - bz*sp;
  *y = by;
  *z = bx*sp + bz*cp;
} /* geocentric2geodetic_vec */

