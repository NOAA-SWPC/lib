/*
 * eci.c
 *
 * Routines related to Earth-Centered Inertial (ECI)
 * coordinates
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "common.h"
#include "eci.h"
#include "julian.h"
#include "ut1.h"

ut1_workspace *eci_ut1_p = NULL;

/*
timet2GRA
  Compute Greenwich right ascension at given time

Inputs: t - timestamp

Return: right ascension angle in radians
*/

static double
timet2GRA(const time_t t)
{
  /* sidereal time, 1 Jan 1970 FK5 */
  const double THGR70 = 1.7321343856509374;

  /* sidereal rate FK5 */
  const double C1 = 1.72027916940703639e-2;
  const double C1P2P = C1 + 2.0 * M_PI;

  /* sidereal acceleration FK5 */
  const double C1DOT = 5.07551419432269442e-15;

  double jd, d70, id70, frac, theta;
  double dut1, t1;

  if (eci_ut1_p == NULL)
    eci_ut1_p = ut1_alloc(UT1_DATA_FILE);

  ut1_get(t, &dut1, eci_ut1_p);

  /* correct for dUT1 */
  t1 = t + dut1;

  /* compute julian day */
  jd = (t1 / 86400.0) + 2440587.5;

  /* number of days since Jan 1 1970 */
  d70 = jd - 2440586.5;

  /* integer part of d70 */
  id70 = (double) ((int) d70);

  /* fractional part of d70 */
  frac = d70 - id70;

  theta = THGR70 + C1 * id70 + C1P2P * frac +
          (d70*d70) * C1DOT;

  return wrap2pi(theta);
}

#if 0
/*
timet2GRA
  Compute Greenwich right ascension at given time

Inputs: t - timestamp

Return: right ascension angle in radians
*/

static double
timet2GRA(const time_t t)
{
  /* sidereal time, 1 Jan 1970 FK5 */
  const double THGR70 = 1.7321343856509374;

  /* sidereal rate FK5 */
  const double C1 = 1.72027916940703639e-2;
  const double C1P2P = C1 + 2.0 * M_PI;

  /* sidereal acceleration FK5 */
  const double C1DOT = 5.07551419432269442e-15;

  double jd = timet2julian(t); /* julian day */
  double d70 = jd - 2440586.5; /* days since Jan 1 1970 */
  double id70 = (double) ((int) d70); /* discard fractional part */
  double frac = d70 - id70;
  double theta;

  theta = THGR70 + C1 * id70 + C1P2P * frac +
          (d70*d70) * C1DOT;

  return wrap2pi(theta);
}
#endif /* 0 */

/*
eci2ecef_pos()
  Transform position vector in ECI coordinates to ECEF

Inputs: t      - timestamp
        r_ECI  - ECI (X,Y,Z) position vector
        r_ECEF - (output) ECEF (X,Y,Z) position vector
*/

int
eci2ecef_pos(const time_t t, const double r_ECI[3], double r_ECEF[3])
{
  int s = 0;
  double jd, GAST;

  /* compute julian day from unix time */
  jd = (t / 86400.0) + 2440587.5;

  /* compute Greenwich apparent sidereal time in radians */
  GAST = julian2GAST(jd);

  /* rotate ECI vector to ECEF */
  r_ECEF[0] = r_ECI[0] * cos(GAST) + r_ECI[1] * sin(GAST);
  r_ECEF[1] = -r_ECI[0] * sin(GAST) + r_ECI[1] * cos(GAST);
  r_ECEF[2] = r_ECI[2];

  return s;
} /* eci2ecef_pos() */

/*
eci2ecef()
  Transform position and velocity vectors in ECI coordinates to ECEF

Inputs: t      - timestamp
        r_ECI  - ECI (X,Y,Z) position vector
        v_ECI  - ECI (X,Y,Z) velocity vector
        r_ECEF - (output) ECEF (X,Y,Z) position vector
        v_ECEF - (output) ECEF (X,Y,Z) velocity vector
*/

int
eci2ecef(const time_t t, const double r_ECI[3], const double v_ECI[3],
         double r_ECEF[3], double v_ECEF[3])
{
  int s = 0;
  const double omega_e = 7.29211585275553e-5; /* average intertial rotation rate of Earth (rad/s) */
  double jd, GAST;
  double cs, sn;

  /* compute julian day from unix time */
  jd = (t / 86400.0) + 2440587.5;

  /* compute Greenwich apparent sidereal time in radians */
  GAST = julian2GAST(jd);

  cs = cos(GAST);
  sn = sin(GAST);

  /* rotate ECI position vector to ECEF */
  r_ECEF[0] = r_ECI[0] * cs + r_ECI[1] * sn;
  r_ECEF[1] = -r_ECI[0] * sn + r_ECI[1] * cs;
  r_ECEF[2] = r_ECI[2];

  /* compute ECEF velocity vector */

  v_ECEF[0] = v_ECI[0] * cs + v_ECI[1] * sn;
  v_ECEF[1] = -v_ECI[0] * sn + v_ECI[1] * cs;
  v_ECEF[2] = v_ECI[2];

  v_ECEF[0] += -omega_e * sn * r_ECI[0] + omega_e * cs * r_ECI[1];
  v_ECEF[1] += -omega_e * cs * r_ECI[0] - omega_e * sn * r_ECI[1];

  return s;
}

/*
ecef2eci_pos()
  Transform position vector in ECEF coordinates to ECI

Inputs: t      - timestamp
        r_ECEF - ECEF (X,Y,Z) position vector
        r_ECI  - (output) ECI (X,Y,Z) position vector
*/

int
ecef2eci_pos(const time_t t, const double r_ECEF[3], double r_ECI[3])
{
  int s = 0;
  double jd, GAST;

  /* compute julian day from unix time */
  jd = (t / 86400.0) + 2440587.5;

  /* compute Greenwich apparent sidereal time in radians */
  GAST = julian2GAST(jd);

  /* rotate ECI vector to ECEF */
  r_ECI[0] = r_ECEF[0] * cos(GAST) - r_ECEF[1] * sin(GAST);
  r_ECI[1] = r_ECEF[0] * sin(GAST) + r_ECEF[1] * cos(GAST);
  r_ECI[2] = r_ECEF[2];

  return s;
} /* ecef2eci_pos() */

/*
eci2sph_pos()
  Transform position vector in ECI coordinates (X,Y,Z) to
geocentric spherical (r,theta,phi)

r     = || r_ECI ||
theta = acos(Z / r)
phi   = ra(X,Y,Z) - ra(Greenwich)

where ra() is the right ascension angle in the inertial
celestial sphere

Inputs: t     - timestamp
        r_ECI - ECI position vector (X,Y,Z)
        r_sph - (output) spherical coordinate vector
                r_sph[0] = radius (same units as r_ECI)
                r_sph[1] = geocentric colatitude theta (rad)
                r_sph[2] = geocentric longitude phi (rad)

Notes:

1) 20 May 2013: this function (and in particular time2GRA()) has
been validated against the SGP4 equivalent functions
*/

int
eci2sph_pos(const time_t t, const double r_ECI[3], double r_sph[3])
{
  int s = 0;
  gsl_vector_const_view v = gsl_vector_const_view_array(r_ECI, 3);
  double r = gsl_blas_dnrm2(&v.vector);     /* r is same in ECI and ECEF */
  double r_asc = atan2(r_ECI[1], r_ECI[0]); /* right ascension of position vector */
  double gr_asc = timet2GRA(t);             /* right ascension of Greenwich */

  r_sph[0] = r;
  r_sph[1] = acos(r_ECI[2] / r);
  r_sph[2] = r_asc - gr_asc;

  return s;
} /* eci2sph_pos() */

/*
sph2eci_pos()
  Transform position vector from spherical coordinates to ECI

Inputs: t     - timestamp
        r     - geocentric radius
        theta - geocentric colatitude (radians)
        phi   - longitude (radians)
        r_ECI - (output) ECI position vector (X,Y,Z)
                same units as r
*/

int
sph2eci_pos(const time_t t, const double r, const double theta,
            const double phi, double r_ECI[3])
{
  int s;
  double r_ECEF[3];

  r_ECEF[0] = r * sin(theta) * cos(phi);
  r_ECEF[1] = r * sin(theta) * sin(phi);
  r_ECEF[2] = r * cos(theta);

  s = ecef2eci_pos(t, r_ECEF, r_ECI);

  return s;
} /* sph2eci_pos() */
