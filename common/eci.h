/*
 * eci.h
 */

#ifndef INCLUDED_eci_h
#define INCLUDED_eci_h

#include <time.h>

/*
 * Prototypes
 */

int eci2ecef_pos(const time_t t, const double r_ECI[3], double r_ECEF[3]);
int ecef2eci_pos(const time_t t, const double r_ECI[3], double r_ECEF[3]);
int eci2sph_pos(const time_t t, const double r_ECI[3], double r_sph[3]);
int sph2eci_pos(const time_t t, const double r, const double theta,
                const double phi, double r_ECI[3]);

#endif /* INCLUDED_eci_h */
