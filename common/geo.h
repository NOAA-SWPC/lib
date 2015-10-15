/*
 * geo.h
 * Patrick Alken
 */

#ifndef INCLUDED_geo_h
#define INCLUDED_geo_h

/*
 * Prototypes
 */

int geo2geodetic(double latc, double lon, double radius, double *latd, double *altd);
void geodetic2geo(double phi, double h, double *latrad, double *r);

#endif /* INCLUDED_geo_h */
