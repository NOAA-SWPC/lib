#include <stdio.h>

#include <gsl/gsl_math.h>

#include "coord.h"
#include "Gdef.h"

/*
get_zenith()
  Calculate the solar zenith angle for given fday, latitude
and longitude

Inputs: fday      - floating point days since 1-Jan-2000 00:00:00 UTC
        latitude  - latitude (radians)
        longitude - longitude (radians)

Return: zenith angle in radians
*/

double
get_zenith(double fday, double latitude, double longitude)
{
  double gse_lon, gse_lat;
  double sun_lon = 0.0;
  double sun_lat = 0.0;
  double delta;
  double az;

  trans(GEO2GSE,
        fday,
        longitude,
        latitude,
        &gse_lon,
        &gse_lat);
  gse_lon *= 180.0 / M_PI;
  gse_lat *= 180.0 / M_PI;
  my_delaz(&gse_lat, &gse_lon, &sun_lat, &sun_lon, &delta, &az);

  return (delta * M_PI / 180.0);
} /* get_zenith() */
