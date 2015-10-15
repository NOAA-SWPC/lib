/*
 * test_geo.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include "geo.h"

int
main()
{
  double latc, rc;
  double lon;

  rc = 100.0 + 6371.2;
  lon = 20.0;

  for (latc = -90.0; latc < 90.0; latc += 1.0)
    {
      double latc_rad = latc * M_PI / 180.0;
      double latd, altd;
      double lat, r;

      geo2geodetic(latc_rad,
                   lon * M_PI / 180.0,
                   rc,
                   &latd,
                   &altd);

       geodetic2geo(latd, altd, &lat, &r);

       gsl_test_rel(lat, latc_rad, 1.0e-7, "lat");
       gsl_test_rel(r, rc, 1.0e-7, "r");
    }

  exit (gsl_test_summary());
}
