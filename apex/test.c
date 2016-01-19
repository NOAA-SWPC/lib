/*
 * test.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include "apex.h"

int
main()
{
  const double R = 6371.2e3;
  apex_workspace *apex_p;
  double qdlat, qdlon;
  double lat, lon, alt;

  apex_p = apex_alloc(2005);

#if 1
  qdlat = 0.0;
  alt = 110.0e3;
  for (qdlon = -180.0; qdlon < 180.0; qdlon += 5.0)
    {
      double theta, alon, alat, qdlat2;

      apex_transform_inv(qdlat, qdlon, alt, &lat, &lon, apex_p);

      theta = M_PI / 2.0 - lat;
      apex_transform(theta, lon, alt + R, &alon, &alat, &qdlat2,
                     NULL, NULL, NULL, apex_p);

      gsl_test_abs(qdlat2, qdlat, 1.0e-1, "qdlat geocentric");
      gsl_test_rel(alon, qdlon, 1.0e-3, "qdlon geocentric");

      apex_transform_inv_geodetic(qdlat, qdlon, alt, &lat, &lon, apex_p);

      theta = M_PI / 2.0 - lat;
      apex_transform_geodetic(theta, lon, alt, &alon, &alat, &qdlat2,
                              NULL, NULL, NULL, apex_p);

      gsl_test_abs(qdlat2, qdlat, 1.0e-3, "qdlat geodetic");
      gsl_test_rel(alon, qdlon, 1.0e-5, "qdlon geodetic");

      printf("%f %f\n",
             lon * 180.0 / M_PI,
             lat * 180.0 / M_PI);
    }
#endif

  apex_free(apex_p);

  exit (gsl_test_summary());
} /* main() */
