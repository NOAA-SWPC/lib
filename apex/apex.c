/*
 * apex.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <common/geo.h>

#include "apex.h"

void makeapxsh_(char *datafilein, float *epochgridin, int *nepochin,
                int *lmaxin, int *mmaxin, int *nmaxin);
void loadapxsh_(char *datafilein, float *epoch);
void apxg2all_(float *glat, float *glon, float *alt, float *hr,
               int *vecflag, float *qlat, float *qlon, float *mlat,
               float *mlon, float *f1, float *f2, float *f,
               float *d1, float *d2, float *d3, float *d,
               float *e1, float *e2, float *e3);
void apxq2g_(float *qlat, float *qlon, float *alt, float *prec,
             float *glat, float *glon, float *error);

apex_workspace *
apex_alloc(int year)
{
  apex_workspace *w;
  int s;
  char filename[8192];
  const size_t start_year = 2000;
  const size_t end_year = 2020;
  size_t i;

  w = calloc(1, sizeof(apex_workspace));
  if (!w)
    return 0;

  w->nepochs = end_year - start_year + 1;
  w->epochgrid = malloc(w->nepochs * sizeof(float));

  for (i = 0; i < w->nepochs; ++i)
    w->epochgrid[i] = (float) (start_year + i);

  w->ntheta = APEX_NLAT;
  w->nphi = APEX_NLON;
  w->nalt = APEX_NALT;

  w->lwk = w->ntheta * w->nphi * w->nalt * 5 + w->ntheta + w->nphi + w->nalt;
  w->wk = malloc(w->lwk * sizeof(float));

  sprintf(filename, "%s", APEX_DATAFILE);

  s = access(filename, F_OK);
  if (s)
    {
      /* file does not exist */
      fprintf(stderr, "apex_alloc: creating apex grid file %s...", filename);
      s = apex_makefile(filename, w);
      fprintf(stderr, "done (s = %d)\n", s);
    }

  fprintf(stderr, "apex_alloc: reading apex grid file %s...", filename);
  s = apex_readfile(filename, year, w);
  fprintf(stderr, "done (s = %d)\n", s);

  return w;
} /* apex_alloc() */

void
apex_free(apex_workspace *w)
{
  if (w->epochgrid)
    free(w->epochgrid);

  if (w->wk)
    free(w->wk);

  free(w);
} /* apex_free() */

/*
apex_makefile()
  Make an apex grid for several epochs and write it to a file

Inputs: filename - file to write
        w        - workspace
*/

int
apex_makefile(const char *filename, apex_workspace *w)
{
  int s = 0;
  int lmax = 5;
  int mmax = 12;
  int nmax = 12;

  makeapxsh_((char *) filename,
             w->epochgrid,
             &(w->nepochs),
             &lmax,
             &mmax,
             &nmax);

  return s;
} /* apex_makefile() */

int
apex_readfile(const char *filename, int year, apex_workspace *w)
{
  int s = 0;
  float epoch = (float) year;

  loadapxsh_((char *) filename, &epoch);

  return s;
} /* apex_readfile() */

/*
apex_transform_geodetic()
  Transform from geodetic to apex coordinates

Inputs: theta   - geodetic colatitude (radians)
        phi     - geodetic longitude (radians)
        alt     - geodetic altitude (m)
        apexlon - (output) apex longitude in degrees
        apexlat - (output) apex latitude in degrees
        qdlat   - (output) quasi-dipole latitude in degrees
        E1      - (output) components (east, north, up) of E1 vector
        E2      - (output) components (east, north, up) of E2 vector
        E3      - (output) components (east, north, up) of E3 vector
        w       - workspace

Notes:
1) Workspace w is not used in this function, but apex_alloc() MUST be called
to load coefficients for appropriate epoch prior to calling this function
*/

int
apex_transform_geodetic(double theta, double phi, double alt,
                        double *apexlon, double *apexlat, double *qdlat,
                        double *E1, double *E2, double *E3, apex_workspace *w)
{
  int s = 0;
  size_t i;
  float lat = 90.0 - theta * 180.0 / M_PI;
  float lon = phi * 180.0 / M_PI;
  float falt = alt * 1.0e-3;
  float hr = 70.0;
  float d1[3], d2[3], d3[3], d[3];
  float e1[3], e2[3], e3[3];
  float f1[2], f2[2], f[3];
  float qlat, qlon, mlat, mlon;
  int vecflag = 1;
  float glon, glat;

  if (lon > 180.0)
    lon -= 360.0;
  if (lon < -180.0)
    lon += 360.0;

  assert(lat >= -90.0 && lat <= 90.0);

  glon = lon;
  glat = lat;

  apxg2all_(&glat,
            &glon,
            &falt,
            &hr,
            &vecflag,
            &qlat,
            &qlon,
            &mlat,
            &mlon,
            f1,
            f2,
            f,
            d1,
            d2,
            d3,
            d,
            e1,
            e2,
            e3);

  *apexlon = (double) mlon;
  *apexlat = (double) mlat;
  *qdlat = (double) qlat;

  for (i = 0; i < 3; ++i)
    {
      if (E1)
        E1[i] = e1[i];

      if (E2)
        E2[i] = e2[i];

      if (E3)
        E3[i] = e3[i];
    }

  return s;
} /* apex_transform_geodetic() */

/*
apex_transform()
  Transform from geocentric to apex coordinates

Inputs: theta   - geocentric colatitude (radians)
        phi     - geocentric longitude (radians)
        r       - geocentric radius (m)
        apexlon - (output) apex longitude in degrees
        apexlat - (output) apex latitude in degrees
        qdlat   - (output) quasi-dipole latitude in degrees
        E1      - (output) components (east, north, up) of E1 vector
        E2      - (output) components (east, north, up) of E2 vector
        E3      - (output) components (east, north, up) of E3 vector
        w       - workspace

Notes:
1) Workspace w is not used in this function, but apex_alloc() MUST be called
to load coefficients for appropriate epoch prior to calling this function
*/

int
apex_transform(double theta, double phi, double r,
               double *apexlon, double *apexlat, double *qdlat,
               double *E1, double *E2, double *E3, apex_workspace *w)
{
  int s = 0;
  double lat = M_PI / 2.0 - theta;
  double latd, altd, thetad;

  geo2geodetic(lat, phi, r * 1.0e-3, &latd, &altd);

  thetad = M_PI / 2.0 - latd;

  s = apex_transform_geodetic(thetad, phi, altd * 1.0e3,
                              apexlon, apexlat, qdlat,
                              E1, E2, E3, w);

  return s;
} /* apex_transform() */

/*
apex_transform_inv_geodetic()
  Perform the inverse transform (magnetic to geodetic coordinates)

Inputs: qdlat - qd latitude in degrees
        qdlon - qd longitude in degrees
        alt   - geodetic altitude in m
        glat  - geodetic latitude (radians)
        glon  - geodetic longitude (radians)
*/

int
apex_transform_inv_geodetic(double qdlat, double qdlon, double alt,
                            double *glat, double *glon,
                            apex_workspace *w)
{
  int s = 0;
  float fqdlat = qdlat;
  float fqdlon = qdlon;
  float falt = alt * 1.0e-3;
  float fglat, fglon;
  float prec = 1.0e-8;
  float err;

  apxq2g_(&fqdlat,
          &fqdlon,
          &falt,
          &prec,
          &fglat,
          &fglon,
          &err);

  *glat = fglat * M_PI / 180.0;
  *glon = fglon * M_PI / 180.0;

  return s;
} /* apex_transform_inv_geodetic() */

/*
apex_transform_inv()
  Perform the inverse transform (magnetic to geocentric coordinates)

Inputs: qdlat - qd latitude in degrees
        qdlon - qd longitude in degrees
        alt   - geodetic altitude in m
        latc  - geocentric latitude (radians)
        lon   - geocentric longitude (radians)
*/

int
apex_transform_inv(double qdlat, double qdlon, double alt,
                   double *latc, double *lon,
                   apex_workspace *w)
{
  int s = 0;
  double latd; /* geodetic latitude */
  double r;    /* geocentric radius */

  s = apex_transform_inv_geodetic(qdlat, qdlon, alt, &latd, lon, w);

  /* convert to geocentric */
  geodetic2geo(latd, alt * 1.0e-3, latc, &r);

  return s;
} /* apex_transform_inv() */
