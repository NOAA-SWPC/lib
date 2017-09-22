/*
 * plot_map.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <common/common.h>

#include "apex.h"
#include "iri.h"

int
map_fixed_alt(const char *filename, const time_t t, const double alt)
{
  int s = 0;
  const size_t nalt = 1;
  iri_workspace *w;
  double lon, lat;
  const double dlon = 5.0;
  const double dlat = 2.0;
  iri_result *result;
  FILE *fp;
  size_t i;
  apex_workspace *apex_p;
  int year = (int) get_year(t);

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "map_fixed_alt: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  apex_p = apex_alloc(year);

  w = iri_alloc(nalt, F107_IDX_FILE);
  result = iri_get_result(0, w);

  i = 1;
  fprintf(fp, "# Field %zu: longitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: electron density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: neutral temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: ion temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: electron temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: O+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: H+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: HE+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: O2+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: NO+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: N+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: cluster ion density (m^-3)\n", i++);

  for (lon = -180.0; lon <= 180.0; lon += dlon)
    {
      double phi = lon * M_PI / 180.0;

      for (lat = -90.0; lat <= 90.0; lat += dlat)
        {
          double theta = M_PI / 2.0 - lat * M_PI / 180.0;
          double alat, alon, qdlat;

          iri_calc(theta, phi, t, alt, 1.0, nalt, w);

          /* compute QD latitude */
          apex_transform_geodetic(theta, phi, alt * 1.0e3, &alon, &alat, &qdlat,
                                  NULL, NULL, NULL, apex_p);

          fprintf(fp, "%.3f %.3f %.3f %.3f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
                  lon,
                  lat,
                  alt,
                  qdlat,
                  result->Ne,
                  result->Tn,
                  result->Ti,
                  result->Te,
                  result->n_Op,
                  result->n_Hp,
                  result->n_HEp,
                  result->n_O2p,
                  result->n_NOp,
                  result->n_Np,
                  result->n_clus);
        }

      fprintf(fp, "\n");
    }

  iri_free(w);
  apex_free(apex_p);
  fclose(fp);

  return s;
} /* map_fixed_alt() */

int
map_fixed_alt_qd(const char *filename, const time_t t, const double alt)
{
  int s = 0;
  const size_t nalt = 1;
  iri_workspace *w;
  double lon, qdlat;
  const double dlon = 5.0;
  const double dlat = 2.0;
  iri_result *result;
  FILE *fp;
  size_t i;
  apex_workspace *apex_p;
  int year = (int) get_year(t);

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "map_fixed_alt: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  apex_p = apex_alloc(year);

  w = iri_alloc(nalt, F107_IDX_FILE);
  result = iri_get_result(0, w);

  i = 1;
  fprintf(fp, "# Field %zu: longitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: electron density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: neutral temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: ion temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: electron temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: O+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: H+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: HE+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: O2+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: NO+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: N+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: cluster ion density (m^-3)\n", i++);

  for (lon = -180.0; lon <= 180.0; lon += dlon)
    {
      for (qdlat = -90.0; qdlat <= 90.0; qdlat += dlat)
        {
          double theta;
          double glat, glon;

          apex_transform_inv_geodetic(qdlat, lon, alt * 1.0e3,
                                      &glat, &glon, apex_p);
          theta = M_PI / 2.0 - glat;

          iri_calc(theta, glon, t, alt, 1.0, nalt, w);

          fprintf(fp, "%.3f %.3f %.3f %.3f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
                  lon,
                  glat * 180.0 / M_PI,
                  alt,
                  qdlat,
                  result->Ne,
                  result->Tn,
                  result->Ti,
                  result->Te,
                  result->n_Op,
                  result->n_Hp,
                  result->n_HEp,
                  result->n_O2p,
                  result->n_NOp,
                  result->n_Np,
                  result->n_clus);
        }

      fprintf(fp, "\n");
    }

  iri_free(w);
  apex_free(apex_p);
  fclose(fp);

  return s;
} /* map_fixed_alt_qd() */

int
map_fixed_phi(const char *filename, const time_t t, const double phi)
{
  int s = 0;
  iri_workspace *w;
  const double altmin = 50.0;
  const double altstp = 2.0;
  const size_t nalt = 300;
  double lat;
  const double dlat = 2.0;
  FILE *fp;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "map_fixed_phi: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  w = iri_alloc(nalt, F107_IDX_FILE);

  i = 1;
  fprintf(fp, "# Field %zu: longitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: electron density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: neutral temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: ion temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: electron temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: O+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: H+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: HE+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: O2+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: NO+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: N+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: cluster ion density (m^-3)\n", i++);

  for (lat = -90.0; lat <= 90.0; lat += dlat)
    {
      double theta = M_PI / 2.0 - lat * M_PI / 180.0;

      iri_calc(theta, phi, t, altmin, altstp, nalt, w);

      for (i = 0; i < nalt; ++i)
        {
          double alt = altmin + i * altstp;
          iri_result *result = iri_get_result(i, w);

          fprintf(fp, "%.3f %.3f %.3f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
                  phi * 180 / M_PI,
                  lat,
                  alt,
                  result->Ne,
                  result->Tn,
                  result->Ti,
                  result->Te,
                  result->n_Op,
                  result->n_Hp,
                  result->n_HEp,
                  result->n_O2p,
                  result->n_NOp,
                  result->n_Np,
                  result->n_clus);
        }

      fprintf(fp, "\n");
    }

  iri_free(w);
  fclose(fp);

  return s;
} /* map_fixed_phi() */

int
plot_season(const char *filename, const time_t t0, const time_t t1,
            const double alt, const double theta, const double phi)
{
  int s = 0;
  iri_workspace *w;
  const size_t nalt = 1;
  const time_t dt = 86400;
  time_t t;
  FILE *fp;
  size_t i;
  iri_result *result;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "map_fixed_phi: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  w = iri_alloc(nalt, F107_IDX_FILE);
  result = iri_get_result(0, w);

  i = 1;
  fprintf(fp, "# Field %zu: timestamp\n", i++);
  fprintf(fp, "# Field %zu: time (decimal year)\n", i++);
  fprintf(fp, "# Field %zu: longitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: season (doy)\n", i++);
  fprintf(fp, "# Field %zu: electron density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: neutral temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: ion temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: electron temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: O+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: H+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: HE+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: O2+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: NO+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: N+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: cluster ion density (m^-3)\n", i++);

  for (t = t0; t <= t1; t += dt)
    {
      iri_calc(theta, phi, t, alt, 1.0, nalt, w);

      fprintf(fp, "%ld %f %.3f %.3f %.3f %.3f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
              t,
              get_year(t),
              phi * 180 / M_PI,
              90.0 - theta * 180.0 / M_PI,
              alt,
              get_season(t),
              result->Ne,
              result->Tn,
              result->Ti,
              result->Te,
              result->n_Op,
              result->n_Hp,
              result->n_HEp,
              result->n_O2p,
              result->n_NOp,
              result->n_Np,
              result->n_clus);
    }

  iri_free(w);
  fclose(fp);

  return s;
} /* plot_season() */

int
plot_lt(const char *filename, const time_t t0, const time_t t1,
        const double alt, const double theta, const double phi)
{
  int s = 0;
  iri_workspace *w;
  const size_t nalt = 1;
  const time_t dt = 86400*5;
  time_t t, tday;
  FILE *fp;
  size_t i;
  iri_result *result;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "map_fixed_phi: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  w = iri_alloc(nalt, F107_IDX_FILE);
  result = iri_get_result(0, w);

  i = 1;
  fprintf(fp, "# Field %zu: time (decimal year)\n", i++);
  fprintf(fp, "# Field %zu: local time (hours)\n", i++);
  fprintf(fp, "# Field %zu: longitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: season (doy)\n", i++);
  fprintf(fp, "# Field %zu: electron density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: neutral temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: ion temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: electron temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: O+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: H+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: HE+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: O2+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: NO+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: N+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: cluster ion density (m^-3)\n", i++);

  /* loop over doy (season) */
  for (tday = t0; tday <= t1; tday += dt)
    {
      /* loop over 1 day */
      for (t = tday; t < tday + 86400; t += 1800)
        {
          iri_calc(theta, phi, t, alt, 1.0, nalt, w);

          fprintf(fp, "%f %.2f %.3f %.3f %.3f %.3f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
                  get_year(t),
                  get_localtime(t, phi),
                  phi * 180 / M_PI,
                  90.0 - theta * 180.0 / M_PI,
                  alt,
                  get_season(t),
                  result->Ne,
                  result->Tn,
                  result->Ti,
                  result->Te,
                  result->n_Op,
                  result->n_Hp,
                  result->n_HEp,
                  result->n_O2p,
                  result->n_NOp,
                  result->n_Np,
                  result->n_clus);
        }
    }

  iri_free(w);
  fclose(fp);

  return s;
} /* plot_lt() */

int
main(void)
{
  time_t t;
  const char *mapfile = "map.dat";

#if 1
  t = 955303200; /* Apr 9 2000 18:00:00 UTC */
  t = 955281600; /* Apr 9 2000 12:00:00 UTC */
#else
  t = 961610400; /* Jun 21 2000 18:00:00 UTC */
#endif

#if 1
  /* fixed altitude, lat/lon map */
  {
    double alt = 400.0;
    fprintf(stderr, "main: writing map to %s...", mapfile);
    map_fixed_alt(mapfile, t, alt);
    fprintf(stderr, "done\n");
  }
#elif 1
  /* fixed altitude, qdlat/lon map */
  {
    double alt = 400.0;
    fprintf(stderr, "main: writing map to %s...", mapfile);
    map_fixed_alt_qd(mapfile, t, alt);
    fprintf(stderr, "done\n");
  }
#elif 0
  /* fixed lon, alt/lat map */
  {
    double phi = -50.0 * M_PI / 180.0;
    fprintf(stderr, "main: writing map to %s...", mapfile);
    map_fixed_phi(mapfile, t, phi);
    fprintf(stderr, "done\n");
  }
#elif 0
  /* fixed position, time/season map */
  {
    double alt = 400.0;
    double phi = -50.0 * M_PI / 180.0;
    double lat = -10.0 * M_PI / 180.0;
    double theta = M_PI / 2.0 - lat;
    const time_t t0 = 946749600;  /* Jan 1 2000 18:00:00 UTC */
    const time_t t1 = 1420135200; /* Jan 1 2015 18:00:00 UTC */
    const char *outfile = "season.dat";

    fprintf(stderr, "main: writing season data to %s...", outfile);
    plot_season(outfile, t0, t1, alt, theta, phi);
    fprintf(stderr, "done\n");
  }
#elif 1
  /* fixed position, local time/season map */
  {
    double alt = 400.0;
    double phi = -50.0 * M_PI / 180.0;
    double lat = -10.0 * M_PI / 180.0;
    double theta = M_PI / 2.0 - lat;
    const time_t t0 = 946684800;  /* Jan 1 2000 00:00:00 UTC */
    const time_t t1 = 978307199; /* Dec 31 2000 23:59:59 UTC */
    const char *outfile = "lt.dat";

    fprintf(stderr, "main: writing LT data to %s...", outfile);
    plot_lt(outfile, t0, t1, alt, theta, phi);
    fprintf(stderr, "done\n");
  }
#endif

  return 0;
} /* main() */
