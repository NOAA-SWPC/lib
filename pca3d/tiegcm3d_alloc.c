/*
 * tiegcm_alloc.c
 *
 * Routines for allocating tiegcm data structures
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cdf.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>

#include "tiegcm3d.h"

tiegcm3d_data *
tiegcm3d_alloc(const size_t nt, const size_t nr, const size_t nlon, const size_t nlat)
{
  tiegcm3d_data *data;

  data = calloc(1, sizeof(tiegcm3d_data));
  if (!data)
    {
      fprintf(stderr, "tiegcm3d_alloc: calloc failed: %s\n",
              strerror(errno));
      return 0;
    }

  tiegcm3d_realloc(nt, nr, nlon, nlat, data);

  data->nt_max = 9000;
  data->workx = malloc(data->nt_max * nlon * nlat * sizeof(double));
  data->worky = malloc(data->nt_max * nlon * nlat * sizeof(double));
  data->workz = malloc(data->nt_max * nlon * nlat * sizeof(double));

  return data;
}

tiegcm3d_data *
tiegcm3d_realloc(const size_t nt, const size_t nr, const size_t nlon, const size_t nlat,
                 tiegcm3d_data *data)
{
  size_t i;

  if (!data)
    {
      return tiegcm3d_alloc(nt, nr, nlon, nlat);
    }

  if (nt <= data->ntot)
    return data; /* nothing to do */

  data->t = realloc(data->t, nt * sizeof(time_t));
  data->year = realloc(data->year, nt * sizeof(double));
  data->doy = realloc(data->doy, nt * sizeof(double));
  data->ut = realloc(data->ut, nt * sizeof(double));

  data->r = realloc(data->r, nr * sizeof(double));
  data->glon = realloc(data->glon, nlon * sizeof(double));
  data->glat = realloc(data->glat, nlat * sizeof(double));

  data->Jr = realloc(data->Jr, nt * nr * nlon * nlat * sizeof(double));
  data->Jt = realloc(data->Jt, nt * nr * nlon * nlat * sizeof(double));
  data->Jp = realloc(data->Jp, nt * nr * nlon * nlat * sizeof(double));

  /* initialize newly allocated memory */
  for (i = data->ntot; i < nt; ++i)
    {
      data->t[i] = 0;
      data->year[i] = 0.0;
      data->doy[i] = 0.0;
      data->ut[i] = 0.0;
    }

  data->ntot = nt;
  data->nr = nr;
  data->nlon = nlon;
  data->nlat = nlat;

  return data;
}

void
tiegcm3d_free(tiegcm3d_data *data)
{
  if (data->t)
    free(data->t);

  if (data->year)
    free(data->year);

  if (data->doy)
    free(data->doy);

  if (data->ut)
    free(data->ut);

  if (data->r)
    free(data->r);

  if (data->glon)
    free(data->glon);

  if (data->glat)
    free(data->glat);

  if (data->Jr)
    free(data->Jr);

  if (data->Jt)
    free(data->Jt);

  if (data->Jp)
    free(data->Jp);

  if (data->workx)
    free(data->workx);

  if (data->worky)
    free(data->worky);

  if (data->workz)
    free(data->workz);

  free(data);
}
