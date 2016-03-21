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

#include "tiegcm.h"

tiegcm_data *
tiegcm_alloc(const size_t nt, const size_t nlon, const size_t nlat)
{
  tiegcm_data *data;

  data = calloc(1, sizeof(tiegcm_data));
  if (!data)
    {
      fprintf(stderr, "tiegcm_alloc: calloc failed: %s\n",
              strerror(errno));
      return 0;
    }

  tiegcm_realloc(nt, nlon, nlat, data);

  data->nt_max = 2000;
  data->workx = malloc(data->nt_max * nlon * nlat * sizeof(double));
  data->worky = malloc(data->nt_max * nlon * nlat * sizeof(double));
  data->workz = malloc(data->nt_max * nlon * nlat * sizeof(double));

  return data;
}

tiegcm_data *
tiegcm_realloc(const size_t nt, const size_t nlon, const size_t nlat,
               tiegcm_data *data)
{
  size_t i;

  if (!data)
    {
      return tiegcm_alloc(nt, nlon, nlat);
    }

  if (nt <= data->ntot)
    return data; /* nothing to do */

  data->t = realloc(data->t, nt * sizeof(time_t));
  data->year = realloc(data->year, nt * sizeof(double));
  data->doy = realloc(data->doy, nt * sizeof(double));
  data->ut = realloc(data->ut, nt * sizeof(double));

  data->glon = realloc(data->glon, nlon * sizeof(double));
  data->glat = realloc(data->glat, nlat * sizeof(double));

  data->Bx = realloc(data->Bx, nt * nlon * nlat * sizeof(double));
  data->By = realloc(data->By, nt * nlon * nlat * sizeof(double));
  data->Bz = realloc(data->Bz, nt * nlon * nlat * sizeof(double));

  /* initialize newly allocated memory */
  for (i = data->ntot; i < nt; ++i)
    {
      data->t[i] = 0;
      data->year[i] = 0.0;
      data->doy[i] = 0.0;
      data->ut[i] = 0.0;
    }

  data->ntot = nt;
  data->nlon = nlon;
  data->nlat = nlat;

  return data;
}

void
tiegcm_free(tiegcm_data *data)
{
  if (data->t)
    free(data->t);

  if (data->year)
    free(data->year);

  if (data->doy)
    free(data->doy);

  if (data->ut)
    free(data->ut);

  if (data->glon)
    free(data->glon);

  if (data->glat)
    free(data->glat);

  if (data->Bx)
    free(data->Bx);

  if (data->By)
    free(data->By);

  if (data->Bz)
    free(data->Bz);

  if (data->workx)
    free(data->workx);

  if (data->worky)
    free(data->worky);

  if (data->workz)
    free(data->workz);

  free(data);
}
