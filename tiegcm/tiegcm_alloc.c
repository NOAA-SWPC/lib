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
tiegcm_alloc(const size_t n)
{
  tiegcm_data *data;

  data = calloc(1, sizeof(tiegcm_data));
  if (!data)
    {
      fprintf(stderr, "tiegcm_alloc: calloc failed: %s\n",
              strerror(errno));
      return 0;
    }

  tiegcm_realloc(n, data);

  data->n = 0;

  return data;
}

tiegcm_data *
tiegcm_realloc(const size_t n, tiegcm_data *data)
{
  size_t i;

  if (!data)
    {
      return tiegcm_alloc(n);
    }

  if (n <= data->ntot)
    return data; /* nothing to do */

  data->t = realloc(data->t, n * sizeof(time_t));
  data->year = realloc(data->year, n * sizeof(double));
  data->doy = realloc(data->doy, n * sizeof(double));
  data->ut = realloc(data->ut, n * sizeof(double));

  /* initialize newly allocated memory */
  for (i = data->ntot; i < n; ++i)
    {
      data->t[i] = 0;
      data->year[i] = 0.0;
      data->doy[i] = 0.0;
      data->ut[i] = 0.0;
    }

  data->ntot = n;

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

  free(data);
}
