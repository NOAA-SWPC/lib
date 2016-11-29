/*
 * grobs.c
 * Routines for reading ground observatory data
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "grobs.h"

grobs_data *
grobs_alloc(const size_t ntot)
{
  grobs_data *data;

  data = calloc(1, sizeof(grobs_data));
  if (!data)
    return 0;

  grobs_realloc(ntot, data);

  data->n = 0;

  return data;
}

grobs_data *
grobs_realloc(const size_t ntot, grobs_data *data)
{
  size_t i;

  if (!data)
    {
      return grobs_alloc(ntot);
    }

  if (ntot <= data->ntot)
    return data; /* nothing to do */

  /*
   * If data is already large (ie: 1 year of data at 1Hz or more),
   * warn the user in case they want to use a more efficient allocation
   * strategy
   */
  if (data->ntot >= 31536000)
    {
      fprintf(stderr, "grobs_realloc: warning: reallocing large structure: %zu ==> %zu\n",
              data->ntot, ntot);
    }

  data->t = realloc(data->t, ntot * sizeof(time_t));
  data->X = realloc(data->X, ntot * sizeof(double));
  data->Y = realloc(data->Y, ntot * sizeof(double));
  data->Z = realloc(data->Z, ntot * sizeof(double));
  data->H = realloc(data->H, ntot * sizeof(double));
  data->D = realloc(data->D, ntot * sizeof(double));
  data->I = realloc(data->I, ntot * sizeof(double));

  for (i = 0; i < ntot; ++i)
    {
      data->t[i] = 0;
      data->X[i] = 0.0;
      data->Y[i] = 0.0;
      data->Z[i] = 0.0;
      data->H[i] = 0.0;
      data->D[i] = 0.0;
      data->I[i] = 0.0;
    }

  data->ntot = ntot;

  return data;
}

void
grobs_free(grobs_data *data)
{
  if (data->t)
    free(data->t);

  if (data->X)
    free(data->X);

  if (data->Y)
    free(data->Y);

  if (data->Z)
    free(data->Z);

  if (data->H)
    free(data->H);

  if (data->D)
    free(data->D);

  if (data->I)
    free(data->I);

  free(data);
}
