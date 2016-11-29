/*
 * wamnet.c
 *
 * Routines for reading IAGA formatted data files
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "grobs.h"

grobs_data *
grobs_wamnet_read(const char *filename, grobs_data *data)
{
  FILE *fp;
  char buf[GROBS_MAX_BUFFER];
  size_t n;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "grobs_wamnet_read: fopen: cannot open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  /* allocate space for 1-minute data */
  if (data == NULL)
    data = grobs_alloc(GROBS_MAX_YEAR * 527040);

  n = data->n;

  while (fgets(buf, sizeof(buf), fp) != NULL)
    {
      int c;
      double fday, dH;

      c = sscanf(buf, "%lf %lf", &fday, &dH);
      if (c < 2)
        continue;

      if (!isfinite(dH))
        continue;

      data->t[n] = fday2timet(fday);
      data->H[n] = dH;

      if (++n >= data->ntot)
        {
          fprintf(stderr, "grobs_wamnet_read: error: ntot too small\n");
          exit(1);
        }
    }

  fclose(fp);

  data->n = n;
  data->glon = -5.625; /* SAM longitude */
  data->glat = 11.394; /* SAM latitude */

  return data;
}
