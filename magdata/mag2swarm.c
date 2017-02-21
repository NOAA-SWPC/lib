/*
 * mag2swarm.c
 *
 * Convert a magdata file to Swarm CDF format
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <satdata/satdata.h>

#include "common.h"
#include "magdata.h"

int
main(int argc, char *argv[])
{
  char *swarm_file = NULL;
  magdata *mdata = NULL;
  satdata_mag *data;
  int c;

  while ((c = getopt(argc, argv, "i:o:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            mdata = magdata_read(optarg, NULL);
            fprintf(stderr, "done (%zu data read)\n", mdata->n);
            break;

          case 'o':
            swarm_file = optarg;
            break;

          default:
            break;
        }
    }

  if (mdata == NULL || swarm_file == NULL)
    {
      fprintf(stderr, "Usage: %s <-i magdata_file> <-o swarm_file>\n", argv[0]);
      exit(1);
    }

  data = magdata_mag2sat(mdata);

  fprintf(stderr, "main: writing Swarm CDF file %s...", swarm_file);
  satdata_swarm_write(0, swarm_file, data);
  fprintf(stderr, "done (%zu data written)\n", data->n);

  magdata_free(mdata);
  satdata_mag_free(data);

  return 0;
}
