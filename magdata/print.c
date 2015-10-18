/*
 * print.c
 *
 * Print contents of magdata file
 *
 * ./print <-i magdata_file>
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <assert.h>
#include <errno.h>

#include "common.h"
#include "magdata.h"

int
main(int argc, char *argv[])
{
  magdata *data = NULL;
  char *outfile = "data.dat";
  struct timeval tv0, tv1;
  int c;

  while ((c = getopt(argc, argv, "i:o:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);

            data = magdata_read(optarg, NULL);
            if (!data)
              {
                fprintf(stderr, "main: error reading %s\n", optarg);
                exit(1);
              }

            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
                    time_diff(tv0, tv1));
            break;

          case 'o':
            outfile = optarg;
            break;

          default:
            break;
        }
    }

  if (!data)
    {
      fprintf(stderr, "Usage: %s <-i magdata_file> [-o output_file]\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: writing output to %s...", outfile);
  magdata_print(outfile, data);
  fprintf(stderr, "done\n");

  magdata_free(data);

  return 0;
} /* main() */
