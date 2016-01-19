/*
 * print.c
 *
 * Print contents of tiegcm data files
 *
 * ./print <-i tiegcm_nc_file>
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

#include <gsl/gsl_math.h>

#include "common.h"
#include "tiegcm.h"

/*
print_data()

Inputs: data - tiegcm data
*/

int
print_data(const tiegcm_data *data)
{
  int s = 0;
  size_t i;

  i = 1;
  printf("# Field %zu: time (decimal year)\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      printf("%ld\n",
             data->t[i]);
    }

  return s;
}

int
main(int argc, char *argv[])
{
  tiegcm_data *data;
  struct timeval tv0, tv1;
  char *infile = NULL;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i tiegcm_nc_file>\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = tiegcm_read(infile, NULL);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
          time_diff(tv0, tv1));

  print_data(data);

  tiegcm_free(data);

  return 0;
}
