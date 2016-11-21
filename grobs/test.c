#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "grobs.h"

int
main(int argc, char *argv[])
{
  int c;
  grobs_data *data = NULL;
  size_t i;

  while ((c = getopt(argc, argv, "i:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            data = iaga_read_HDZF(optarg, NULL);
            break;

          case '?':
          default:
            printf("usage: %s [options]\n", argv[0]);
            exit(1);
            break;
        }
    }

  if (data == NULL)
    {
      printf("usage: %s <-i iaga_file>\n", argv[0]);
      exit(1);
    }

  for (i = 0; i < data->n; ++i)
    {
      printf("%ld %f %f %f\n",
             data->t[i],
             data->X[i],
             data->Y[i],
             data->Z[i]);
    }

  return 0;
}
