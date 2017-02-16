/*
 * print_eph.c
 *
 * Print Swarm ephemeris data
 *
 * ./print_eph <-i swarm_index_file>
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

#include <satdata/satdata.h>

#include "common.h"

/*
print_data()

Inputs: down_sample - number of samples to throw out (>= 1)
                      (ie: if this is 5, every 5th sample is kept and
                       the rest discarded)
        data        - satellite data input
*/

int
print_data(const int down_sample, const satdata_mag *data)
{
  int s = 0;
  size_t i;

  i = 1;
  printf("# Field %zu: year\n", i++);
  printf("# Field %zu: month\n", i++);
  printf("# Field %zu: day\n", i++);
  printf("# Field %zu: hour\n", i++);
  printf("# Field %zu: minute\n", i++);
  printf("# Field %zu: second\n", i++);
  printf("# Field %zu: longitude (degrees)\n", i++);
  printf("# Field %zu: latitude (degrees)\n", i++);
  printf("# Field %zu: geocentric radius (km)\n", i++);

  for (i = 0; i < data->n; i += down_sample)
    {
      long year, month, day, hour;
      long min, sec, msec;

      EPOCHbreakdown(data->t[i], &year, &month, &day, &hour, &min, &sec, &msec);

      printf("%04ld %02ld %02ld %02ld %02ld %02ld %9.4f %8.4f %9.4f\n",
             year,
             month,
             day,
             hour,
             min,
             sec,
             data->longitude[i],
             data->latitude[i],
             data->r[i]);
    }

  return s;
}

int
main(int argc, char *argv[])
{
  satdata_mag *data;
  struct timeval tv0, tv1;
  int c;
  char *infile = NULL;
  int down_sample = 1;

  while ((c = getopt(argc, argv, "i:d:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'd':
            down_sample = atoi(optarg);
            break;

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i swarm_index_file> [-d down_sample]\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);
  fprintf(stderr, "downsample factor = %d\n", down_sample);

  fprintf(stderr, "Reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = satdata_swarm_read_idx(infile, 0);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
          time_diff(tv0, tv1));

  print_data(down_sample, data);

  satdata_mag_free(data);

  return 0;
}
