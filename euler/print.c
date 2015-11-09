/*
 * print.c
 *
 * Print Euler angles in different formats
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "magdata.h"

#include "euler.h"

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --input_file | -i file         - Euler angle input file\n");
  fprintf(stderr, "\t --swarm_file | -s file         - Swarm formatted output file\n");
  fprintf(stderr, "\t --start | -a md2000_start      - Starting MJD2000 day\n");
  fprintf(stderr, "\t --end | -b md2000_end          - Ending MJD2000 day\n");
} /* print_help() */

int
main(int argc, char *argv[])
{
  char *infile = NULL;
  char *swarm_file = NULL;
  struct timeval tv0, tv1;
  euler_workspace *euler_workspace_p;
  double fday_start = 5090.0;
  double fday_end = 5715.0;
  double fday_step = 10.0;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "start", required_argument, NULL, 'a' },
          { "end", required_argument, NULL, 'b' },
          { "input_file", required_argument, NULL, 'i' },
          { "swarm_file", required_argument, NULL, 's' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:b:i:s:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            fday_start = atof(optarg);
            break;

          case 'b':
            fday_end = atof(optarg);
            break;

          case 'i':
            infile = optarg;
            break;

          case 's':
            swarm_file = optarg;
            break;

          default:
            print_help(argv);
            exit(1);
            break;
        }
    }

  if (infile == NULL)
    {
      print_help(argv);
      exit(1);
    }

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);
  euler_workspace_p = euler_read(infile);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu angles read, %g seconds)\n",
          euler_workspace_p->n, time_diff(tv0, tv1));

  if (swarm_file)
    {
      fprintf(stderr, "main: writing Euler angles to %s...", swarm_file);
      euler_write_swarm(fday_start, fday_end, fday_step, swarm_file, euler_workspace_p);
      fprintf(stderr, "done\n");
    }

  return 0;
}
