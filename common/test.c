#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "common.h"

int
main(int argc, char *argv[])
{
  int c;
  double phi;
  double t;
  int year;
  int season;
  double fday;
  double lt;
  time_t ts;

  phi = 90.0 * M_PI / 180.0;
  t = 12.0;
  year = 2002;
  season = 79;

  while ((c = getopt(argc, argv, "p:t:s:y:")) != (-1))
    {
      switch (c)
        {
          case 'p':
            phi = strtod(optarg, NULL) * M_PI / 180.0;
            break;

          case 't':
            t = strtod(optarg, NULL);
            break;

          case 's':
            season = strtol(optarg, NULL, 0);
            break;

          case 'y':
            year = strtol(optarg, NULL, 0);
            break;

          case '?':
          default:
            printf("usage: %s [options]\n", argv[0]);
            exit(1);
            break;
        }
    }

  /*fday = get_fday(year, phi, t, season);*/

  ts = fday2timet(2521.97322);
  printf("ts = %d, %s\n", ts, ctime(&ts));
    
  return 0;
} /* main() */
