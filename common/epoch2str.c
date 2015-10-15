#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <satdata/satdata.h>

#include "common.h"

int
main(int argc, char *argv[])
{
  double epoch = 0.0;
  time_t t;

  if (argc > 1)
    epoch = atof(argv[1]);

  t = satdata_epoch2timet(epoch);

  putenv("TZ=GMT");
  printf("%s", ctime(&t));
    
  return 0;
} /* main() */
