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
  double fday = 0.0;
  time_t t;

  if (argc > 1)
    fday = atof(argv[1]);

  t = fday2timet(fday);

  putenv("TZ=GMT");
  printf("%s", ctime(&t));
    
  return 0;
} /* main() */
