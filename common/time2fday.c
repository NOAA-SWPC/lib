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
  double fday;
  time_t t = 0;

  if (argc > 1)
    t = atol(argv[1]);

  fday = time2fday(t);
  printf("fday = %f\n", fday);
    
  return 0;
} /* main() */
