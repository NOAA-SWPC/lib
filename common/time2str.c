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
  time_t t = 0;

  if (argc > 1)
    t = atol(argv[1]);

  putenv("TZ=GMT");
  printf("year = %f\n", get_year(t));
  printf("%s", ctime(&t));
    
  return 0;
} /* main() */
