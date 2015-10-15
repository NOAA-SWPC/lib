#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_math.h>

#include "common.h"

int
main(int argc, char *argv[])
{
  time_t t = 0;
  char buf[2048];

  if (argc > 1)
    t = atol(argv[1]);

  putenv("TZ=GMT");

  sprintf(buf, "%s", ctime(&t));
  buf[strlen(buf) - 1] = '\0';
  printf("%s", buf);
    
  return 0;
} /* main() */
