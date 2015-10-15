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
  struct tm *tm_p;

  if (argc > 1)
    t = atol(argv[1]);

  putenv("TZ=GMT");

  tm_p = gmtime(&t);

  sprintf(buf, "%04d%02d%02d_%02d%02d%02d",
          tm_p->tm_year + 1900,
          tm_p->tm_mon + 1,
          tm_p->tm_mday,
          tm_p->tm_hour,
          tm_p->tm_min,
          tm_p->tm_sec);
  printf("%s", buf);
    
  return 0;
} /* main() */
