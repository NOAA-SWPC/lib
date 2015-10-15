#include <stdio.h>
#include <stdlib.h>

#include "common.h"

int
main(int argc, char *argv[])
{
  int year, doy, month, day;

  if (argc < 3)
    {
      fprintf(stderr, "Usage: %s <year> <doy>\n", argv[0]);
      exit(1);
    }

  year = atoi(argv[1]);
  doy = atoi(argv[2]);

  doy2md(year, doy, &month, &day);

  printf("%02d %02d\n", month, day);

  return 0;
} /* main() */
