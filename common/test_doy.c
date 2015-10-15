#include <stdio.h>
#include <stdlib.h>

#include "common.h"

int
main(int argc, char *argv[])
{
  int year, doy, month, day;

  year = 2000;
  doy = 31 + 29;

  if (argc > 1)
    year = atoi(argv[1]);
  if (argc > 2)
    doy = atoi(argv[2]);

  doy2md(year, doy, &month, &day);

  printf("%d/%d ==> %d/%d/%d\n", year, doy, month, day, year);

  return 0;
} /* main() */
