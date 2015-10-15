/*
 * test.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "common.h"
#include "magdata.h"

int
main()
{
  magdata *data;

  data = magdata_alloc(100, R_EARTH_KM);

  magdata_free(data);

  return 0;
}
