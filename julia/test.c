/*
 * test.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <common/common.h>
#include <common/geo.h>

#include "julia.h"

int
main(int argc, char *argv[])
{
  julia_data *data;
  size_t i;

  /* compute geocentric latitude of JULIA */
  {
    double latd = JULIA_GEODETIC_LAT * M_PI / 180.0;
    double latc, r;

    geodetic2geo(latd, 0.0, &latc, &r);
    fprintf(stderr, "main: JULIA geocentric latitude = %f [deg]\n",
            latc * 180.0 / M_PI);
    fprintf(stderr, "main: JULIA geocentric radius = %f [km]\n", r);
  }

  fprintf(stderr, "main: reading all avg JULIA files...");
  data = julia_read_avg_idx(JULIA_AVG_IDX_FILE, NULL);
  fprintf(stderr, "done (%zu data read)\n", data->n);

  fprintf(stderr, "main: reading all non-avg JULIA files...");
  julia_read_idx(JULIA_IDX_FILE, data);
  fprintf(stderr, "done (%zu data total)\n", data->n);

  fprintf(stderr, "main: sorting JULIA data...");
  julia_sort(data);
  fprintf(stderr, "done\n");

  for (i = 0; i < data->n; ++i)
    {
      double lt = get_localtime(data->array[i].t, JULIA_LON_RAD);
      time_t t = data->array[i].t;

      if (t < 1397692800)
        continue;

      printf("%ld %f %f %f\n",
             data->array[i].t,
             lt,
             data->array[i].v_zonal,
             data->array[i].v_vert);
    }

  julia_free(data);

  return 0;
} /* main() */
