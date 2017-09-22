/*
 * test.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <common/common.h>

#include "efi.h"

int
main(int argc, char *argv[])
{
  efi_data *data;
  size_t i;

  fprintf(stderr, "main: reading all EFI files...");
  data = efi_read_idx("efi.idx", NULL);
  fprintf(stderr, "done (%zu data read)\n", data->n);

  fprintf(stderr, "main: sorting EFI data...");
  efi_sort(data);
  fprintf(stderr, "done\n");

  for (i = 0; i < data->n; ++i)
    {
      double lt = get_localtime(data->array[i].t, data->array[i].longitude);
      printf("%f %ld %f %f %f %f\n",
             lt,
             data->array[i].t,
             data->array[i].longitude,
             data->array[i].glat,
             data->array[i].galt,
             data->array[i].E_phi);
    }

  efi_free(data);

  return 0;
} /* main() */
