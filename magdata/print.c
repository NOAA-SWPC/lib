/*
 * print.c
 *
 * Print contents of magdata file
 *
 * ./print <-i magdata_file>
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <assert.h>
#include <errno.h>

#include "common.h"
#include "euler.h"
#include "magdata.h"

/*
magdata_euler_apply()
  Apply Euler angle rotations to magdata dataset

Inputs: data  - magdata
        w     - workspace

Return: success or error

Notes:
1) data->B_nec is replaced with data->B_vfm vectors rotated
with Euler angles and satellite quaternions
*/

int
magdata_euler_apply(magdata *data, const euler_workspace *w)
{
  int s = 0;
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      double *q = &(data->q[4 * i]);
      double B_vfm[3], B_nec[3];
      double t = data->t[i];

      B_vfm[0] = data->Bx_vfm[i];
      B_vfm[1] = data->By_vfm[i];
      B_vfm[2] = data->Bz_vfm[i];

      s += euler_vfm2nec_t(t, q, B_vfm, B_nec, w);

      data->Bx_nec[i] = B_nec[0];
      data->By_nec[i] = B_nec[1];
      data->Bz_nec[i] = B_nec[2];
    }

  return s;
} /* magdata_euler_apply() */

int
main(int argc, char *argv[])
{
  magdata *data = NULL;
  euler_workspace *euler_p = NULL;
  char *outfile = "data.dat";
  struct timeval tv0, tv1;
  int c;

  while ((c = getopt(argc, argv, "i:o:e:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);

            data = magdata_read(optarg, NULL);
            if (!data)
              {
                fprintf(stderr, "main: error reading %s\n", optarg);
                exit(1);
              }

            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
                    time_diff(tv0, tv1));
            break;

          case 'o':
            outfile = optarg;
            break;

          case 'e':
            fprintf(stderr, "main: reading Euler file %s...", optarg);
            euler_p = euler_read(optarg);
            fprintf(stderr, "done\n");

          default:
            break;
        }
    }

  if (!data)
    {
      fprintf(stderr, "Usage: %s <-i magdata_file> [-o output_file] [-e euler_file]\n",
              argv[0]);
      exit(1);
    }

  if (euler_p)
    {
      fprintf(stderr, "main: applying Euler angles to data...");
      magdata_euler_apply(data, euler_p);
      fprintf(stderr, "done\n");
    }

  fprintf(stderr, "main: writing output to %s...", outfile);
  magdata_print(outfile, data);
  fprintf(stderr, "done\n");

  magdata_free(data);

  return 0;
} /* main() */
