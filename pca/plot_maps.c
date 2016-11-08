/*
 * plot_maps.c
 *
 * Print lat/lon maps of magnetic field and current flow due to
 * first 3 PCs in time domain
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>

#include "common.h"
#include "pca.h"

int
main(int argc, char *argv[])
{
  pca_workspace *pca_workspace_p;
  char buf[2048];
  size_t i;
  double r = R_EARTH_KM + 0.0; /* radius for magnetic field maps */

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            r = R_EARTH_KM + atof(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s [-a altitude]\n", argv[0]);
            break;
        }
    }

  fprintf(stderr, "main: radius = %g [km]\n", r);

  pca_workspace_p = pca_alloc();
  
  for (i = 0; i < 40; ++i)
    {
      sprintf(buf, "maps/map%02zu.dat", i);

      fprintf(stderr, "main: printing principle component %zu map to %s...",
              i, buf);

      pca_print_pc_map(buf, r, i, pca_workspace_p);

      fprintf(stderr, "done\n");
    }

  pca_free(pca_workspace_p);

  return 0;
}
