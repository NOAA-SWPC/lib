/*
 * extrapolate_epoch.c
 *
 * Read in coefficient file, extrapolate to new epoch using SV/SA,
 * and write new output file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include "common.h"
#include "msynth.h"

int
main(int argc, char *argv[])
{
  msynth_workspace *w = NULL;
  int c;
  size_t nmax = 0;
  double epoch = -1.0;
  char *outfile = NULL;
  double data_start = -1.0;
  double data_end = -1.0;

  while ((c = getopt(argc, argv, "c:m:w:o:a:h:p:n:i:e:y:z:")) != (-1))
    {
      switch (c)
        {
          case 'e':
            epoch = atof(optarg);
            break;

          case 'c':
            w = msynth_read(optarg);
            break;

          case 'a':
            w = msynth_arnaud_read(optarg);
            break;

          case 'm':
            w = msynth_mf7_read(optarg);
            break;

          case 'w':
            w = msynth_wmm_read(optarg);
            break;

          case 'h':
            w = msynth_chaos_read(optarg);
            break;

          case 'p':
            w = msynth_pomme_read(optarg);
            break;

          case 'i':
            w = msynth_igrf_read_mf(optarg);
            msynth_igrf_read_sv("cof/IGRF12/SV-2015-2020-E.cof", w);
            break;

          case 'n':
            nmax = (size_t) atoi(optarg);
            break;

          case 'o':
            outfile = optarg;
            break;

          case 'y':
            data_start = atof(optarg);
            break;

          case 'z':
            data_end = atof(optarg);
            break;
        }
    }

  if (!w)
    {
      fprintf(stderr, "Usage: %s [-c coef_file] [-m mf7_file] [-w wmm_file] [-a arnaud_file] [-h chaos_file] [-p pomme_file] [-i igrf12_mf_candidate] [-n nmax] [-e epoch] [-y data_start] [-z data_end] [-o output_file]\n", argv[0]);
      exit(1);
    }

  if (epoch < 0.0)
    {
      fprintf(stderr, "main: error: no epoch specified\n");
      exit(1);
    }

  msynth_set_data_interval(data_start, data_end, w);

  fprintf(stderr, "main: epoch = %g\n", epoch);
  if (data_start > 0.0)
    fprintf(stderr, "main: data_start = %g\n", data_start);
  if (data_end > 0.0)
    fprintf(stderr, "main: data_end = %g\n", data_end);

  fprintf(stderr, "main: extrapolating model coefficients to epoch %g...", epoch);
  msynth_extrapolate_g(epoch, w);
  fprintf(stderr, "done\n");

  if (outfile)
    {
      fprintf(stderr, "main: writing coefficients to %s...", outfile);
      msynth_write(outfile, epoch, w);
      fprintf(stderr, "done\n");
    }

  msynth_free(w);

  return 0;
} /* main() */
