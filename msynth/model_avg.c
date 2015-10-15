/*
 * model_avg.c
 *
 * Compute average of two models
 *
 * Usage: ./print_spectrum [args]
 * [-a coef_file]            - ASCII coefficient file
 * [-b coef_file]            - ASCII coefficient file
 * [-o output_file]
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
model_average(msynth_workspace *msynth1,
              const msynth_workspace *msynth2)
{
  const size_t nmax = msynth1->nmax;
  size_t n;

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;
      int m;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = msynth_nmidx(n, m, msynth1);
          double gnm1 = msynth1->c[cidx];
          double dgnm1 = msynth1->c[cidx + msynth1->sv_offset];
          double ddgnm1 = msynth1->c[cidx + msynth1->sa_offset];
          double gnm2 = msynth2->c[cidx];
          double dgnm2 = msynth2->c[cidx + msynth2->sv_offset];
          double ddgnm2 = msynth2->c[cidx + msynth2->sa_offset];

          msynth1->c[cidx] = 0.5*(gnm1 + gnm2);
          msynth1->c[cidx + msynth1->sv_offset] = 0.5*(dgnm1 + dgnm2);
          msynth1->c[cidx + msynth1->sa_offset] = 0.5*(ddgnm1 + ddgnm2);
        }
    }

  return 0;
} /* model_average() */

int
main(int argc, char *argv[])
{
  msynth_workspace *msynth1 = NULL;
  msynth_workspace *msynth2 = NULL;
  char *outfile = NULL;
  int c;

  while ((c = getopt(argc, argv, "a:b:o:")) != (-1))
    {
      switch (c)
        {
          case 'a':
            msynth1 = msynth_read(optarg);
            break;

          case 'b':
            msynth2 = msynth_read(optarg);
            break;

          case 'o':
            outfile = optarg;
            break;
        }
    }

  if (!msynth1 && !msynth2)
    {
      fprintf(stderr, "Usage: %s [-a coef_file] [-b coef_file] [-o output_file]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: computing model average...");
  model_average(msynth1, msynth2);
  fprintf(stderr, "done\n");

  if (outfile)
    {
      double epoch = msynth1->epochs[0];

      fprintf(stderr, "main: writing %s...", outfile);
      msynth_write(outfile, epoch, msynth1);
      fprintf(stderr, "done\n");
    }

  msynth_free(msynth1);
  msynth_free(msynth2);

  return 0;
} /* main() */
