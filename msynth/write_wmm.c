/*
 * write_wmm.c
 *
 * Write official formatted WMM coefficient file from
 * existing coefficient file
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
  msynth_workspace *w_sv = NULL;
  char *outfile = "wmm_output.cof";
  int c;

  while ((c = getopt(argc, argv, "c:w:o:a:v:")) != (-1))
    {
      switch (c)
        {
          case 'c':
            fprintf(stderr, "main: reading coefficient file %s...", optarg);
            w = msynth_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'a':
            fprintf(stderr, "main: reading Arnaud coefficient file %s...", optarg);
            w = msynth_arnaud_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'w':
            fprintf(stderr, "main: reading WMM coefficient file %s...", optarg);
            w = msynth_wmm_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'v':
            fprintf(stderr, "main: reading WMM coefficient file %s...", optarg);
            w_sv = msynth_wmm_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'o':
            outfile = optarg;
            break;
        }
    }

  if (!w)
    {
      fprintf(stderr, "Usage: %s [-c coef_file] [-w wmm_file] [-v wmm_file_sv] [-a arnaud_file] [-o output_file]\n", argv[0]);
      exit(1);
    }

  if (w_sv)
    {
      fprintf(stderr, "main: replacing SV coefficients...");
      msynth_wmm_replace_sv(w_sv, w);
      fprintf(stderr, "done\n");
    }

  fprintf(stderr, "main: writing WMM file %s...", outfile);
  msynth_wmm_write(outfile, w);
  fprintf(stderr, "done\n");

  msynth_free(w);

  return 0;
} /* main() */
