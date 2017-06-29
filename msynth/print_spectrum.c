/*
 * print_spectrum.c
 *
 * Print model spectrum
 *
 * Usage: ./print_spectrum [args]
 * [-c coef_file]            - ASCII coefficient file
 * [-m mf7_file]
 * [-w wmm_file]
 * [-a arnaud_file]
 * [-b NGDC720_file]
 * [-f EMM_crust_file]
 * [-h swarm_file]
 * [-p pomme_file]
 * [-i igrf12_mf_candidate]
 * [-g ipgp_file]
 * [-z]
 * [-n nmax]
 * [-e epoch]
 * [-o output_file]
 * [-d]                      - Compute model difference with CHAOS
 *
 * Note: if two coefficient files are specified, the model difference is computed
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
  msynth_workspace *msynth1 = NULL;
  msynth_workspace *msynth2 = NULL;
  char *outfile = "spectrum.s";
  char *corrfile = "corr.dat";
  int print_azim = 0;
  int c;
  size_t nmax = 0;
  double epoch = -1.0;

  while ((c = getopt(argc, argv, "a:b:c:e:f:g:m:w:o:h:p:zn:i:dt:")) != (-1))
    {
      switch (c)
        {
          case 'c':
            if (!msynth1)
              msynth1 = msynth_read(optarg);
            else
              msynth2 = msynth_read(optarg);
            break;

          case 'a':
            if (!msynth1)
              msynth1 = msynth_arnaud_read(optarg);
            else
              msynth2 = msynth_arnaud_read(optarg);
            break;

          case 'b':
            if (!msynth1)
              msynth1 = msynth_ngdc720_read(optarg);
            else
              msynth2 = msynth_ngdc720_read(optarg);
            break;

          case 'f':
            if (!msynth1)
              msynth1 = msynth_emm_read(optarg);
            else
              msynth2 = msynth_emm_read(optarg);
            break;

          case 'm':
            if (!msynth1)
              msynth1 = msynth_mf7_read(optarg);
            else
              msynth2 = msynth_mf7_read(optarg);
            break;

          case 'w':
            if (!msynth1)
              msynth1 = msynth_wmm_read(optarg);
            else
              msynth2 = msynth_wmm_read(optarg);
            break;

          case 'h':
            if (!msynth1)
              msynth1 = msynth_swarm_read(optarg);
            else
              msynth2 = msynth_swarm_read(optarg);
            break;

          case 'p':
            if (!msynth1)
              msynth1 = msynth_pomme_read(optarg);
            else
              msynth2 = msynth_pomme_read(optarg);
            break;

          case 'g':
            if (!msynth1)
              msynth1 = msynth_ipgp_read(optarg);
            else
              msynth2 = msynth_ipgp_read(optarg);
            break;

          case 'i':
            if (!msynth1)
              {
                msynth1 = msynth_igrf_read_mf(optarg);
                msynth_igrf_read_sv("cof/IGRF12/SV-2015-2020-E.cof", msynth1);
              }
            else
              {
                msynth2 = msynth_igrf_read_mf(optarg);
                msynth_igrf_read_sv("cof/IGRF12/SV-2015-2020-E.cof", msynth2);
              }
            break;

          case 't':
            if (!msynth1)
              msynth1 = msynth_tgcm_read(optarg);
            else
              msynth2 = msynth_tgcm_read(optarg);
            break;

          case 'z':
            print_azim = 1;
            break;

          case 'n':
            nmax = (size_t) atoi(optarg);
            break;

          case 'e':
            epoch = atof(optarg);
            break;

          case 'o':
            outfile = optarg;
            break;

          case 'd':
            msynth2 = msynth_swarm_read(MSYNTH_CHAOS_FILE);
            break;
        }
    }

  if (!msynth1)
    {
      fprintf(stderr, "Usage: %s [-c coef_file] [-m mf7_file] [-w wmm_file] [-a arnaud_file] [-b NGDC720_file] [-f EMM_crust_file] [-h swarm_shc_file] [-p pomme_file] [-t tgcm_file] [-i igrf12_mf_candidate] [-g ipgp_file] [-z] [-n nmax] [-e epoch] [-o output_file] [-d]\n", argv[0]);
      exit(1);
    }

  if (epoch < 0.0)
    epoch = msynth1->epochs[0];

  if (msynth2)
    {
      /* compute model difference */
      fprintf(stderr, "main: computing coefficient difference...");
      w = msynth_diff(epoch, msynth1, msynth2);
      fprintf(stderr, "done\n");
    }
  else
    w = msynth_copy(msynth1);

  if (nmax != 0)
    msynth_set(1, nmax, w);

  fprintf(stderr, "main: epoch = %g\n", epoch);
  fprintf(stderr, "main: spectrum nmax = %zu\n", w->eval_nmax);

  fprintf(stderr, "main: printing spectrum to %s...", outfile);
  if (print_azim)
    msynth_print_spectrum_m(outfile, w);
  else
    msynth_print_spectrum(outfile, epoch, w);
  fprintf(stderr, "done\n");

  if (msynth2)
    {
      fprintf(stderr, "main: printing degree correlation to %s...", corrfile);
      msynth_print_correlation(corrfile, msynth1, msynth2);
      fprintf(stderr, "done\n");
    }

#if 0
  /*XXX*/
  {
    const char *sfile = "new.txt";
    msynth_swarm_write(sfile, msynth1);
  }
#endif

  msynth_free(w);
  msynth_free(msynth1);

  if (msynth2)
    msynth_free(msynth2);

  return 0;
} /* main() */
