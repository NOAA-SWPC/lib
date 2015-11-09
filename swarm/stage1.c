/*
 * stage1.c
 * Patrick Alken
 *
 * 1. Read Swarm CDF data file
 * 2. Compute main and external fields
 * 3. Create new CDF with main/external field outputs
 *
 * Usage:
 *
 * ./stage1 <-i swarm_cdf_file>
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>

#include <satdata/satdata.h>

#include "chaos.h"
#include "common.h"
#include "msynth.h"
#include "track.h"

int
main(int argc, char *argv[])
{
  satdata_mag *data_in = NULL;
  satdata_mag *data_out;
  satdata_lp *lp_data = NULL;
  struct timeval tv0, tv1;
  int c;
  char *outfile = NULL;
  msynth_workspace *msynth_core;
  int down_sample = 1;
  int use_chaos = 0;
  int status;

  while ((c = getopt(argc, argv, "a:ci:o:d:l:")) != (-1))
    {
      switch (c)
        {
          case 'c':
            use_chaos = 1;
            break;

          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);

            data_in = satdata_swarm_read(optarg, NULL);
            if (!data_in)
              {
                fprintf(stderr, "main: error reading %s\n", optarg);
                exit(1);
              }

            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data_in->n,
                    time_diff(tv0, tv1));
            break;

          case 'a':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);

            data_in = satdata_swarm_asmv_read(optarg, NULL);
            if (!data_in)
              {
                fprintf(stderr, "main: error reading %s\n", optarg);
                exit(1);
              }

            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data_in->n,
                    time_diff(tv0, tv1));
            break;

          case 'o':
            outfile = optarg;
            break;

          case 'd':
            down_sample = atoi(optarg);
            break;

          case 'l':
            fprintf(stderr, "main: reading Langmuir data from %s...", optarg);
            gettimeofday(&tv0, NULL);

            lp_data = satdata_swarm_lp_read(optarg, 0);

            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu data read, %g seconds)\n", lp_data->n,
                    time_diff(tv0, tv1));

            break;

          default:
            break;
        }
    }

  if (!data_in || !outfile)
    {
      fprintf(stderr, "Usage: %s <-i swarm_cdf_file> <-a swarm_asmv_cdf_file> <-o output_cdf_file> [-d down_sample] [-c]\n",
              argv[0]);
      exit(1);
    }

  {
    size_t nF;

    fprintf(stderr, "main: replacing missing scalar data from vector...");
    satdata_swarm_filter_instrument(0, data_in);
    nF = satdata_mag_fill_F(data_in);
    fprintf(stderr, "done (%zu missing scalar data replaced)\n", nF);
  }

  if (lp_data)
    {
      fprintf(stderr, "main: filling Ne data...");
      satdata_mag_fill_ne(data_in, lp_data);
      fprintf(stderr, "done\n");
    }

  if (use_chaos)
    chaos_init();

  fprintf(stderr, "output file = %s\n", outfile);

  if (data_in->n == 0)
    {
      fprintf(stderr, "main: no data to process\n");
      exit(1);
    }

  data_out = satdata_mag_alloc(data_in->n);
  msynth_core = msynth_read(MSYNTH_BOUMME_FILE);

  fprintf(stderr, "Computing magnetic field model (downsampling by %d)...", down_sample);
  gettimeofday(&tv0, NULL);

  if (use_chaos)
    status = track_synth_chaos(down_sample, data_in, data_out, msynth_core);
  else
    status = track_synth(down_sample, data_in, data_out, msynth_core);

  if (status)
    {
      /* don't produce output file if error */
      fprintf(stderr, "main: error computing residuals\n");
      exit(1);
    }
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, %zu records kept)\n", time_diff(tv0, tv1), data_out->n);

  fprintf(stderr, "Writing %s...", outfile);
  gettimeofday(&tv0, NULL);
  satdata_swarm_write(outfile, data_out);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records written, %g seconds)\n", data_out->n,
          time_diff(tv0, tv1));

  if (use_chaos)
    chaos_term();

  satdata_mag_free(data_in);
  satdata_mag_free(data_out);
  msynth_free(msynth_core);

  return 0;
} /* main() */
