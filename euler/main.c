/*
 * main.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>

#include <satdata/satdata.h>

#include "common.h"
#include "oct.h"
#include "magdata.h"

#include "euler_calc.h"
#include "euler.h"

euler_workspace *
main_proc(const double period, magdata *data)
{
  int s = 0;
  const double dt = period * 86400000.0; /* period in ms */
  double t0 = data->t[0];
  double t1 = data->t[data->n - 1];
  const size_t nangles = (t1 - t0) / dt + 1;
#if 0
  size_t flags = EULER_FLG_ROTSC | EULER_FLG_ZYZ | EULER_FLG_RINV;
#else
  size_t flags = EULER_FLG_ROTSC | EULER_FLG_ZYX;
#endif
  double t;
  double x_data[3] = { 0.0, 0.0, 0.0 };
  gsl_vector_view xv = gsl_vector_view_array(x_data, 3);
  size_t i;
  euler_workspace *euler_p;

  fprintf(stderr, "main_proc: number of angle sets: %zu\n", nangles);

  euler_p = euler_alloc(flags);

  i = 0;
  for (t = t0; t <= t1; t += dt, i++)
    {
      /* clear previous discard flags */
      magdata_clear(data);

      /* flag data outside [t, t + dt] */
      magdata_flag_t(t, t + dt, data);

      /* initial guess */
      gsl_vector_set_zero(&xv.vector);

      s = euler_calc_proc(flags, &xv.vector, data);
      if (s)
        continue;

      /* store these angles in workspace */
      euler_add(t + 0.5 * dt, &xv.vector, euler_p);
    }

  /* cleanup */
  magdata_clear(data);

  return euler_p;
} /* main_proc() */

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options] sat.dat\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --output_file | -o file         - Euler angle output file\n");
  fprintf(stderr, "\t --euler | -p period             - Euler bin size in days\n");
  fprintf(stderr, "\t --residual_file | -r file       - residual output file\n");
} /* print_help() */

int
main(int argc, char *argv[])
{
  char *outfile = NULL;
  char *resfile = NULL;
  magdata *data;
  double euler_period = 3.0; /* set to 0 for single set of angles */
  struct timeval tv0, tv1;
  euler_workspace *euler_workspace_p;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "residual_file", required_argument, NULL, 'r' },
          { "output_file", required_argument, NULL, 'o' },
          { "maxit", required_argument, NULL, 'n' },
          { "euler", required_argument, NULL, 'p' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "o:p:r:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'o':
            outfile = optarg;
            break;

          case 'p':
            euler_period = atof(optarg);
            break;

          case 'r':
            resfile = optarg;
            break;

          default:
            print_help(argv);
            exit(1);
            break;
        }
    }

  if (argc - optind == 0)
    {
      print_help(argv);
      exit(1);
    }

  if (outfile)
    fprintf(stderr, "main: output coefficient file = %s\n", outfile);
  if (resfile)
    fprintf(stderr, "main: residual output file = %s\n", resfile);

  fprintf(stderr, "main: reading %s...", argv[optind]);
  gettimeofday(&tv0, NULL);
  data = magdata_read(argv[optind], NULL);
  gettimeofday(&tv1, NULL);

  if (!data)
    exit(1);

  fprintf(stderr, "done (%zu data total, %g seconds)\n",
          data->n, time_diff(tv0, tv1));

  magdata_init(data);
  magdata_calc(data);

  euler_workspace_p = main_proc(euler_period, data);

  if (resfile)
    {
      fprintf(stderr, "main: writing residuals to %s...", resfile);
      euler_calc_print_residuals(resfile, 0, data, euler_workspace_p);
      fprintf(stderr, "done\n");
    }

  if (outfile)
    {
      fprintf(stderr, "main: writing Euler angles to %s...", outfile);
      euler_write(outfile, euler_workspace_p);
      fprintf(stderr, "done\n");
    }

  magdata_free(data);

  return 0;
} /* main() */
