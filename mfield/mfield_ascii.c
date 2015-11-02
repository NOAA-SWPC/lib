/*
 * mfield_ascii.c
 *
 * Rewrite binary coefficient in ascii format
 *
 * usage: mfield_ascii -c dmsp_coef_file -o ascii_coef_file -e new_epoch
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include "mfield.h"

int
write_coef_block(const char *filename, const double epoch, mfield_workspace *w)
{
  int s = 0;
  FILE *fp;
  size_t n;
  gsl_vector *c = w->c_copy;
  int iout;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "write_coef_block: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  /* convert coefficients to physical time units */
  mfield_coeffs(1, w->c, c, w);

  fprintf(fp, "{%.1f, %.1f, %.1f, %.1f, \"HDGM2014_Core_%g\",\n",
          2014.0,
          epoch,
          epoch,
          epoch + 1.0,
          epoch);

  /* g_nm */
  fprintf(fp, "{0,");
  iout = 0;
  for (n = 1; n <= w->nmax_mf; ++n)
    {
      int m, ni = (int) n;

      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);
          double gnm = gsl_vector_get(c, cidx); 

          fprintf(fp, "%6.4f", gnm);
          ++iout;
          if (n < w->nmax_mf || m < ni)
            {
              putc(',', fp);
              if ((iout & 7) == 0)
                putc('\n', fp);
            }
          else
            fprintf(fp, "},\n");
        }
    }

  /* h_nm */
  fprintf(fp, "{0,");
  iout = 0;
  for (n = 1; n <= w->nmax_mf; ++n)
    {
      int m, ni = (int) n;

      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, -m);
          double gnm = gsl_vector_get(c, cidx); 

          if (m == 0)
            fprintf(fp, "%6.4f", 0.0);
          else
            fprintf(fp, "%6.4f", gnm);

          ++iout;
          if (n < w->nmax_mf || m < ni)
            {
              putc(',', fp);
              if ((iout & 7) == 0)
                putc('\n', fp);
            }
          else
            fprintf(fp, "},\n");
        }
    }

  /* dg_nm */
  fprintf(fp, "{0,");
  iout = 0;
  for (n = 1; n <= w->nmax_mf; ++n)
    {
      int m, ni = (int) n;

      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);
          double dgnm = gsl_vector_get(c, cidx + w->sv_offset); 

          fprintf(fp, "%6.4f", dgnm);
          ++iout;
          if (n < w->nmax_mf || m < ni)
            {
              putc(',', fp);
              if ((iout & 7) == 0)
                putc('\n', fp);
            }
          else
            fprintf(fp, "},\n");
        }
    }

  /* dh_nm */
  fprintf(fp, "{0,");
  iout = 0;
  for (n = 1; n <= w->nmax_mf; ++n)
    {
      int m, ni = (int) n;

      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, -m);
          double dgnm = gsl_vector_get(c, cidx + w->sv_offset); 

          if (m == 0)
            fprintf(fp, "%6.4f", 0.0);
          else
            fprintf(fp, "%6.4f", dgnm);

          ++iout;
          if (n < w->nmax_mf || m < ni)
            {
              putc(',', fp);
              if ((iout & 7) == 0)
                putc('\n', fp);
            }
          else
            fprintf(fp, "},\n");
        }
    }

  fprintf(fp, "%zu,%zu,1},", w->nmax_mf, w->nmax_mf);

  fclose(fp);

  return s;
} /* write_coef_block() */

/* write Stefan's my_ files */
int
write_coef_block_stefan(const char *filename, const double epoch, mfield_workspace *w)
{
  int s = 0;
  FILE *fp;
  size_t n;
  gsl_vector *c = w->c_copy;
  int iout;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "write_coef_block_stefan: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  /* convert coefficients to physical time units */
  mfield_coeffs(1, w->c, c, w);

  /* g_nm */
  fprintf(fp, "{{");
  iout = 0;
  for (n = 1; n <= w->nmax_mf; ++n)
    {
      int m, ni = (int) n;

      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);
          double gnm = gsl_vector_get(c, cidx); 

          fprintf(fp, "%6.4f", gnm);
          ++iout;
          putc(',', fp);
          if ((iout & 7) == 0)
            putc('\n', fp);

          if (m > 0)
            {
              cidx = mfield_coeff_nmidx(n, -m);
              gnm = gsl_vector_get(c, cidx); 

              fprintf(fp, "%6.4f", gnm);
              ++iout;
              if (n < w->nmax_mf || m < ni)
                {
                  putc(',', fp);
                  if ((iout & 7) == 0)
                    putc('\n', fp);
                }
              else
                fprintf(fp, "},\n");
            }
        }
    }

  /* dg_nm */
  fprintf(fp, "{");
  iout = 0;
  for (n = 1; n <= w->nmax_mf; ++n)
    {
      int m, ni = (int) n;

      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);
          double dgnm = gsl_vector_get(c, cidx + w->sv_offset); 

          fprintf(fp, "%6.4f", dgnm);
          ++iout;
          putc(',', fp);
          if ((iout & 7) == 0)
            putc('\n', fp);

          if (m > 0)
            {
              cidx = mfield_coeff_nmidx(n, -m);
              dgnm = gsl_vector_get(c, cidx + w->sv_offset); 

              fprintf(fp, "%6.4f", dgnm);
              ++iout;
              if (n < w->nmax_mf || m < ni)
                {
                  putc(',', fp);
                  if ((iout & 7) == 0)
                    putc('\n', fp);
                }
              else
                fprintf(fp, "}},\n");
            }
        }
    }

  fclose(fp);

  return s;
} /* write_coef_block_stefan() */

int
print_spectrum(const char *filename, mfield_workspace *w)
{
  size_t n;
  FILE *fp = fopen(filename, "w");

  fprintf(stderr, "print_spectrum: writing spectrum to %s...", filename);
  for (n = 1; n <= w->nmax_mf; ++n)
    {
      double gn = mfield_spectrum(n, w);
      double dgn = mfield_spectrum_sv(n, w);
      double ddgn = mfield_spectrum_sa(n, w);
      fprintf(fp, "%zu %.12e %.12e %.12e\n", n, gn, dgn, ddgn);
    }
  fprintf(stderr, "done\n");

  fclose(fp);

  return 0;
} /* print_spectrum() */

int
main(int argc, char *argv[])
{
  int c;
  char *coeffile = NULL;
  char *outfile = NULL;
  char *blockfile = NULL;
  char *specfile = NULL;
  double new_epoch = -1.0;
  mfield_workspace *mfield_workspace_p;
  int write_delta = 0;

  while ((c = getopt(argc, argv, "c:o:e:b:s:u")) != (-1))
    {
      switch (c)
        {
          case 'c':
            coeffile = optarg;
            break;

          case 'o':
            outfile = optarg;
            break;

          case 'e':
            new_epoch = atof(optarg);
            break;

          case 'b':
            blockfile = optarg;
            break;

          case 's':
            specfile = optarg;
            break;

          case 'u':
            write_delta = 1;

          default:
            break;
        }
    }

  if (!coeffile || !outfile)
    {
      fprintf(stderr, "usage: %s [-c coef_file] [-o ascii_coef_file] [-e new_epoch] [-b coef_block_file] [-s spectrum_file] [-u]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: loading coefficients from %s...", coeffile);
  mfield_workspace_p = mfield_read(coeffile);
  fprintf(stderr, "done\n");

  if (new_epoch > 0.0)
    {
      fprintf(stderr, "main: extrapolating coefficients from %g to %g\n",
              mfield_workspace_p->epoch, new_epoch);
      mfield_new_epoch(new_epoch, mfield_workspace_p);
      fprintf(stderr, "done\n");
    }
  else
    new_epoch = mfield_workspace_p->epoch;

  fprintf(stderr, "main: writing coefficients to %s...", outfile);
  mfield_write_ascii(outfile, new_epoch, write_delta, mfield_workspace_p);
  fprintf(stderr, "done\n");

  if (blockfile)
    {
      fprintf(stderr, "main: writing coefficient block to %s...", blockfile);
#if 1
      write_coef_block(blockfile, new_epoch, mfield_workspace_p);
#else
      write_coef_block_stefan(blockfile, new_epoch, mfield_workspace_p);
#endif
      fprintf(stderr, "done\n");
    }

  if (specfile)
    print_spectrum(specfile, mfield_workspace_p);

  mfield_free(mfield_workspace_p);

  return 0;
}
