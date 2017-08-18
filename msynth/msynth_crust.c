/*
 * msynth_mf7.c
 * 
 * Routines to read a POMME .cof file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include "msynth.h"

static msynth_workspace *msynth_crust_read(const char *filename, const size_t nmax,
                                           const double epoch);

msynth_workspace *
msynth_mf7_read(const char *filename)
{
  const size_t nmax = 133;
  const double epoch = 2008.83;

  return msynth_crust_read(filename, nmax, epoch);
} /* msynth_mf7_read() */

msynth_workspace *
msynth_ngdc720_read(const char *filename)
{
  const size_t nmax = 740;
  const double epoch = 2010.0;

  return msynth_crust_read(filename, nmax, epoch);
} /* msynth_ngdc720_read() */

msynth_workspace *
msynth_emm_read(const char *filename)
{
  const size_t nmax = 790;
  const double epoch = 2017.0;

  return msynth_crust_read(filename, nmax, epoch);
} /* msynth_emm_read() */

static msynth_workspace *
msynth_crust_read(const char *filename, const size_t nmax, const double epoch)
{
  FILE *fp;
  msynth_workspace *w = NULL;
  char buffer[2048];

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "msynth_crust_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  w = msynth_alloc(nmax, 1, &epoch);
  msynth_set(16, nmax, w);

  while (fgets(buffer, sizeof(buffer), fp) != NULL)
    {
      int c;
      size_t cidx;
      size_t n;
      int m;
      double gnm, hnm;

      c = sscanf(buffer, "%zu %d %lf %lf", &n, &m, &gnm, &hnm);
      if (c != 4)
        continue;

      /* quality control (n,m) values */
      if (n > nmax)
        {
          fprintf(stderr, "msynth_crust_read: error: n = %zu\n", n);
          return w;
        }
      else if ((size_t) m > n)
        {
          fprintf(stderr, "msynth_crust_read: error: m(%d) > n(%zu)\n", m, n);
          return w;
        }

      cidx = msynth_nmidx(n, m, w);
      w->c[cidx] = gnm;
      w->c[cidx + w->sv_offset] = 0.0;
      w->c[cidx + w->sa_offset] = 0.0;

      if (m > 0)
        {
          cidx = msynth_nmidx(n, -m, w);
          w->c[cidx] = hnm;
          w->c[cidx + w->sv_offset] = 0.0;
          w->c[cidx + w->sa_offset] = 0.0;
        }
    }

  fclose(fp);

  return w;
} /* msynth_crust_read() */

int
msynth_crust_write(const char *filename, const msynth_workspace *w)
{
  size_t n;
  int m;
  FILE *fp;
  int iepoch = (int) w->epochs[0];

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "msynth_crust_write: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  /* print header information */
  fprintf(fp, "%% Crustal magnetic field model coefficients\n");
  fprintf(fp, "%% nmax:  %zu\n", w->nmax);
  fprintf(fp, "%% epoch: %.4f\n", w->epochs[0]);
  fprintf(fp, "%% radius: %.1f\n", w->R);
  fprintf(fp, "%% %3s %5s %10s %10s\n",
          "n",
          "m",
          "gnm (nT)",
          "hnm (nT)");

  for (n = 1; n <= w->nmax; ++n)
    {
      int ni = (int) n;

      for (m = 0; m <= ni; ++m)
        {
          double gnm, hnm;
          size_t cidx;

          cidx = msynth_nmidx(n, m, w);
          gnm = w->c[cidx];

          if (m == 0)
            {
              hnm = 0.0;
            }
          else
            {
              cidx = msynth_nmidx(n, -m, w);
              hnm = w->c[cidx];
            }

          fprintf(fp, "%5zu %5d %10.4f %10.4f\n",
                  n,
                  m,
                  gnm,
                  hnm);
        }
    }

  fclose(fp);

  return 0;
}
