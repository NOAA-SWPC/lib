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
