/*
 * msynth_wmm.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include "msynth.h"

msynth_workspace *
msynth_wmm_read(const char *filename)
{
  FILE *fp;
  const size_t nmax = 12;
  int gotepoch = 0;
  msynth_workspace *w = NULL;
  char buffer[2048], scratch[2048];

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "msynth_wmm_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  while (fgets(buffer, 2048, fp) != NULL)
    {
      int n, m, c;
      double gnm, hnm, dgnm, dhnm;
      size_t cidx;

      if (!gotepoch)
        {
          double epoch;
          int mon, day, year;

          c = sscanf(buffer, "%lf %s %d/%d/%d",
                     &epoch, scratch, &mon, &day, &year);
          if (c < 5)
            continue;

          w = msynth_alloc(nmax, 1, &epoch);
          gotepoch = 1;
        }
      else
        {
          assert(w != NULL);

          c = sscanf(buffer, "%d %d %lf %lf %lf %lf",
                     &n, &m, &gnm, &hnm, &dgnm, &dhnm);
          if (c < 6)
            continue;

          assert(n <= 12);
          assert(m <= n);

          cidx = msynth_nmidx(n, m, w);
          w->c[cidx] = gnm;
          w->c[cidx + w->sv_offset] = dgnm;
          w->c[cidx + w->sa_offset] = 0.0;

          if (m != 0)
            {
              cidx = msynth_nmidx(n, -m, w);
              w->c[cidx] = hnm;
              w->c[cidx + w->sv_offset] = dhnm;
              w->c[cidx + w->sa_offset] = 0.0;
            }
        }
    }

  fclose(fp);

  return w;
} /* msynth_wmm_read() */

int
msynth_wmm_write(const char *filename, const msynth_workspace *w)
{
  const size_t nmax = 12;
  const size_t nmax_sv = 12;
  size_t n;
  int m;
  FILE *fp;
  int iepoch = (int) w->epochs[0];
  time_t t;
  struct tm *tm_p;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "msynth_wmm_write: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  t = time(0);
  tm_p = gmtime(&t);

  fprintf(fp, "%10.1f %10s WMM-%d %6s %02d/%02d/%4d\n",
          w->epochs[0],
          "",
          iepoch,
          "",
          tm_p->tm_mon + 1,
          tm_p->tm_mday,
          tm_p->tm_year + 1900);

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;

      for (m = 0; m <= ni; ++m)
        {
          double gnm, hnm, dgnm, dhnm;
          size_t cidx;

          cidx = msynth_nmidx(n, m, w);
          gnm = w->c[cidx];
          dgnm = w->c[cidx + w->sv_offset];

          if (m == 0)
            {
              hnm = 0.0;
              dhnm = 0.0;
            }
          else
            {
              cidx = msynth_nmidx(n, -m, w);
              hnm = w->c[cidx];
              dhnm = w->c[cidx + w->sv_offset];
            }

          if (n > nmax_sv)
            dgnm = dhnm = 0.0;

          fprintf(fp, "%3zu %2d %9.1f %9.1f %10.1f %10.1f\n",
                  n,
                  m,
                  gnm,
                  hnm,
                  dgnm,
                  dhnm);
        }
    }

  /* needed for some military customers */
  fprintf(fp, "999999999999999999999999999999999999999999999999\n");
  fprintf(fp, "999999999999999999999999999999999999999999999999\n");

  fclose(fp);

  return 0;
} /* msynth_wmm_write() */

int
msynth_wmm_replace_sv(const msynth_workspace *w_sv, msynth_workspace *w)
{
  const size_t nmax = 12;
  size_t n;

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;
      int m;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = msynth_nmidx(n, m, w);
          w->c[cidx + w->sv_offset] = w_sv->c[cidx + w_sv->sv_offset];
        }
    }

  return 0;
} /* msynth_wmm_replace_sv() */
