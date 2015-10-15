/*
 * msynth_pomme.c
 * 
 * Routines to read a POMME .cof file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include "msynth.h"

msynth_workspace *
msynth_pomme_read(const char *filename)
{
  FILE *fp;
  const size_t nmax = 130;
  msynth_workspace *w = NULL;
  char buffer[2048];
  int got_epoch = 0;
  int got_nm = 0;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "msynth_pomme_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  while (fgets(buffer, sizeof(buffer), fp) != NULL)
    {
      int c;

      if (!got_epoch)
        {
          double epoch;

          c = sscanf(buffer, "%lf", &epoch);
          if (c == 1)
            {
              got_epoch = 1;
              w = msynth_alloc(nmax, 1, &epoch);
            }
        }
      else if (!got_nm)
        {
          char n, m;

          c = sscanf(buffer, "%*[ \t\n]%c %c", &n, &m);
          if (c == 2 && n == 'n' && m == 'm')
            got_nm = 1; /* ready to read in coefficients */
        }
      else
        {
          size_t cidx;
          size_t n;
          int m;
          double gnm, hnm;
          double dgnm = 0.0, dhnm = 0.0;
          double ddgnm = 0.0, ddhnm = 0.0;

          c = sscanf(buffer, "%zu %d %lf %lf %lf %lf %lf %lf",
                     &n, &m, &gnm, &hnm, &dgnm, &dhnm,
                     &ddgnm, &ddhnm);
          if (c != 8 && c != 4)
            continue;

          /* quality control (n,m) values */
          if (n > nmax)
            {
              fprintf(stderr, "msynth_pomme_read: error: n = %zu\n", n);
              return w;
            }
          else if ((size_t) m > n)
            {
              fprintf(stderr, "msynth_pomme_read: error: m(%d) > n(%zu)\n", m, n);
              return w;
            }

          cidx = msynth_nmidx(n, m, w);
          w->c[cidx] = gnm;
          w->c[cidx + w->sv_offset] = dgnm;
          w->c[cidx + w->sa_offset] = ddgnm;

          if (m > 0)
            {
              cidx = msynth_nmidx(n, -m, w);
              w->c[cidx] = hnm;
              w->c[cidx + w->sv_offset] = dhnm;
              w->c[cidx + w->sa_offset] = ddhnm;
            }
        }
    }

  fclose(fp);

  if (!got_epoch || !got_nm)
    {
      fprintf(stderr, "msynth_pomme_read: badly formatted file\n");

      if (w)
        msynth_free(w);

      return 0;
    }

  return w;
} /* msynth_pomme_read() */
