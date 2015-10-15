/*
 * msynth_arnaud.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include "msynth.h"

msynth_workspace *
msynth_arnaud_read(const char *filename)
{
  FILE *fp;
  const size_t nmax = 20;
  int gotepoch = 0;
  msynth_workspace *w = NULL;
  char buffer[2048], scratch[2048];

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "msynth_arnaud_read: unable to open %s: %s\n",
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

          c = sscanf(buffer, "%lf %s",
                     &epoch, scratch);
          if (c < 2)
            continue;

          if (strcasecmp(scratch, "epoch") != 0)
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

          assert(n <= nmax);
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
} /* msynth_arnaud_read() */
