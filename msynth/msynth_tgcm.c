/*
 * msynth_tgcm.c
 *
 * Reading routine for TIEGCM equivalent current SH file
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
msynth_tgcm_read(const char *filename)
{
  FILE *fp;
  const size_t nmax = 180;
  msynth_workspace *w;
  char buffer[2048];

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "msynth_tgcm_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  w = msynth_alloc(nmax, 1, NULL);

  while (fgets(buffer, 2048, fp) != NULL)
    {
      int n, m, c;
      double gnm, hnm;
      size_t cidx;

      c = sscanf(buffer, "%d %d %lf %lf",
                 &m, &n, &gnm, &hnm);
      if (c < 4)
        continue;

      assert(m <= n);

      gnm *= sqrt(1.0 / (n + 0.5));
      hnm *= sqrt(1.0 / (n + 0.5));

      cidx = msynth_nmidx(n, m, w);
      w->c[cidx] = gnm * 1.0e9; /* convert to nT */
      w->c[cidx + w->sv_offset] = 0.0;
      w->c[cidx + w->sa_offset] = 0.0;

      if (m != 0)
        {
          cidx = msynth_nmidx(n, -m, w);
          w->c[cidx] = hnm * 1.0e9; /* convert to nT */
          w->c[cidx + w->sv_offset] = 0.0;
          w->c[cidx + w->sa_offset] = 0.0;
        }
    }

  fclose(fp);

  return w;
}
