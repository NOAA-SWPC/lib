/*
 * msynth_emm.c
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

/*
msynth_emm_write()
  Write EMM main field coefficient file

Inputs: filename - output file
        year     - model epoch
        w        - workspace
*/

int
msynth_emm_write(const char *filename, const double year,
                 const msynth_workspace *w)
{
  const size_t epoch_idx = msynth_epoch_idx(year, w);
  const double epoch = w->epochs[epoch_idx];
  const int iepoch = (int) epoch;
  const double *g = w->c + epoch_idx * w->p;
  const size_t nmax = 15;
  size_t n;
  int m;
  FILE *fp;
  time_t t;
  struct tm *tm_p;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "msynth_emm_write: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  t = time(0);
  tm_p = gmtime(&t);

  fprintf(fp, "%10.1f %13s EMM%d %12s %02d/%02d/%4d\n",
          epoch,
          "",
          iepoch,
          "",
#if 0
          tm_p->tm_mon + 1,
          tm_p->tm_mday,
#else
          1,
          1,
#endif
          tm_p->tm_year + 1900);

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;

      for (m = 0; m <= ni; ++m)
        {
          double gnm, hnm;
          size_t cidx;

          cidx = msynth_nmidx(n, m, w);
          gnm = g[cidx];

          if (m == 0)
            {
              hnm = 0.0;
            }
          else
            {
              cidx = msynth_nmidx(n, -m, w);
              hnm = g[cidx];
            }

          fprintf(fp, "%3zu %3d %12.4f %12.4f\n",
                  n,
                  m,
                  gnm,
                  hnm);
        }
    }

  fclose(fp);

  return 0;
} /* msynth_emm_write() */

/*
msynth_emm_write_sv()
  Write EMM SV coefficient file

Inputs: filename - output file
        year     - model epoch
        w        - workspace
*/

int
msynth_emm_write_sv(const char *filename, const double year,
                    const msynth_workspace *w)
{
  const size_t epoch_idx = msynth_epoch_idx(year, w);
  const double *g = w->c + epoch_idx * w->p;
  const size_t nmax = 15;
  size_t n;
  int m;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "msynth_emm_write_sv: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;

      for (m = 0; m <= ni; ++m)
        {
          double dgnm, dhnm;
          size_t cidx;

          cidx = msynth_nmidx(n, m, w);
          dgnm = g[cidx + w->sv_offset];

          if (m == 0)
            {
              dhnm = 0.0;
            }
          else
            {
              cidx = msynth_nmidx(n, -m, w);
              dhnm = g[cidx + w->sv_offset];
            }

          fprintf(fp, "%3zu %2d %12.4f %12.4f\n",
                  n,
                  m,
                  dgnm,
                  dhnm);
        }
    }

  fclose(fp);

  return 0;
} /* msynth_emm_write_sv() */
