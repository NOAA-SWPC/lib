/*
 * msynth_bggm.c
 * 
 * Routines to read a BGGM formatted coefficient file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include "msynth.h"

/*
msynth_bggm_read()
  Read BGGM coefficient file
*/

msynth_workspace *
msynth_bggm_read(const char *filename)
{
  FILE *fp;
  const size_t nmax = 133;
  const size_t nsnapshot = 40;
  size_t n = 1;
  int m = 0;
  msynth_workspace *w = NULL;
  char buffer[8196];
  size_t nepoch = 0;  /* number of model epochs found so far */
  int read_block = 0; /* currently reading coefficient block? */
  size_t i;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "msynth_bggm_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  w = msynth_alloc(nmax, nsnapshot, NULL);

  while (fgets(buffer, sizeof(buffer), fp) != NULL)
    {
      int c;
      double epoch;

      if (*buffer == '#' || *buffer == '\0')
        {
          read_block = 0;
          continue;
        }

      c = sscanf(buffer, "M %lf", &epoch);
      if (c == 1)
        {
          /* start of new coefficient block */

          if (epoch < 1945.0)
            {
              fprintf(stderr, "msynth_bggm_read: error: epoch = %f\n", epoch);
              continue;
            }

          n = 1;
          m = 0;
          w->epochs[nepoch++] = epoch;
          read_block = 1;

          if (nepoch >= nsnapshot)
            {
              fprintf(stderr, "msynth_bggm_read: error: nsnapshot too small: %zu\n", nsnapshot);
              break;
            }
        }
      else if (read_block == 1)
        {
          double a[9];
          size_t i;
          size_t epoch_idx = nepoch - 1;

          c = sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                     &a[0],
                     &a[1],
                     &a[2],
                     &a[3],
                     &a[4],
                     &a[5],
                     &a[6],
                     &a[7],
                     &a[8]);
          if (c < 9)
            {
              read_block = 0;
              continue;
            }


          if (n > nmax)
            {
              fprintf(stderr, "msynth_bggm_read: error: nmax too small: %zu\n", nmax);
              break;
            }

          for (i = 0; i < 9; ++i)
            {
              size_t cidx = msynth_nmidx(n, m, w);
              w->c[epoch_idx * w->p + cidx] = a[i];

              /* update n, m*/

              if (m == 0)
                ++m;
              else if (m > 0)
                m = -m; /* just read g(n,m), next is h(n,m) */
              else
                {
                  m = -m;
                  ++m;
                }

              if (abs(m) > (int) n)
                {
                  /* read in next SH degree coefficients */
                  ++n;
                  m = 0;
                }
            }
        }
    } /* while (fgets(buffer, sizeof(buffer), fp) != NULL) */

  fclose(fp);

  if (nepoch == 0)
    {
      fprintf(stderr, "msynth_bggm_read: error parsing file %s\n", filename);
      return w;
    }

  w->n_epochs = nepoch;

  /* now compute secular variations from main field coefficients */
  for (i = 0; i < nepoch - 1; ++i)
    {
      double *cptr = w->c + i * w->p;
      double *cptr_next = w->c + (i + 1) * w->p;
      double dt = w->epochs[i + 1] - w->epochs[i];
      size_t j;

      for (j = 0; j < w->nnm; ++j)
        cptr[j + w->sv_offset] = (cptr_next[j] - cptr[j]) / dt;
    }

  return w;
}
