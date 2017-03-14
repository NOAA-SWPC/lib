/*
 * msynth_swarm.c
 * 
 * Routines to read a SWARM SHC formatted coefficient file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include "msynth.h"

/*
msynth_swarm_read()
  Read SWARM coefficient file
*/

msynth_workspace *
msynth_swarm_read(const char *filename)
{
  FILE *fp;
  size_t nmin = 0;
  size_t nmax = 0;
  size_t num_epochs = 0;
  int gotheader = 0;
  int gotepochs = 0;
  msynth_workspace *w = NULL;
  char buffer[8196];
  double epochs[MSYNTH_MAX_SNAPSHOT];
  size_t i;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "msynth_swarm_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  while (fgets(buffer, sizeof(buffer), fp) != NULL)
    {
      int c;
      int offset = 0;
      char *str;

      if (*buffer == '#')
        continue;

      if (!gotheader)
        {
          size_t dummy;

          c = sscanf(buffer, "%zu %zu %zu %zu %zu",
                     &nmin,
                     &nmax,
                     &num_epochs,
                     &dummy,
                     &dummy);
          if (c < 5)
            continue;

          gotheader = 1;
        }
      else if (!gotepochs)
        {
          double epoch;

          /* read in epochs */
          str = buffer + offset;
          i = 0;
          while ((c = sscanf(str, "%lf %n", &epoch, &offset)) == 1)
            {
              epochs[i] = epoch;
              str += offset;

              ++i;

              if (i == num_epochs)
                break; /* done */

              if (i >= MSYNTH_MAX_SNAPSHOT)
                {
                  fprintf(stderr, "msynth_swarm_read: MSYNTH_MAX_SNAPSHOT not large enough\n");
                  return 0;
                }
            }

          if (i != num_epochs)
            {
              fprintf(stderr, "msynth_swarm_read: error reading epochs\n");
              return 0;
            }

          /* allocate workspace */
          w = msynth_alloc(nmax, num_epochs, epochs);
          msynth_set(nmin, nmax, w);

          gotepochs = 1;
        }
      else
        {
          size_t n;
          int m;
          size_t cidx;
          double gnm;

          assert(w != NULL);

          str = buffer + offset;
          c = sscanf(buffer, "%zu %d %n", &n, &m, &offset);
          if (c < 2)
            continue;

          cidx = msynth_nmidx(n, m, w);

          /* now read in coefficients for all snapshots */
          str = buffer + offset;
          i = 0;
          while ((c = sscanf(str, "%lf %n", &gnm, &offset)) == 1)
            {
              double *cptr = w->c + i * w->p;

              cptr[cidx] = gnm;
              cptr[cidx + w->sv_offset] = 0.0; /* computed later */
              cptr[cidx + w->sa_offset] = 0.0;

              str += offset;
              ++i;
            }

          if (i != num_epochs)
            {
              fprintf(stderr, "msynth_swarm_read: error in file format\n");
              break;
            }
        }
    } /* while (fgets(buffer, sizeof(buffer), fp) != NULL) */

  fclose(fp);

  if (num_epochs == 0)
    {
      fprintf(stderr, "msynth_swarm_read: error parsing file %s\n", filename);
      return w;
    }

  /* now compute secular variations from main field coefficients */
  for (i = 0; i < num_epochs - 1; ++i)
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
