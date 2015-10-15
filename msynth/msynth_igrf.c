/*
 * msynth_igrf.c
 * 
 * Routines to read a IGRF coefficient file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include "msynth.h"

/*
msynth_igrf_read()
  Read IGRF coefficient file
*/

msynth_workspace *
msynth_igrf_read(const char *filename)
{
  FILE *fp;
  const size_t nmax = 13;
  size_t num_epochs = 0;
  msynth_workspace *w = NULL;
  char buffer[2048];
  double epochs[MSYNTH_MAX_SNAPSHOT];
  size_t i;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "msynth_igrf_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  while (fgets(buffer, sizeof(buffer), fp) != NULL)
    {
      size_t n;
      int c, m;
      char gh;
      size_t cidx;
      int offset;
      char *str;

      if (num_epochs == 0)
        {
          char nc, mc;
          int epoch, prev_epoch = 0;

          c = sscanf(buffer, "g/h %c %c %n", &nc, &mc, &offset);
          if (c < 2 || nc != 'n' || mc != 'm')
            continue;

          /* read in epochs */
          str = buffer + offset;
          while ((c = sscanf(str, "%04d.0 %n", &epoch, &offset)) == 1)
            {
              if (prev_epoch == epoch)
                break; /* end of the line (SV column) */

              epochs[num_epochs++] = (double) epoch;
              str += offset;
              prev_epoch = epoch;
            }

          /* allocate workspace */
          w = msynth_alloc(nmax, num_epochs, epochs);

          continue;
        }

      assert(w != NULL);

      c = sscanf(buffer, "%c %zu %d %n", &gh, &n, &m, &offset);
      if (c < 3)
        continue;

      if (gh == 'g')
        cidx = msynth_nmidx(n, m, w);
      else if (gh == 'h')
        cidx = msynth_nmidx(n, -m, w);
      else
        continue;

      /* read in main field coefficients for all epochs */
      str = buffer + offset;
      for (i = 0; i < num_epochs; ++i)
        {
          double gnm;
          double *cptr = w->c + i * w->p;

          c = sscanf(str, "%lf %n", &gnm, &offset);
          if (c < 1)
            {
              fprintf(stderr, "msynth_igrf_read: error in file format\n");
              break;
            }

          str += offset;

          /* store coefficient */
          cptr[cidx] = gnm;
          cptr[cidx + w->sa_offset] = 0.0;

          if (i == num_epochs - 1)
            {
              double dgnm;

              c = sscanf(str, "%lf %n", &dgnm, &offset);
              if (c < 1)
                {
                  fprintf(stderr, "msynth_igrf_read: error in file format\n");
                  break;
                }

              cptr[cidx + w->sv_offset] = dgnm;
            }
          else
            cptr[cidx + w->sv_offset] = 0.0; /* computed later */
        }
    } /* while (fgets(buffer, sizeof(buffer), fp) != NULL) */

  fclose(fp);

  if (num_epochs == 0)
    {
      fprintf(stderr, "msynth_igrf_read: error parsing file %s\n", filename);
      return w;
    }

  /* now compute secular variations from main field coefficients */
  for (i = 0; i < num_epochs - 1; ++i)
    {
      double *cptr = w->c + i * w->p;
      double *cptr_next = w->c + (i + 1) * w->p;
      size_t j;

      for (j = 0; j < w->nnm; ++j)
        cptr[j + w->sv_offset] = (cptr_next[j] - cptr[j]) / 5.0;
    }

  return w;
} /* msynth_igrf_read() */

/* read IGRF-12 candidate main field model */
msynth_workspace *
msynth_igrf_read_mf(const char *filename)
{
  FILE *fp;
  const size_t nmax = 13;
  const double epoch = 2015.0;
  msynth_workspace *w = msynth_alloc(nmax, 1, &epoch);
  char buffer[2048];

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "msynth_igrf_read_mf: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  while (fgets(buffer, 2048, fp) != NULL)
    {
      int n, m, c;
      double gnm, hnm;
      size_t cidx;

      c = sscanf(buffer, "%d %d %lf %lf",
                 &n, &m, &gnm, &hnm);
      if (c < 4)
        continue;

      assert(n <= 13);
      assert(m <= n);

      cidx = msynth_nmidx(n, m, w);
      w->c[cidx] = gnm;
      w->c[cidx + w->sv_offset] = 0.0;
      w->c[cidx + w->sa_offset] = 0.0;

      if (m != 0)
        {
          cidx = msynth_nmidx(n, -m, w);
          w->c[cidx] = hnm;
          w->c[cidx + w->sv_offset] = 0.0;
          w->c[cidx + w->sa_offset] = 0.0;
        }
    }

  fclose(fp);

  return w;
} /* msynth_igrf_read_mf() */

/* read IGRF-12 candidate SV model */
int
msynth_igrf_read_sv(const char *filename, msynth_workspace *w)
{
  FILE *fp;
  char buffer[2048];

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "msynth_igrf_read_sv: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  while (fgets(buffer, 2048, fp) != NULL)
    {
      int n, m, c;
      double dgnm, dhnm;
      size_t cidx;

      c = sscanf(buffer, "%d %d %lf %lf",
                 &n, &m, &dgnm, &dhnm);
      if (c < 4)
        continue;

      assert(n <= 13);
      assert(m <= n);

      cidx = msynth_nmidx(n, m, w);
      w->c[cidx + w->sv_offset] = dgnm;

      if (m != 0)
        {
          cidx = msynth_nmidx(n, -m, w);
          w->c[cidx + w->sv_offset] = dhnm;
        }
    }

  fclose(fp);

  return 0;
} /* msynth_igrf_read_sv() */

/* write candidate main field model in requested format */
int
msynth_igrf_write(const char *filename, const char *desc,
                  const msynth_workspace *w)
{
  const size_t nmax = 13;
  size_t n;
  int m;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "msynth_igrf_write: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  fprintf(fp, "# National Geophysical Data Center, Boulder, CO\n");
  fprintf(fp, "# Candidate for %s\n", desc);
  fprintf(fp, "# %c  %c %9s %9s %20s %20s\n",
          'n',
          'm',
          "gnm",
          "hnm",
          "uncertainty_gnm",
          "uncertainty_hnm");

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;

      for (m = 0; m <= ni; ++m)
        {
          double gnm, hnm;
          size_t cidx = msynth_nmidx(n, m, w);

          gnm = w->c[cidx];

          if (m == 0)
            hnm = 0.0;
          else
            {
              cidx = msynth_nmidx(n, -m, w);
              hnm = w->c[cidx];
            }

          fprintf(fp, "%3zu %2d %9.2f %9.2f %20.2f %20.2f\n",
                  n,
                  m,
                  gnm,
                  hnm,
                  0.0,
                  0.0);
        }
    }

  fclose(fp);

  return 0;
} /* msynth_igrf_write() */

/* write candidate SV field model in requested format */
int
msynth_igrf_sv_write(const char *filename, const char *desc,
                     const msynth_workspace *w)
{
  const size_t nmax = 8;
  size_t n;
  int m;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "msynth_igrf_sv_write: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  fprintf(fp, "# National Geophysical Data Center, Boulder, CO\n");
  fprintf(fp, "# Candidate for %s\n", desc);
  fprintf(fp, "# %c  %c %9s %9s %20s %20s\n",
          'n',
          'm',
          "gnm",
          "hnm",
          "uncertainty_gnm",
          "uncertainty_hnm");

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;

      for (m = 0; m <= ni; ++m)
        {
          double dgnm, dhnm;
          size_t cidx = msynth_nmidx(n, m, w);

          dgnm = w->c[cidx + w->sv_offset];

          if (m == 0)
            dhnm = 0.0;
          else
            {
              cidx = msynth_nmidx(n, -m, w);
              dhnm = w->c[cidx + w->sv_offset];
            }

          fprintf(fp, "%3zu %2d %9.2f %9.2f %20.2f %20.2f\n",
                  n,
                  m,
                  dgnm,
                  dhnm,
                  0.0,
                  0.0);
        }
    }

  fclose(fp);

  return 0;
} /* msynth_igrf_sv_sv_write() */
