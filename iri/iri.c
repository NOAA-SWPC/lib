/*
 * iri.c
 * Patrick Alken
 *
 * Wrapper routines for IRI model
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_math.h>

#include <indices/indices.h>

#include "iri.h"

static int iri_call(float latitude, float longitude, int year, int month,
                    int day, float hour, float hbeg, float hstp,
                    size_t nstp, float f107, iri_workspace *w);

iri_workspace *
iri_alloc(size_t nalt, const char *f107_datadir)
{
  iri_workspace *w;
  size_t i;

  w = calloc(1, sizeof(iri_workspace));
  if (!w)
    {
      fprintf(stderr, "iri_alloc: calloc failed: %s\n", strerror(errno));
      return 0;
    }

  w->iri_results = malloc(nalt * sizeof(iri_result));
  if (!w->iri_results)
    {
      fprintf(stderr, "iri_alloc: malloc failed: %s\n", strerror(errno));
      iri_free(w);
      return 0;
    }

  w->f107_workspace_p = f107_alloc(f107_datadir);
  if (!w->f107_workspace_p)
    {
      fprintf(stderr, "iri_alloc: f107_alloc failed\n");
      iri_free(w);
      return 0;
    }

  w->nalt = nalt;

  /* initialize IRI model */
  /*initialize_();*/

  /* initialize output arrays */
  for (i = 0; i < IRI_NR_OUTF*IRI_NC_OUTF; ++i)
    w->outf[i] = 0.0;
  for (i = 0; i < IRI_N_OARR; ++i)
    w->oarr[i] = 0.0;

  /* initialize flags parameter */

  for (i = 0; i < IRI_N_FLAGS; ++i)
    w->jf[i] = 1;

  w->jf[21] = 0; /* ion densities in /m^3 */
  w->jf[24] = 0; /* use F10.7 in oarr[40] */
  w->jf[25] = 0; /* no storm updating */
  w->jf[31] = 0; /* use F10.7A in oarr[45] */
  w->jf[32] = 0; /* turn off auroral boundary model */

  /* 0 = geographic, 1 = geomagnetic coords */
  w->jmag = 0;

  w->f107_override = -1.0;

  return w;
} /* iri_alloc() */

void
iri_free(iri_workspace *w)
{
  if (w->iri_results)
    free(w->iri_results);

  if (w->f107_workspace_p)
    f107_free(w->f107_workspace_p);

  free(w);
} /* iri_free() */

void
iri_f107_override(double f107, iri_workspace *w)
{
  w->f107_override = f107;
} /* iri_f107_override() */

/*
iri_calc()
  Calculate IRI parameters

Inputs: theta  - geographic colatitude (radians)
        phi    - geographic longitude (radians)
        t      - timestamp
        altmin - minimum altitude in km
        altstp - altitude step size in km
        nalt   - number of altitude steps
        w      - iri workspace

Return: success or error

Notes: IRI parameters are stored in w->iri_results

For altitudes below about 120km, IRI fails to find Tn, Ti, Te values
(they are set to -1). A similar problem exists for the densities, but
it appears safe to set them to 0 at these altitudes since they are
very small. It is not clear that the temperatures are small enough to
set to 0 however.
*/

int
iri_calc(double theta, double phi, time_t t, double altmin,
         double altstp, size_t nalt, iri_workspace *w)
{
  int s = 0;
  const int n_iri = 100; /* maximum altitudes IRI can compute at once */
  size_t i, j;
  int year, month, day;
  float hour;
  float flat = 90.0 - theta * 180.0 / M_PI,
        flon = phi * 180.0 / M_PI;
  double f107;
  struct tm *tm_p;
  float *outfv = w->outf;
  size_t nblocks; /* number of altitude blocks of n_iri steps */
  float hbeg;

  if (nalt > w->nalt)
    {
      fprintf(stderr, "iri_calc: specified nalt exceeds previous value (%zu,%zu)\n",
              nalt, w->nalt);
      return 1;
    }

  if (w->f107_override > 0.0)
    f107 = w->f107_override;
  else
    f107_get(t, &f107, w->f107_workspace_p);

  putenv("TZ=GMT");
  tm_p = gmtime(&t);

  year = tm_p->tm_year + 1900;
  month = tm_p->tm_mon + 1;
  day = tm_p->tm_mday;
  hour = (float) tm_p->tm_hour +
         (float) tm_p->tm_min / 60.0 +
         (float) tm_p->tm_sec / 3600.0;

  nblocks = (nalt - 1) / n_iri + 1;
  hbeg = altmin;

  for (i = 0; i < nblocks; ++i)
    {
      size_t nstp;
      
      if ((i + 1) * n_iri < nalt)
        nstp = n_iri;
      else
        nstp = nalt - i * n_iri;

      iri_call(flat,
               flon,
               year,
               month,
               day,
               hour,
               hbeg,
               (float)altstp,
               nstp,
               (float) f107,
               w);

      hbeg += (float) (altstp * nstp);

      /* save results into w->iri_results */

      for (j = 0; j < nstp; ++j)
        {
          size_t idx = i * n_iri + j;

          w->iri_results[idx].Ne = (double) outfv[IRI_FIDX(0, j)];
          w->iri_results[idx].Tn = (double) outfv[IRI_FIDX(1, j)];
          w->iri_results[idx].Ti = (double) outfv[IRI_FIDX(2, j)];
          w->iri_results[idx].Te = (double) outfv[IRI_FIDX(3, j)];

          if (w->iri_results[idx].Ne < 0.0)
            {
              /*
               * Ne is very small at these altitudes (< 80km) so just
               * set it to 0
               */
              w->iri_results[idx].Ne = 0.0;
            }

          w->iri_results[idx].n_Op = outfv[IRI_FIDX(4, j)];   /* O+ */
          w->iri_results[idx].n_Hp = outfv[IRI_FIDX(5, j)];   /* H+ */
          w->iri_results[idx].n_HEp = outfv[IRI_FIDX(6, j)];  /* HE+ */
          w->iri_results[idx].n_O2p = outfv[IRI_FIDX(7, j)];  /* O2+ */
          w->iri_results[idx].n_NOp = outfv[IRI_FIDX(8, j)];  /* NO+ */
          w->iri_results[idx].n_clus = outfv[IRI_FIDX(9, j)]; /* cluster */
          w->iri_results[idx].n_Np = outfv[IRI_FIDX(10, j)];  /* N+ */

          if (w->iri_results[idx].n_Op < 0.0)
            w->iri_results[idx].n_Op = 0.0;

          if (w->iri_results[idx].n_Hp < 0.0)
            w->iri_results[idx].n_Hp = 0.0;

          if (w->iri_results[idx].n_HEp < 0.0)
            w->iri_results[idx].n_HEp = 0.0;

          if (w->iri_results[idx].n_O2p < 0.0)
            w->iri_results[idx].n_O2p = 0.0;

          if (w->iri_results[idx].n_NOp < 0.0)
            w->iri_results[idx].n_NOp = 0.0;

          if (w->iri_results[idx].n_clus < 0.0)
            w->iri_results[idx].n_clus = 0.0;

          if (w->iri_results[idx].n_Np < 0.0)
            w->iri_results[idx].n_Np = 0.0;
        }
    }

  return s;
} /* iri_calc() */

/*
iri_get_result()
*/

iri_result *
iri_get_result(size_t idx, iri_workspace *w)
{
  return (&(w->iri_results[idx]));
} /* iri_get_result() */

/*
iri_call()
*/

static int
iri_call(float latitude, float longitude, int year, int month, int day,
         float hour, float hbeg, float hstp, size_t nstep, float f107,
         iri_workspace *w)
{
  int s = 0;
  float *outfv = w->outf;
  float *oarrv = w->oarr;
  float hend;
  int mmdd;

  hend = hbeg + nstep * hstp;

  if (longitude < 0.0)
    longitude += 360.0;

  hour += 25.0;

  oarrv[40] = f107;
  oarrv[45] = f107;

  mmdd = month * 100 + day;

  iri_sub_(w->jf,
           &(w->jmag),
           &latitude,
           &longitude,
           &year,
           &mmdd,
           &hour,
           &hbeg,
           &hend,
           &hstp,
           outfv,
           oarrv);

  return s;
} /* iri_call() */
