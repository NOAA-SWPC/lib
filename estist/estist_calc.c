/*
 * estist_calc.c
 *
 * Convert Dst into external (Est) and internal (Ist) parts
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <indices/indices.h>

#include "estist_calc.h"

estist_calc_workspace *
estist_calc_alloc()
{
  estist_calc_workspace *w;

  w = calloc(1, sizeof(estist_calc_workspace));

  w->dst_workspace_p = dst_alloc(DST_IDX_FILE);

  w->model = 1;

  return w;
} /* estist_calc_alloc() */

void
estist_calc_free(estist_calc_workspace *w)
{
  if (w->dst_workspace_p)
    dst_free(w->dst_workspace_p);

  free(w);
} /* estist_calc_free() */

/*
estist_calc()
  Deparate a DST time series into Est/Ist components

Inputs: ndst - number of DST measurements
        dst  - array of DST time series
        est  - (output) where to store Est time series
        ist  - (output) where to store Ist time series
        w    - workspace

Notes: about 1 year of Dst data is needed to get an accurate separation
of Est/Ist. Since Dst measurements are provided every hour, ndst should
therefore be at least:

ndst >= 365 (days/yr) * 24 (hrs/day)
*/

int
estist_calc(int ndst, double dst[], double est[], double ist[],
            estist_calc_workspace *w)
{
  int s = 0;

  weidelt_dst_(&(w->model),
               &ndst,
               dst,
               est,
               ist);

  return s;
} /* estist_calc() */

/*
estist_calc_get()
  Compute Est/Ist for a given timestamp
*/

int
estist_calc_get(const time_t t, double *est, double *ist,
                estist_calc_workspace *w)
{
  int s = 0;
  const size_t ndst = ESTIST_NDST;
  size_t i;
  time_t t1 = t - (ndst - 1) * 3600; /* start time for Dst series */

  for (i = 0; i < ndst; ++i)
    {
      s += dst_get(t1, &(w->dst[i]), w->dst_workspace_p);
      if (s)
        {
          fprintf(stderr, "estist_calc_get: dst not found: %ld\n", t1);
          return s;
        }

      /* advance 1 hour */
      t1 += 3600;
    }

  /* perform Est/Ist separation */
  s += estist_calc(ndst, w->dst, w->est, w->ist, w);

  *est = w->est[ndst - 1];
  *ist = w->ist[ndst - 1];

  return s;
} /* estist_calc_get() */
