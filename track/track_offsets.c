/*
 * track_offsets.c
 *
 * Attempt to remove track offsets due to magnetospheric
 * contamination
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>

#include "curvefit.h"

#include "track.h"

static int track_calc_offset(const double qdlat_min, const double qdlat_max,
                             const double *qdlat, double *B, const size_t n);

int
track_fix_offsets(const satdata_mag *data, track_workspace *w)
{
  int s = 0;
  size_t i;
  const double qdlat_min = 25.0;
  const double qdlat_max = 60.0;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      double *qdlat = &(data->qdlat[tptr->start_idx]);

      /* ignore flagged tracks */
      if (tptr->flags)
        continue;

      track_calc_offset(qdlat_min, qdlat_max, qdlat, tptr->Bx, tptr->n);
      track_calc_offset(qdlat_min, qdlat_max, qdlat, tptr->Bz, tptr->n);
    }

  return s;
} /* track_calc_offsets() */

/*
track_calc_offset()
  Fit a robust polynomial to track data and subtract polynomial from
track. The polynomial is fit to data satisfying:

|qdlat| \in [qdlat_min,qdlat_max]

in order to exclude the equatorial region

Inputs: qdlat_min - minimum QD latitude (degrees)
        qldat_max - maximum QD latitude (degrees)
        qdlat     - array of QD latitudes (degrees)
        B         - (input/output) track magnetic field data (nT)
                    On output, the polynomial is subtracted from B
        n         - size of qdlat and B arrays

Return: success/error
*/

static int
track_calc_offset(const double qdlat_min, const double qdlat_max,
                  const double *qdlat, double *B, const size_t n)
{
  int s = 0;
  const gsl_multifit_robust_type *robust_t = gsl_multifit_robust_bisquare;
  const curvefit_type *curve_t = curvefit_poly;
  const size_t p = 1;
  curvefit_workspace *curvefit_p;
  double *x, *y;
  size_t npts = 0;
  size_t i;

  x = malloc(n * sizeof(double));
  y = malloc(n * sizeof(double));

  /* store data in arrays, excluding equatorial region */
  for (i = 0; i < n; ++i)
    {
      if (fabs(qdlat[i]) < qdlat_min || fabs(qdlat[i]) > qdlat_max)
        continue;

      x[npts] = qdlat[i];
      y[npts] = B[i];
      ++npts;
    }

  curvefit_p = curvefit_alloc(robust_t, curve_t, npts, p);

  curvefit(1, x, y, curvefit_p);

  for (i = 0; i < n; ++i)
    {
      double cval = curvefit_eval(qdlat[i], curvefit_p);
      B[i] -= cval;
    }

  curvefit_free(curvefit_p);
  free(x);
  free(y);

  return s;
} /* track_calc_offset() */
