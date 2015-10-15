/*
 * msynth_grid.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "msynth.h"
#include "msynth_grid.h"

/*
msynth_grid_alloc()
  Allocate msynth grid workspace

Inputs: nphi     - number of equally spaced longitude nodes from [0,2pi]
        msynth_p - msynth workspace containing model coefficients

Return: pointer to new workspace
*/

msynth_grid_workspace *
msynth_grid_alloc(const size_t nphi, const msynth_workspace *msynth_p)
{
  msynth_grid_workspace *w;

  w = calloc(1, sizeof(msynth_grid_workspace));
  if (!w)
    return 0;

  w->msynth_workspace_p = msynth_p;

  return w;
} /* msynth_grid_alloc() */

void
msynth_grid_free(msynth_grid_workspace *w)
{
  free(w);
} /* msynth_grid_free() */

/*
msynth_grid_calc()
  Synthesize magnetic field values at all longitude nodes for
a given r,theta

Inputs: r     - radius (km)
        theta - colatitude (radians)
        B     - (output) nphi-by-3 matrix
                B(:,1) = B_x (nT)
                B(:,2) = B_y (nT)
                B(:,3) = B_z (nT)
        w     - workspace

Return: success or error
*/

int
msynth_grid_calc(const double r, const double theta, gsl_matrix *B,
                 msynth_grid_workspace *w)
{
  int s = 0;

  return s;
} /* msynth_grid_calc() */
