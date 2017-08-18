/*
 * msynth_grid.h
 */

#ifndef INCLUDED_msynth_grid_h
#define INCLUDED_msynth_grid_h

#ifndef INCLUDED_msynth_h
#include "msynth.h"
#define INCLUDED_msynth_h
#endif

#include <fftw3.h>

typedef struct
{
  double *fft_x;      /* real output data, X component */
  double *fft_y;      /* real output data, Y component */
  double *fft_z;      /* real output data, Z component */
  size_t nphi;        /* size of transform (number of longitude grid points) */
  size_t nmin;        /* minimum spherical harmonic degree */
  size_t nmax;        /* maximum spherical harmonic degree */
  fftw_plan fftw_p;
  const msynth_workspace *msynth_workspace_p;
} msynth_grid_workspace;

/*
 * Prototypes
 */

msynth_grid_workspace *msynth_grid_alloc(const size_t nphi, const msynth_workspace *msynth_p);
void msynth_grid_free(msynth_grid_workspace *w);
int msynth_grid_calc(const double r, const double theta, gsl_matrix *B,
                     msynth_grid_workspace *w);
int msynth_grid_calc2(const double r, const double theta, gsl_matrix *B,
                      msynth_grid_workspace *w);

#endif /* INCLUDED_msynth_grid_h */
