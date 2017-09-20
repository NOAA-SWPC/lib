/*
 * msynth_grid.h
 */

#ifndef INCLUDED_msynth_grid_h
#define INCLUDED_msynth_grid_h

#include <complex.h>

#ifndef INCLUDED_msynth_h
#include "msynth.h"
#define INCLUDED_msynth_h
#endif

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

typedef struct
{
  size_t nphi;                /* size of transform (number of longitude grid points) */
  size_t nmin;                /* minimum spherical harmonic degree */
  size_t nmax;                /* maximum spherical harmonic degree */
  complex double *fnm;        /* complex gauss coefficients */
  complex double *fnm_mirror; /* (-1)^{n+m} * f_{nm}, for computing mirror points */
  double *fft_r;              /* real output data, r component */
  double *fft_t;              /* real output data, theta component */
  double *fft_p;              /* real output data, phi component */
  double *rterm;              /* rterm[n] = (a/r)^{n+2}, size nmax + 1 */
  fftw_plan fftw_p;
  msynth_workspace *msynth_workspace_p;
} msynth_grid_workspace;

/*
 * Prototypes
 */

msynth_grid_workspace * msynth_grid_alloc(const size_t nphi, const size_t nmin, const size_t nmax,
                        const double * g);
void msynth_grid_free(msynth_grid_workspace *w);
int msynth_grid_init_r(const double r, msynth_grid_workspace * w);
int msynth_grid_init_theta(const double theta, msynth_grid_workspace * w);
int msynth_grid_calc(const int mirror, const double r, const double theta, gsl_matrix *B,
                     msynth_grid_workspace *w);

#endif /* INCLUDED_msynth_grid_h */
