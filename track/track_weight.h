/*
 * weight.h
 */

#ifndef INCLUDED_track_weight_h
#define INCLUDED_track_weight_h

#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_matrix.h>

#define WEIGHT_IDX(i,j,w)   (CIDX2((i), ((w)->nphi), (j), ((w)->ntheta)))

typedef struct
{
  size_t ntheta;
  size_t nphi;
  gsl_matrix *weight; /* weight grid */
  gsl_histogram2d *hist_p;
} track_weight_workspace;

/*
 * Prototypes
 */

track_weight_workspace *track_weight_alloc(const size_t ntheta, const size_t nphi);
void track_weight_free(track_weight_workspace *w);
int track_weight_reset(track_weight_workspace *w);
int track_weight_add_data(const double theta, const double phi, track_weight_workspace *w);
int track_weight_calc(track_weight_workspace *w);
int track_weight_get(const double phi, const double theta, double *weight, track_weight_workspace *w);
int track_weight_n(const double phi, const double theta, size_t *n, track_weight_workspace *w);
int track_weight_write(const char *filename, track_weight_workspace *w);

#endif /* INCLUDED_track_weight_h */
