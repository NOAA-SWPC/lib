/*
 * bin.h
 */

#ifndef INCLUDED_bin_h
#define INCLUDED_bin_h

#include <gsl/gsl_rstat.h>

#define MAX_DATA_PER_BIN          10000

typedef struct
{
  size_t nx;
  double xmin;
  double xmax;

  gsl_rstat_workspace **bins;
  size_t nbins;

  /* for correlations */
  double **z1;
  double **z2;
  size_t *n;
} bin_workspace;

/*
 * Prototypes
 */

bin_workspace *bin_alloc(const double xmin, const double xmax, const size_t nx);
void bin_free(bin_workspace *w);
int bin_add_element(const double x, const double y, bin_workspace *w);
int bin_add_element_corr(double x, double data1, double data2, bin_workspace *w);
double bin_correlation(const double x, const bin_workspace *w);
double bin_mean(const double x, bin_workspace *w);
double bin_sd(const double x, bin_workspace *w);
size_t bin_n(const double x, bin_workspace *w);
int bin_xval(const size_t i, double *x, bin_workspace *w);
int bin_reset(bin_workspace *w);

#endif /* INCLUDED_bin_h */
