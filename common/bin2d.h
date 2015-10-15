/*
 * bin2d.h
 */

#ifndef INCLUDED_bin2d_h
#define INCLUDED_bin2d_h

#include <gsl/gsl_rstat.h>

#define MAX_DATA_PER_BIN          10000

typedef struct
{
  size_t nx;
  size_t ny;
  double xmin;
  double xmax;
  double ymin;
  double ymax;

  gsl_rstat_workspace **bins;
  size_t nbins;

  /* for correlations */
  double **z1;
  double **z2;
  size_t *n;
} bin2d_workspace;

#define CIDX2(i,Ni,j,Nj)      ((j) + (Nj)*(i))
#define BIN2D_IDX(x, y, w)    (CIDX2((x), (w)->nx, (y), (w)->ny))

/*
 * Prototypes
 */

bin2d_workspace *bin2d_alloc(double xmin, double xmax, size_t nx,
                             double ymin, double ymax, size_t ny);
void bin2d_free(bin2d_workspace *w);
int bin2d_add_element(const double x, const double y, const double data, bin2d_workspace *w);
double bin2d_mean(const double x, const double y, const bin2d_workspace *w);
double bin2d_median(const double x, const double y, const bin2d_workspace *w);
double bin2d_sd(const double x, const double y, const bin2d_workspace *w);
size_t bin2d_n(const double x, const double y, const bin2d_workspace *w);
int bin2d_xyval(const size_t i, const size_t j,
                double *x, double *y, const bin2d_workspace *w);
int bin2d_print(const char *filename, const bin2d_workspace *w);
int bin2d_add_element_corr(double x, double y, double data1, double data2, bin2d_workspace *w);
double bin2d_correlation(const double x, const double y, const bin2d_workspace *w);

#endif /* INCLUDED_bin2d_h */
