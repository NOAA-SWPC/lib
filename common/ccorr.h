/*
 * ccorr.h
 */

#ifndef INCLUDED_ccorr_h
#define INCLUDED_ccorr_h

#define CCORR_MAX       1000

typedef struct
{
  double u[CCORR_MAX];
  double v[CCORR_MAX];
  size_t n;
} ccorr_bin;

typedef struct
{
  size_t nx;
  double xmin;
  double xmax;

  ccorr_bin *bins;
  size_t nbins;
} ccorr_workspace;

/*
 * Prototypes
 */

ccorr_workspace *ccorr_alloc(const double xmin, const double xmax, const size_t nx);
void ccorr_free(ccorr_workspace *w);
int ccorr_add(const double x, const double u, const double v, ccorr_workspace *w);
double ccorr_r(const double x, ccorr_workspace *w);
size_t ccorr_n(const double x, ccorr_workspace *w);
int ccorr_xval(const size_t i, double *x, ccorr_workspace *w);
int ccorr_reset(ccorr_workspace *w);

#endif /* INCLUDED_ccorr_h */
