/*
 * histogramNd.h
 * Patrick Alken
 */

#ifndef INCLUDED_histogramNd_h
#define INCLUDED_histogramNd_h

struct histogramNd_workspace
{
  size_t n;          /* number of dimensions */
  size_t *N;         /* number of bins in each dimension */
  double *bins;      /* bin data */
  size_t total_bins; /* total number of bins */
  
  size_t *loc;       /* location of an element in histogram */

  double **range;    /* ranges of each bin for each dimension */
};

/*
 * Prototypes
 */

struct histogramNd_workspace *histogramNd_alloc(size_t n, size_t *N);
void histogramNd_free(struct histogramNd_workspace *h);
void histogramNd_set_ranges_uniform(struct histogramNd_workspace *h,
                                    double *min, double *max);
void histogramNd_increment(struct histogramNd_workspace *h, double *x);
void histogramNd_accumulate(struct histogramNd_workspace *h, double *x,
                            double weight);
int histogramNd_get_by_x(struct histogramNd_workspace *h,
                         double *x, double *w);

#endif /* INCLUDED_histogramNd_h */
