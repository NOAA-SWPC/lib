/*
 * peak.h
 */

#ifndef INCLUDED_peak_h
#define INCLUDED_peak_h

typedef struct
{
  gsl_vector *deriv;  /* first derivative of input data */
  gsl_vector *sderiv; /* smoothed first derivative of input data */
  size_t n;           /* number of data points to process */
  size_t p;           /* number of model parameters for gaussian fit */
  size_t pidx;        /* index into data arrays of center of current peak */
  size_t idx;         /* index into data arrays of current search location */

  gsl_vector *c;      /* parameters of gaussian fit */
  size_t idx0;        /* starting index of data for gaussian fit */
  size_t idx1;        /* ending index of data for gaussian fit */
} peak_workspace;

/*
 * Prototypes
 */

peak_workspace *peak_alloc(const size_t n);
void peak_free(peak_workspace *w);
int peak_init(const size_t smooth_window, const gsl_vector *x,
              const gsl_vector *y, peak_workspace *w);
int peak_find(const int minmax, const double minslope, const double minheight,
              const gsl_vector *x, const gsl_vector *y, peak_workspace *w);
int peak_gaussian(const size_t fit_width, const gsl_vector *x,
                  const gsl_vector *y, peak_workspace *w);
double peak_eval(const gsl_vector *c, const double x);
double peak_deriv(const size_t i, const peak_workspace *w);
double peak_sderiv(const size_t i, const peak_workspace *w);

#endif /* INCLUDED_peak_h */
