/*
 * curvefit.h
 */

#ifndef INCLUDED_curvefit_h
#define INCLUDED_curvefit_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>

typedef struct
{
  const char *name;
  size_t size;
  int (*alloc) (void *vstate, const size_t n, const size_t p);
  int (*init) (void *vstate, const double *x, const double *y);
  int (*design_row) (void *vstate, const double x, gsl_vector *v);
  void (*free) (void *vstate);
} curvefit_type;

typedef struct
{
  gsl_multifit_robust_workspace *robust_p;
  gsl_matrix *A;   /* design matrix */
  gsl_vector *f;   /* RHS vector */
  gsl_vector *c;   /* polynomial coefficients */
  gsl_matrix *cov; /* covariance matrix */
  gsl_vector *x;   /* standardized x values */
  gsl_vector *r;   /* residuals */
  size_t n;        /* number of observations */
  size_t p;        /* number of coefficients */
  const curvefit_type *type;
  double chisq;
  double mean;     /* mean for center/scaling */
  double sigma;    /* standard deviation for center/scaling */
  void *state;
} curvefit_workspace;

/*
 * Prototypes
 */
curvefit_workspace *curvefit_alloc(const gsl_multifit_robust_type *robust_T,
                                   const curvefit_type *curve_T,
                                   const size_t n, const size_t p);
void curvefit_free(curvefit_workspace *w);
int curvefit(const int standardize, const double *x, const double *y,
             curvefit_workspace *w);
int curvefit_init(const int standardize, const double *x, const double *y,
                  curvefit_workspace *w);
int curvefit_solve(curvefit_workspace *w);
double curvefit_eval(const double x, curvefit_workspace *w);
double curvefit_residual(const curvefit_workspace *w);

extern const curvefit_type * curvefit_bspline;
extern const curvefit_type * curvefit_poly;

#endif /* INCLUDED_curvefit_h */
