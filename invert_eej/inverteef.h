/*
 * inverteef.h
 * Patrick Alken
 */

#ifndef INCLUDED_inverteef_h
#define INCLUDED_inverteef_h

#include <stdio.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>

/* define to invert for electric field */
#define PDE_INVERT_E_FIELD

/* define to allow a DC shift in inversion */
#define PDE_ALLOW_DC_SHIFT

/*
 * define to constrain inversion to have solution match satellite
 * profile at 90 degrees
 */
#define PDE_USE_LSE

#define PROF_LAT_MIN         (-20.0)
#define PROF_LAT_MAX         (20.0)
#define PROF_LAT_STEP        (0.5)
#define PROF_LAT_N           ((size_t) ((PROF_LAT_MAX - PROF_LAT_MIN) / PROF_LAT_STEP + 1))

#define PROF_THETA_MIN       (90.0 - PROF_LAT_MAX)
#define PROF_THETA_MAX       (90.0 - PROF_LAT_MIN)

typedef struct
{
  size_t ntheta;
  double theta_min;
  double theta_max;
} inverteef_parameters;

typedef struct
{
  size_t ntheta;
  double theta_min;
  double theta_max;
  double dtheta;

  gsl_vector *J_pde_E; /* corresponding pde profile for E_0 */
  gsl_vector *J_pde_u; /* corresponding pde profile for winds */
  gsl_vector *J_pde;  /* total PDE modeled profile */
  gsl_vector *J_rhs;  /* RHS for inversion */
  gsl_vector *J_diff; /* J_pde - J_sat */
  gsl_matrix *X;      /* least squares design matrix */
  gsl_vector *coeffs; /* least squares coefficients */
  gsl_matrix *cov;    /* covariance matrix */
  gsl_matrix *B;      /* LSE constraint matrix */
  gsl_vector *d;      /* LSE rhs */
  gsl_vector *r;      /* fit residuals */
  double chisq;
  double Rsq;
  double R;
  gsl_multifit_linear_workspace *multifit_p;
  double E_scale;     /* E field scaling factor */
  double J_DC;        /* DC current shift */
  double J_champ;     /* CHAMP peak current value */

  double RelErr;      /* relative error between model and satellite profile */
} inverteef_workspace;

/*
 * Prototypes
 */

inverteef_workspace *inverteef_alloc(inverteef_parameters *params);
void inverteef_free(inverteef_workspace *w);
int inverteef_calc(gsl_vector *J_sat, gsl_vector *J_lat_E,
                   gsl_vector *J_lat_u, inverteef_workspace *w);

#endif /* INCLUDED_inverteef_h */
