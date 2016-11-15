/*
 * secs1d.h
 */

#ifndef INCLUDED_secs1d_h
#define INCLUDED_secs1d_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

/* mu_0 in units of: nT / (kA km^{-1}) */
#define SECS1D_MU_0                  (400.0 * M_PI)

/* recommended lmax for secs1d_alloc() */
#define SECS1D_LMAX                  200

typedef struct
{
  size_t n;         /* total number of measurements in system */
  size_t p;         /* p_df + p_cf */
  size_t p_df;      /* number of divergence free 1D SECS */
  size_t p_cf;      /* number of curl free 1D SECS */
  size_t nmax;      /* maximum number of measurements in LS system */

  double R_iono;    /* radius of ionosphere (km) */

  size_t lmax;      /* maximum order for legendre functions */
  double *theta0;   /* pole locations */
  double *Pl1;      /* P_{l,1}(cos(theta)) legendre functions */
  double *Pltheta0; /* P_l(cos(theta0)) functions */
  double *Pltheta;  /* P_l(cos(theta)) functions */

  gsl_matrix *X;    /* LS matrix */
  gsl_vector *c;    /* solution vector */
  gsl_vector *rhs;  /* rhs vector */
  gsl_vector *wts;  /* weight vector */
} secs1d_workspace;

/*
 * Prototypes
 */
secs1d_workspace *secs1d_alloc(const size_t lmax, const double R_iono, const size_t npoles);
void secs1d_free(secs1d_workspace *w);
int secs1d_df_green(const double r, const double theta, const double theta0,
                    double B[3], secs1d_workspace *w);
int secs1d_df_green_J(const double theta, const double theta0, double K[3], secs1d_workspace *w);

#endif /* INCLUDED_secs1d_h */
