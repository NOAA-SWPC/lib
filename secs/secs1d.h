/*
 * secs1d.h
 */

#ifndef INCLUDED_secs1d_h
#define INCLUDED_secs1d_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>

#ifndef INCLUDED_track_h
#include "track.h"
#define INCLUDED_track_h
#endif

/* mu_0 in units of: nT / (kA km^{-1}) */
#define SECS1D_MU_0                  (400.0 * M_PI)

/* recommended lmax for secs1d_alloc() */
#define SECS1D_LMAX                  200

/* maximum QD latitude for fitting SECS */
#define SECS1D_QDMAX                 (40.0)

/* flags */
#define SECS1D_FLG_FIT_DF            (1 << 0) /* fit divergence-free SECS */
#define SECS1D_FLG_FIT_CF            (1 << 1) /* fit curl-free SECS */

typedef struct
{
  size_t n;         /* total number of measurements in system */
  size_t p;         /* p_df + p_cf */
  size_t nmax;      /* maximum number of measurements in LS system */
  size_t flags;     /* SECS1D_FLG_xxx */
  size_t npoles;    /* number of poles */
  double dtheta;    /* spacing of 1D SECS poles in radians */

  size_t df_offset; /* offset in 'c' of DF coefficients */
  size_t cf_offset; /* offset in 'c' of CF coefficients */

  double R_iono;    /* radius of ionosphere (km) */

  size_t lmax;      /* maximum order for legendre functions */
  double *theta0;   /* pole locations */
  double *Pl1theta; /* P_{l,1}(cos(theta)) legendre functions, size lmax + 1 */
  gsl_matrix *Pltheta0; /* P_l(cos(theta0)) functions, npoles-by-(lmax + 1) */
  double *Pltheta;  /* P_l(cos(theta)) functions */

  gsl_matrix *X;    /* LS matrix */
  gsl_vector *rhs;  /* rhs vector */
  gsl_vector *wts;  /* weight vector */
  gsl_matrix *cov;  /* covariance matrix */
  gsl_matrix *L;    /* regularization matrix */
  gsl_vector *Ltau; /* regularization matrix Householder scalars */
  gsl_matrix *M;    /* workspace matrix */
  gsl_matrix *Xs;   /* standard form X */
  gsl_vector *bs;   /* standard form b */
  gsl_vector *cs;   /* standard form c */

  /*
   * solution vector is organized as:
   * c = [ c_df ; c_cf ]
   * where c_df and c_cf are both length 'npoles'
   */
  gsl_vector *c;    /* solution vector */

  gsl_multifit_linear_workspace *multifit_p;
} secs1d_workspace;

/*
 * Prototypes
 */
secs1d_workspace *secs1d_alloc(const size_t flags, const size_t lmax,
                               const double R_iono, const double pole_spacing);
void secs1d_free(secs1d_workspace *w);
int secs1d_add_track(const track_data *tptr, const satdata_mag *data,
                     secs1d_workspace *w);
int secs1d_fit(secs1d_workspace *w);
int secs1d_eval_B(const double r, const double theta,
                  double B[3], secs1d_workspace *w);
int secs1d_eval_J(const double r, const double theta,
                  double J[3], secs1d_workspace *w);
int secs1d_print_track(const int header, FILE *fp, const track_data *tptr,
                       const satdata_mag *data, secs1d_workspace *w);
int secs1d_green_df_init(const double theta, secs1d_workspace *w);
int secs1d_green_df(const double r, const double theta, const size_t pole_idx,
                    double B[3], secs1d_workspace *w);
int secs1d_green_cf(const double r, const double theta, const double theta0,
                    double B[3], secs1d_workspace *w);
int secs1d_green_df_J(const double theta, const double theta0, double K[3], secs1d_workspace *w);
int secs1d_green_cf_J(const double r, const double theta, const double theta0, double K[3], secs1d_workspace *w);

#endif /* INCLUDED_secs1d_h */
