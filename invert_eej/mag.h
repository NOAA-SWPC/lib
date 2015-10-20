/*
 * mag.h
 */

#ifndef INCLUDED_mag_h
#define INCLUDED_mag_h

#include <stdio.h>
#include <satdata/satdata.h>
#include <indices/indices.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>

#include "apex.h"
#include "green.h"
#include "log.h"
#include "inverteef.h"
#include "mageq.h"
#include "pde.h"
#include "track.h"

/* local time window for equator crossings */
#define MAG_LT_MIN             (6.0)
#define MAG_LT_MAX             (18.0)

/* where to store log files */
#define LOG_DIR                "log"

/* maximum latitude spacing to allow for data gaps */
#define MAG_DLAT               (2.0)

/* maximum kp index for track processing */
#define MAG_MAX_KP             (20.0)

/* maximum QD latitude for track processing */
#define MAG_MAX_QD_LATITUDE    (65.0)

/* maximum spherical harmonic degree for internal Sq filter model */
#define MAG_SQFILT_INT_DEG     (12)
#define MAG_SQFILT_INT_ORD     (0)

/* maximum spherical harmonic degree for external Sq filter model */
#define MAG_SQFILT_EXT_DEG     (1)
#define MAG_SQFILT_EXT_ORD     (1)

/* maximum data in one track */
#define MAG_MAX_TRACK          4000

/****************************
 * EEJ inversion parameters *
 ****************************/

/* number of current arcs representing EEJ */
#define MAG_EEJ_NCURR          81

/* qd latitude range for current flow */
#define MAG_EEJ_QDLAT_MIN      (-20.0)
#define MAG_EEJ_QDLAT_MAX      (20.0)

typedef struct
{
  double lt_min;                  /* minimum local time (hours) */
  double lt_max;                  /* maximum local time (hours) */
  double lon_min;                 /* minimum longitude (degrees) */
  double lon_max;                 /* maximum longitude (degrees) */
  int profiles_only;              /* compute profiles only (no EEF) */
} mag_params;

/* store 1 satellite track */
typedef struct
{
  double t[MAG_MAX_TRACK];        /* timestamp (CDF_EPOCH) */
  double theta[MAG_MAX_TRACK];    /* colatitude in radians */
  double phi[MAG_MAX_TRACK];      /* longitude in radians */
  double r[MAG_MAX_TRACK];        /* radius in km */
  double thetaq[MAG_MAX_TRACK];   /* QD colatitude in radians */
  double qdlat[MAG_MAX_TRACK];    /* QD latitude in degrees */
  double lat_deg[MAG_MAX_TRACK];  /* geocentric latitude in degrees */
  double F[MAG_MAX_TRACK];        /* scalar field measurement in nT */
  double Bx_int[MAG_MAX_TRACK];   /* Bx_internal NEC (nT) */
  double By_int[MAG_MAX_TRACK];   /* By_internal NEC (nT) */
  double Bz_int[MAG_MAX_TRACK];   /* Bz_internal NEC (nT) */
  double F_int[MAG_MAX_TRACK];    /* F_internal (nT) */
  double dF_ext[MAG_MAX_TRACK];   /* b_int . B_ext (nT) */
  double F1[MAG_MAX_TRACK];       /* F_sat - F_int - dF_ext (nT) */
  double F2[MAG_MAX_TRACK];       /* F^(1) - b . (M + K) (nT) */
  double Sq_int[MAG_MAX_TRACK];   /* internal Sq model b . M (nT) */
  double Sq_ext[MAG_MAX_TRACK];   /* external Sq model b . K (nT) */
  double F2_fit[MAG_MAX_TRACK];   /* fit to F^(2) from line current model */
  double phi_eq;                  /* longitude of equator crossing (rad) */
  double t_eq;                    /* time of equator crossing (CDF_EPOCH) */
  size_t n;                       /* number of data in this track */
} mag_track;

typedef struct
{
  double qd_alt;    /* altitude at which QD map is made (km) */
  size_t nlon;      /* number of longitude segments per current arc */
  size_t p;         /* number of model parameters (number of line currents) */

  double qdlat_min; /* minimum QD latitude for arc currents */
  double dqdlat;    /* distance in deg between arc currents */
  double curr_dist_km; /* approx. distance in km between arc currents */

  gsl_matrix *Jx;   /* current density greens functions */
  gsl_matrix *Jy;
  gsl_matrix *Jz;

  /* ECEF Cartesian positions of all current segment midpoints (km) */
  gsl_matrix *mid_pos_x;
  gsl_matrix *mid_pos_y;
  gsl_matrix *mid_pos_z;

  double *mid_lon; /* approx. geocentric longitudes of segment midpts (deg) */

  gsl_matrix *X;     /* least squares matrix */
  gsl_matrix *cov;   /* covariance matrix */
  gsl_vector *S;     /* current strength coefficients (S_j) */
  gsl_vector *rhs;   /* right hand side vector (F^(2)) */
  gsl_matrix *L;     /* regularization matrix */
  gsl_vector *tau;   /* for QR decomposition */
  double rnorm;      /* residual norm || y - X c || */
  double snorm;      /* solution norm || L c || */
  gsl_multifit_linear_workspace *multifit_p;

  /* parameters for L-curve analysis */
  size_t nreg;       /* number of regularization parameters on L-curve */
  gsl_vector *rho;   /* vector of residual norms */
  gsl_vector *eta;   /* vector of solution norms */
  gsl_vector *reg_param; /* vector of regularization parameters */
  size_t reg_idx;    /* index of optimal regularization parameter */

  apex_workspace *apex_workspace_p;
  mageq_workspace *mageq_workspace_p;
} mag_eej_workspace;

typedef struct
{
  gsl_matrix *X;     /* least squares matrix */
  gsl_matrix *cov;   /* covariance matrix */
  gsl_vector *c;     /* coefficient vector */
  gsl_vector *rhs;   /* right hand side vector */
  size_t p;          /* number of coefficients for Sq model */
  size_t n;          /* number of data for LS fit */
  size_t ntot;       /* total data size allocated */
  double rnorm;      /* residual norm || y - X c || */
  double snorm;      /* solution norm || L c || */

  size_t nmax_int;   /* maximum internal spherical harmonic degree */
  size_t mmax_int;   /* maximum internal spherical harmonic order */
  size_t nmax_ext;   /* maximum external spherical harmonic degree */
  size_t mmax_ext;   /* maximum external spherical harmonic order */
  size_t p_int;      /* number of internal coefficients */
  size_t p_ext;      /* number of external coefficients */
  size_t int_offset; /* offset in 'c' of internal coefficients */
  size_t ext_offset; /* offset in 'c' of external coefficients */

  /*
   * indexing for coefficients, base[n] gives the offset in c for
   * all coefficients of degree n
   */
  size_t *base_int;
  size_t *base_ext;

  gsl_vector *work_p; /* workspace of size p */

  size_t nreg;       /* number of regularization parameters for L-curve analysis */
  gsl_vector *rho;   /* vector of residual norms */
  gsl_vector *eta;   /* vector of solution norms */
  gsl_vector *reg_param; /* vector of regularization parameters */
  size_t reg_idx;    /* index of optimal regularization parameter */

  gsl_multifit_linear_workspace *multifit_workspace_p;
} mag_sqfilt_workspace;

typedef struct
{
  mag_track track;  /* stores current track for processing */

  double EEJ[MAG_EEJ_NCURR]; /* EEJ current density (A/m) */

  double EEF;                /* final EEF value (mV/m) */
  double RelErr;             /* relative error output parameter */

  log_workspace *log_general;
  log_workspace *log_profile;
  log_workspace *log_F2;
  log_workspace *log_Sq_Lcurve;
  log_workspace *log_Sq_Lcorner;
  log_workspace *log_EEJ;
  log_workspace *log_EEJ_Lcurve;
  log_workspace *log_EEJ_Lcorner;
  log_workspace *log_PDE;
  log_workspace *log_model;
  log_workspace *log_EEF;

  mag_sqfilt_workspace *sqfilt_workspace_p;
  mag_eej_workspace *eej_workspace_p;
  pde_workspace *pde_workspace_p;
  inverteef_workspace *inverteef_workspace_p;

  const char *output_file; /* output file */
  FILE *fp_output;

  green_workspace *green_workspace_p;
  kp_workspace *kp_workspace_p;
} mag_workspace;

/*
 * Prototypes
 */

mag_workspace *mag_alloc(const int year, const char *output_dir, const char *log_dir);
void mag_free(mag_workspace *w);
int mag_proc(const mag_params *params, track_workspace *track_p,
             satdata_mag *data, mag_workspace *w);

/* mag_log.c */
int mag_log_profile(const int header, const size_t ntrack,
                    const double kp, const int dir,
                    const mag_workspace *w);
int mag_log_F2(const int header, const mag_workspace *w);
int mag_log_Sq_Lcurve(const int header, const mag_workspace *w);
int mag_log_Sq_Lcorner(const int header, const mag_workspace *w);
int mag_log_EEJ(const int header, const mag_workspace *w);
int mag_log_EEJ_Lcurve(const int header, const mag_workspace *w);
int mag_log_EEJ_Lcorner(const int header, const mag_workspace *w);
int mag_log_PDE(const int header, const mag_workspace *w);
int mag_log_model(const int header, const mag_workspace *w);
int mag_log_EEF(const int header, const time_t t, const double phi,
                const double kp, const mag_workspace *w);

/* mag_sqfilt.c */
mag_sqfilt_workspace *mag_sqfilt_alloc(const size_t nmax_int, const size_t mmax_int,
                                       const size_t nmax_ext, const size_t mmax_ext);
void mag_sqfilt_free(mag_sqfilt_workspace *w);
int mag_sqfilt(mag_workspace *mag_p, mag_sqfilt_workspace *w);

/* mag_eej.c */
mag_eej_workspace *mag_eej_alloc(const int year, const size_t ncurr,
                                 const double qdlat_min,
                                 const double qdlat_max);
void mag_eej_free(mag_eej_workspace *w);
int mag_eej_proc(mag_track *track, double *J, mag_eej_workspace *w);

#endif /* INCLUDED_mag_h */
