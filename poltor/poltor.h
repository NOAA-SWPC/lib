/*
 * poltor.h
 */

#ifndef INCLUDED_poltor_h
#define INCLUDED_poltor_h

#include <complex.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>

#include "bin2d.h"
#include "lls.h"
#include "magdata.h"
#include "track_weight.h"

/* use synthetic data for testing */
#define POLTOR_SYNTH_DATA            0

/* apply spatial weighting to data */
#define POLTOR_SPATIAL_WEIGHTS       1

/* mu_0 in units of: nT / (kA km^{-1}) */
#define POLTOR_MU_0                  (400.0 * M_PI)

/* types of coefficients for poltor_nmidx() */
#define POLTOR_IDX_PINT              (1 << 0)
#define POLTOR_IDX_PEXT              (1 << 1)
#define POLTOR_IDX_PSH               (1 << 2)
#define POLTOR_IDX_TOR               (1 << 3)

/* maximum number of possible current shells */
#define POLTOR_MAX_SHELL             2

/* maximum number of vector measurements (X,Y,Z) to fold into A^H A matrix at once */
#define POLTOR_BLOCK_SIZE            100000

#if !POLTOR_SYNTH_DATA
#define POLTOR_WEIGHT_X              (1.0)
#define POLTOR_WEIGHT_Y              (1.0)
#define POLTOR_WEIGHT_Z              (8.0)
#define POLTOR_WEIGHT_DX             (8.0)
#define POLTOR_WEIGHT_DY             (8.0)
#define POLTOR_WEIGHT_DZ             (8.0)
#else
#define POLTOR_WEIGHT_X              (1.0)
#define POLTOR_WEIGHT_Y              (1.0)
#define POLTOR_WEIGHT_Z              (1.0)
#define POLTOR_WEIGHT_DX             (1.0)
#define POLTOR_WEIGHT_DY             (1.0)
#define POLTOR_WEIGHT_DZ             (1.0)
#endif

#define POLTOR_FLG_QD_HARMONICS      (1 << 0) /* use QD coordinates for spherical harmonics */

typedef struct
{
  double R;          /* reference radius (km) */
  double b;          /* radius of internal current shell in km */
  double d;          /* radius of current shell inside satellite shell in km */
  double rmin;       /* bottom of toroidal current shell containing F-region currents (km) */
  double rmax;       /* top of toroidal current shell containing F-region currents (km) */;
  size_t nmax_int;   /* maximum internal spherical harmonic degree */
  size_t mmax_int;   /* maximum internal spherical harmonic order */
  size_t nmax_ext;   /* maximum external spherical harmonic degree */
  size_t mmax_ext;   /* maximum external spherical harmonic order */
  size_t nmax_sh;    /* maximum spherical harmonic degree for shell B_pol */
  size_t mmax_sh;    /* maximum spherical harmonic order for shell B_pol */
  size_t nmax_tor;   /* maximum spherical harmonic degree for B_tor */
  size_t mmax_tor;   /* maximum spherical harmonic order for B_tor */

  size_t shell_J;    /* order of Taylor series expansion of q_{nm}(r) for shell B_pol */

  double alpha_int;  /* damping parameter for internal coefficients */
  double alpha_sh;   /* damping parameter for poloidal shell coefficients */
  double alpha_tor;  /* damping parameter for toroidal shell coefficients */

  size_t flags;      /* POLTOR_FLG_xxx */

  magdata *data;     /* satellite data */
} poltor_parameters;

typedef struct
{
  double R;          /* reference radius (km) */
  double b;          /* radius of internal toroidal current shell (km) */
  double d;          /* radius of toroidal current shell inside satellite shell (km) */
  double rmin;       /* bottom of toroidal current shell containing F-region currents (km) */
  double rmax;       /* top of toroidal current shell containing F-region currents (km) */;
  double rmid;       /* midpoint of toroidal current shell (km) */
  size_t nmax_int;   /* maximum internal spherical harmonic degree */
  size_t mmax_int;   /* maximum internal spherical harmonic order */
  size_t nmax_ext;   /* maximum external spherical harmonic degree */
  size_t mmax_ext;   /* maximum external spherical harmonic order */
  size_t nmax_sh;    /* maximum external spherical harmonic degree for B_pol^{sh} */
  size_t mmax_sh;    /* maximum external spherical harmonic order for B_pol^{sh} */
  size_t nmax_tor;   /* maximum spherical harmonic degree for B_tor */
  size_t mmax_tor;   /* maximum spherical harmonic order for B_tor */
  size_t shell_J;    /* order of Taylor series expansion of q_{nm}(r) for shell B_pol */
  magdata *data;     /* satellite data */
  double alpha_int;  /* damping parameter for internal coefficients */
  double alpha_sh;   /* damping parameter for poloidal shell coefficients */
  double alpha_tor;  /* damping parameter for toroidal shell coefficients */

  size_t nmax_max;   /* max(nmax_int,nmax_ext,nmax_sh,nmax_tor) */
  size_t mmax_max;   /* max(mmax_int,mmax_ext,mmax_sh,mmax_tor) */

  poltor_parameters params;

  double *Pnm;       /* legendre array */
  double *dPnm;      /* legendre derivative array */
  double *Pnm2;      /* legendre array for along-track gradient */
  double *dPnm2;     /* legendre derivative array for along-track gradient */

  complex double *Ynm;     /* Pnm * exp(i m phi) */
  complex double *dYnm;    /* dPnm * exp(i m phi) */
  complex double *Ynm2;    /* Pnm * exp(i m phi2) */
  complex double *dYnm2;   /* dPnm * exp(i m phi2) */

  gsl_matrix_complex *A;   /* least squares matrix */
  gsl_vector_complex *rhs; /* rhs vector */
  gsl_vector *weights;     /* weights vector */
  gsl_vector *L;           /* Tikhonov regularization matrix = diag(L) */

  /*
   * coefficient vector is partitioned as
   * c = [ c_pint | c_pext | c_psh | c_tor ]
   * for B_pol^i, B_pol^e, B_pol^sh, B_tor
   *
   * c_psh is further partioned as
   * c_psh = [ c_psh^{(0)} | c_psh^{(1)} | ... | c_psh^{(J)} ]
   * where J = shell_J
   */
  gsl_vector_complex *c;

  gsl_vector_complex *residuals; /* residuals: rhs - A*c */
  double chisq;      /* residual sum of squares */
  int dof;           /* n - p degrees of freedom */

  size_t n;          /* total data in LS system */
  size_t p;          /* number of coefficients in LS system */
  size_t nblock;     /* max data to fold into normal matrix at a time */

  size_t nnm_sh;     /* number of (n,m) pairs for B_pol^{sh} */

  size_t p_pint;     /* number of coefficients for B_pol^i */
  size_t p_pext;     /* number of coefficients for B_pol^e */
  size_t p_psh;      /* number of coefficients for B_pol^{sh} */
  size_t p_tor;      /* number of coefficients for B_tor */

  size_t pint_offset; /* offset in 'c' of B_pol^i coefficients */
  size_t pext_offset; /* offset in 'c' of B_pol^e coefficients */
  size_t psh_offset;  /* offset in 'c' of B_pol^{sh} coefficients */
  size_t tor_offset;  /* offset in 'c' of B_tor coefficients */

  /*
   * Indexing for coefficients, base[n] gives the offset in the coef
   * vector for all the coefficients of degree n
   */
  size_t *base_int;
  size_t *base_ext;
  size_t *base_sh;
  size_t *base_tor;

  /* L-curve parameters */
  gsl_vector *reg_param;  /* regularization parameters */
  gsl_vector *rho;        /* residual norms */
  gsl_vector *eta;        /* solution norms */
  size_t nreg;            /* number of points on L-curve */

  size_t flags;       /* POLTOR_FLG_xxx */

  lls_complex_workspace *lls_workspace_p;
  gsl_integration_cquad_workspace *cquad_workspace_p;
} poltor_workspace;

/*
 * Prototypes
 */

poltor_workspace *poltor_alloc(const poltor_parameters *params);
void poltor_free(poltor_workspace *w);
int poltor_calc(poltor_workspace *w);
int poltor_solve(poltor_workspace *w);
int poltor_eval_J_shell(const double r, const double theta, const double phi,
                        double J[3], poltor_workspace *w);
int poltor_eval_K_tor(const double b, const double theta, const double phi,
                      double J[3], poltor_workspace *w);
int poltor_eval_K_ext(const double r, const double theta, const double phi,
                      double J[3], poltor_workspace *w);
int poltor_eval_chi_int(const double b, const double theta, const double phi,
                        double *chi, poltor_workspace *w);
int poltor_eval_chi_sh(const double r, const double theta, const double phi,
                       double *chi, poltor_workspace *w);
int poltor_eval_chi_ext(const double r, const double theta, const double phi,
                        double *chi, poltor_workspace *w);
int poltor_eval_B(const double r, const double theta, const double phi,
                  double B[4], poltor_workspace *w);
int poltor_eval_B_all(const double r, const double theta, const double phi,
                      double B_int[4], double B_ext[4], double B_sh[4], double B_tor[4],
                      poltor_workspace *w);
int poltor_write(const char *filename, poltor_workspace *w);
poltor_workspace *poltor_read(const char *filename);
gsl_complex poltor_get(const size_t cidx, const poltor_workspace *w);
double poltor_spectrum_int(const size_t n, poltor_workspace *w);
double poltor_spectrum_ext(const size_t n, poltor_workspace *w);
double poltor_spectrum_sh(const size_t n, poltor_workspace *w);
double poltor_spectrum_tor(const size_t n, poltor_workspace *w);
int poltor_print_spectrum(const char *filename, poltor_workspace *w);
size_t poltor_nmidx(const size_t type, const size_t n, const int m, poltor_workspace *w);
size_t poltor_jnmidx(const size_t j, const size_t n, const int m, poltor_workspace *w);
int poltor_regularize(gsl_vector *L, poltor_workspace *w);
int poltor_build_ls(const int fold, poltor_workspace *w);
double poltor_theta(const size_t idx, const poltor_workspace *w);
double poltor_theta_ns(const size_t idx, const poltor_workspace *w);

/* poltor_shell.c */
int poltor_shell_An(const size_t n, const size_t j, const double r,
                    double *result, poltor_workspace *w);
int poltor_shell_Bn(const size_t n, const size_t j, const double r,
                    double *result, poltor_workspace *w);

/* poltor_synth.c */
int poltor_synth(poltor_workspace *w);
int poltor_synth_print(poltor_workspace *w);

#endif /* INCLUDED_poltor_h */
