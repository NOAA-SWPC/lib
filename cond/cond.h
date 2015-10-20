/*
 * cond.h
 * Patrick Alken
 */

#ifndef INCLUDED_cond_h
#define INCLUDED_cond_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>

#include "iri.h"
#include "msis.h"
#include "msynth.h"

typedef struct
{
  double n_n;     /* neutral density in cm^{-3} */
  double n_e;     /* electron density in cm^{-3} */
  double v_en;    /* eletron-neutral collision frequency in 1/s */
  double v_ei;    /* electron-ion collision frequency in 1/s */
  double v_e;     /* v_en + v_ei */
  double v_e_sq;  /* v_e^2 */
  double T_e;     /* electron temperature from IRI (K) */
  double T_i;     /* ion temperature from IRI (K) */
  double T_n_iri; /* neutral temperature from IRI (K) */
  double T_n_msis; /* neutral temperature from MSIS (K) */
  double w_e;     /* electron gyro-frequency */
  double w_e_sq;  /* w_e^2 */

  double n_Op;    /* O+ density in cm^{-3} */
  double n_Hp;    /* H+ density in cm^{-3} */
  double n_HEp;   /* HE+ density in cm^{-3} */
  double n_O2p;   /* O2+ density in cm^{-3} */
  double n_NOp;   /* NO+ density in cm^{-3} */
  double n_Np;    /* N+ density in cm^{-3} */
  double n_N2p;   /* N2+ density in cm^{-3} */

  double n_He;    /* He density in cm^{-3} */
  double n_O;     /* O density in cm^{-3} */
  double n_O2;    /* O2 density in cm^{-3} */
  double n_N2;    /* N2 density in cm^{-3} */
  double n_Ar;    /* Ar density in cm^{-3} */
  double n_H;     /* H density in cm^{-3} */
  double n_N;     /* N density in cm^{-3} */

  double v_Op;    /* O+ frequency in 1/s */
  double v_Hp;    /* H+ frequency in 1/s */
  double v_HEp;   /* HE+ frequency in 1/s */
  double v_O2p;   /* O2+ frequency in 1/s */
  double v_NOp;   /* NO+ frequency in 1/s */
  double v_Np;    /* N+ frequency in 1/s */

  double v_O;     /* O frequency in 1/s */
  double v_O2;    /* O2 frequency in 1/s */
  double v_N2;    /* N2 frequency in 1/s */

  double v_Op_sq;  /* v_Op^2 */
  double v_Hp_sq;  /* v_Hp^2 */
  double v_HEp_sq; /* v_HEp^2 */
  double v_O2p_sq; /* v_O2p^2 */
  double v_NOp_sq; /* v_NOp^2 */
  double v_Np_sq;  /* v_Np^2 */

  double v_O_sq;   /* v_O^2 */
  double v_O2_sq;  /* v_O2^2 */
  double v_N2_sq;  /* v_N2^2 */

  /* gyro-frequencies */
  double w_Op;     /* w_Op */
  double w_Hp;     /* w_Hp */
  double w_HEp;    /* w_HEp */
  double w_O2p;    /* w_O2p */
  double w_NOp;    /* w_NOp */
  double w_Np;     /* w_Np */

  double w_O;      /* w_O */
  double w_O2;     /* w_O2 */
  double w_N2;     /* w_N2 */

  double w_Op_sq;  /* w_Op^2 */
  double w_Hp_sq;  /* w_Hp^2 */
  double w_HEp_sq; /* w_HEp^2 */
  double w_O2p_sq; /* w_O2p^2 */
  double w_NOp_sq; /* w_NOp^2 */
  double w_Np_sq;  /* w_Np^2 */

  double w_O_sq;   /* w_O^2 */
  double w_O2_sq;  /* w_O2^2 */
  double w_N2_sq;  /* w_N2^2 */

  double B;        /* magnetic field in T */
  double height;   /* altitude in km */

  /* conductivity outputs */
  double sigma_0;
  double sigma_p;
  double sigma_h;
} cond_result;

typedef struct
{
  msynth_workspace *msynth_workspace_p;
  iri_workspace *iri_workspace_p;
  msis_workspace *msis_workspace_p;

  double e_sq;  /* e^2 */

  /*
   * "factor-4": factor with which to multiply the electron collision
   * frequency in conductivity computations
   */
  double alpha;

  /* scaling factors for error analysis in IRI/MSIS */
  double iri_n_scale;
  double iri_T_scale;

  cond_result *cond_results; /* results for each altitude step */
  size_t nalt;               /* number of altitude steps */
} cond_workspace;

#define COND_MASS_O   (16.0 * GSL_CONST_MKSA_MASS_PROTON)
#define COND_MASS_H   (GSL_CONST_MKSA_MASS_PROTON)
#define COND_MASS_HE  (4.0 * GSL_CONST_MKSA_MASS_PROTON)
#define COND_MASS_O2  (32.0 * GSL_CONST_MKSA_MASS_PROTON)
#define COND_MASS_NO  (30.01 * GSL_CONST_MKSA_MASS_PROTON)
#define COND_MASS_N   (14.01 * GSL_CONST_MKSA_MASS_PROTON)
#define COND_MASS_N2  (28.02 * GSL_CONST_MKSA_MASS_PROTON)
#define COND_MASS_AR  (39.95 * GSL_CONST_MKSA_MASS_PROTON)

/*
 * Prototypes
 */

cond_workspace *cond_alloc(size_t nalt, const char *f107_datafile);
void cond_free(cond_workspace *w);
void cond_f107_override(double f107, double f107a, cond_workspace *w);
void cond_set_alpha(double alpha, cond_workspace *w);
int cond_set_error_scale(double iri_n_scale, double iri_T_scale,
                         cond_workspace *w);
int cond_calc(double theta, double phi, time_t t,
              double altmin, double altstp, size_t nalt,
              cond_workspace *w);
cond_result *cond_get_result(size_t idx, cond_workspace *w);

#endif /* INCLUDED_cond_h */
