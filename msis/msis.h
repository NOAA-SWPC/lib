/*
 * msis.h
 * Patrick Alken
 */

#ifndef INCLUDED_msis_h
#define INCLUDED_msis_h

#include <indices/indices.h>

/* output array sizes */
#define MSIS_N_D         9
#define MSIS_N_T         2

typedef struct
{
  double n_He;  /* He density in m^{-3} */
  double n_O;   /* O density in m^{-3} */
  double n_N2;  /* N2 density in m^{-3} */
  double n_O2;  /* O2 density in m^{-3} */
  double n_Ar;  /* Ar density in m^{-3} */
  double n_H;   /* H density in m^{-3} */
  double n_N;   /* N density in m^{-3} */
  double n_n;   /* total neutral density in m^{-3} */
  double T_n;   /* neutral temperature (K) */
} msis_result;

typedef struct
{
  float D[MSIS_N_D]; /* D output array */
  float T[MSIS_N_T]; /* T output array */

  /* array containing IRI parameters for each altitude */
  msis_result *msis_results;

  size_t nalt; /* size of msis_results array */

  int day_of_year;
  double ap[7]; /* AP indices */

  double f107_override;
  double f107a_override;

  f107_workspace *f107_workspace_p;
  kp_workspace *kp_workspace_p;
} msis_workspace;

/*
 * Prototypes
 */

msis_workspace *msis_alloc(size_t nalt, const char *f107_datadir);
void msis_free(msis_workspace *w);
void msis_f107_override(double f107, double f107a, msis_workspace *w);
int msis_calc(double theta, double phi, time_t t, double altmin,
              double altstp, size_t nalt, msis_workspace *w);
msis_result *msis_get_result(size_t idx, msis_workspace *w);

/* MSIS prototypes */

void gtd7_(int *iyd, float *sec, float *alt, float *glat, float *glong,
           float *stl, float *f107a, float *f107, float *ap,
           int *mass, float *d, float *t);

#endif /* INCLUDED_msis_h */
