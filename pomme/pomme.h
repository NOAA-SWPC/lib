/*
 * pomme.h
 * Patrick Alken
 */

#ifndef INCLUDED_pomme_h
#define INCLUDED_pomme_h

#include <time.h>
#include <indices/indices.h>

#ifndef INCLUDED_pomme_mod_h
#include "pomme_mod.h"
#define INCLUDED_pomme_mod_h
#endif

#ifndef INCLUDED_POMME8_ext_h
#include "POMME8_ext.h"
#define INCLUDED_POMME8_ext_h
#endif

#ifndef INCLUDED_estist2_h
#include "estist2.h"
#define INCLUDED_estist2_h
#endif

#define R_EARTH_KM        (6371.2)

typedef struct
{
  double r_earth; /* Earth radius in km */

  tcontrol control;
  f107_workspace *f107_workspace_p;
  estist_workspace *estist_workspace_p;
  estist2_workspace *estist2_workspace_p;
  ace_workspace *ace_workspace_p;
} pomme_workspace;

typedef struct
{
  double r_earth;  /* Earth radius in km */
  const char *pomme_file_coeffs;
  const char *pomme_file_sm2geo;
  const char *pomme_file_gsm2geo;
  size_t ndeg;     /* SH degree of internal field */
  size_t flags;    /* POMME_xxx flags */
} pomme_parameters;

/* POMME flags */
#define POMME_NO_COF_READ         (1 << 0) /* don't read cof file */

/* spherical harmonic degree for main field */
#define POMME_MAIN_FIELD_DEG      (133)

/*
 * Prototypes
 */

extern int pomme_indx(int n, int m);

pomme_workspace *pomme_alloc(pomme_parameters *params);
pomme_workspace *pomme_alloc_default(void);
void pomme_free(pomme_workspace *w);
int pomme_set_deg(size_t ndeg, pomme_workspace *w);
int pomme_set_radius(double radius, pomme_workspace *w);
int pomme_calc(double theta, double phi, time_t t,
               double alt, double *B, pomme_workspace *w);
int pomme_calc_int(double theta, double phi, time_t t,
                   double alt, double *B, pomme_workspace *w);
int pomme_calc_ext(double theta, double phi, time_t t,
                   double alt, double *B, pomme_workspace *w);
int pomme_calc_ext_indices(double theta, double phi, time_t t, double alt,
                           double E_st, double I_st, double IMF_By, double Em,
                           double f107, double B[4], pomme_workspace *w);
int pomme_calc_sph(double theta, double phi, time_t t,
                   double alt, double B[4], pomme_workspace *w);
int pomme_call(double theta, double phi, time_t t, double alt,
               double E_st, double I_st, double IMF_By, double f107,
               double Em, double B[3], pomme_workspace *w);
int pomme_get_indices(time_t t, double *E_st, double *I_st, double *IMF_By,
                      double *Em, double *f107, pomme_workspace *w);

#endif /* INCLUDED_pomme_h */
