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

#include "estist_calc.h"

#define R_EARTH_KM        (6371.2)

typedef struct
{
  double r_earth; /* Earth radius in km */

  tcontrol control;
  f107_workspace *f107_workspace_p;
  estist_workspace *estist_workspace_p;
  estist_calc_workspace *estist_calc_workspace_p;
  ace_workspace *ace_workspace_p;
} pomme_workspace;

typedef struct
{
  double r_earth;  /* Earth radius in km */
  const char *coef_file;
  const char *sm2geo_file;
  const char *gsm2geo_file;
  const char *dst_file;
  const char *f107_file;
  const char *ace_file;
  size_t ndeg;     /* SH degree of internal field */
} pomme_parameters;

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
int pomme_calc_ext2(const double theta, const double phi, const time_t t, const double alt,
                    const double E_st, const double I_st, double B[4], pomme_workspace *w);
int pomme_calc_ext_indices(double theta, double phi, time_t t, double alt,
                           double E_st, double I_st, double IMF_By, double Em,
                           double f107, double B[4], pomme_workspace *w);
int pomme_calc_sph(double theta, double phi, time_t t,
                   double alt, double B[4], pomme_workspace *w);
int pomme_call(double theta, double phi, time_t t, double alt,
               double E_st, double I_st, double IMF_By, double f107,
               double Em, double B[3], pomme_workspace *w);
int pomme_get_indices(const int calc_estist, time_t t, double *E_st, double *I_st, double *IMF_By,
                      double *Em, double *f107, pomme_workspace *w);

#endif /* INCLUDED_pomme_h */
