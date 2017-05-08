/*
 * eph.h
 */

#ifndef INCLUDED_eph_h
#define INCLUDED_eph_h

#include <gsl/gsl_interp.h>

#include "hermite.h"

#ifndef INCLUDED_eph_data_h
#include "eph_data.h"
#define INCLUDED_eph_data_h
#endif

typedef struct
{
  const eph_data *data;
  hermite_workspace *hermite_x;
  hermite_workspace *hermite_y;
  hermite_workspace *hermite_z;
  gsl_interp_accel *acc;

  double *t;     /* ephemeris timestamps in seconds */
  size_t degree; /* degree of Hermite interpolation */
  size_t n;      /* number of ephemeris data points */
} eph_workspace;

/*
 * Prototypes
 */

eph_workspace *eph_alloc(const eph_data *data);
void eph_free(eph_workspace *w);
int eph_interp(const double t, double r_ECI[3], double v_ECI[3],
               eph_workspace *w);
int eph_interp_sph(const double t, double r_sph[3], eph_workspace *w);

#endif /* INCLUDED_eph_h */
