/*
 * pomme_wrap.c
 *
 * This module provides a C interface to the POMME model
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_math.h>

#include <indices/indices.h>

#include "common.h"
#include "estist_calc.h"
#include "pomme.h"

void pomme_geodetic(int control[MAXCONTROL], double fday,
                    double geod_alt, double geod_lat, double lon,
                    double est, double ist, double imf_by,
                    double f107, double Em, double *geod_x,
                    double *geod_y, double *geod_z, double *h,
                    double *f, double *incl, double *decl);
void pomme_geocentric(int control[MAXCONTROL], double fday,
                      double r, double lat, double lon,
                      double est, double ist, double imf_by,
                      double f107, double Em, double *x,
                      double *y, double *z);
void pomme_internal_geocentric(int limit_deg, double fyear, double rad,
                               double theta, double phi, double *mx,
                               double *my, double *mz, double *cx,
                               double *cy, double *cz);
void pomme_external_geocentric(int control[MAXCONTROL], double fday,
                               double rrel, double theta, double phi,
                               double est, double ist, double imf_by,
                               double f107, double Em, double *x,
                               double *y, double *z);

/*
pomme_alloc()
  Allocate pomme workspace

Inputs: params - pomme parameters
*/

pomme_workspace *
pomme_alloc(pomme_parameters *params)
{
  pomme_workspace *w;

  if (params->ndeg > MAX_0_DEG)
    {
      fprintf(stderr, "pomme_alloc: error: ndeg exceeds MAX_0_DEG\n");
      return 0;
    }

  w = calloc(1, sizeof(pomme_workspace));
  if (!w)
    {
      fprintf(stderr, "pomme_alloc: malloc failed: %s\n", strerror(errno));
      return 0;
    }

  w->f107_workspace_p = f107_alloc(params->f107_file);

  if (params->ace_file)
    w->ace_workspace_p = ace_alloc(params->ace_file);

  if (params->dst_file != NULL)
    w->estist_calc_workspace_p = estist_calc_alloc(params->dst_file);
  else
    w->estist_workspace_p = estist_alloc(ESTIST_IDX_FILE);

  w->control[POS_0] = params->ndeg; /* internal field coefficients */
  w->control[POS_1] = 15;
  w->control[POS_SM] = DEG_SM;
  w->control[POS_SM_IND] = DEG_SM;
  w->control[POS_F107] = DEG_F107;
  w->control[POS_F107_IND] = DEG_F107;
  w->control[POS_EST] = 1;
  w->control[POS_IST] = 1;
  w->control[POS_EST_FAC] = 1;
  w->control[POS_GSM] = DEG_GSM;
  w->control[POS_GSM_IND] = DEG_GSM;
  w->control[POS_IMFBY] = 1;
  w->control[POS_IMFBY_IND] = 1;
  w->control[POS_EM] = DEG_EM;
  w->control[POS_EM_IND] = DEG_EM;

  w->r_earth = params->r_earth;

  return w;
} /* pomme_alloc() */

/* default settings */
pomme_workspace *
pomme_alloc_default(void)
{
  pomme_workspace *w;
  pomme_parameters params;

  params.r_earth = R_EARTH_KM;
  params.ndeg = POMME_MAIN_FIELD_DEG;
  params.dst_file = NULL;
  params.f107_file = F107_IDX_FILE;
  params.ace_file = ACE_IDX_FILE;

  w = pomme_alloc(&params);

  return w;
} /* pomme_alloc_default() */

int
pomme_set_deg(size_t ndeg, pomme_workspace *w)
{
  w->control[POS_0] = ndeg;

  return 0;
}

/*
pomme_set_radius()
  Set Earth reference radius (in km)
*/

int
pomme_set_radius(double radius, pomme_workspace *w)
{
  w->r_earth = radius;
  return 0;
}

void
pomme_free(pomme_workspace *w)
{
  if (w->f107_workspace_p)
    f107_free(w->f107_workspace_p);

  if (w->estist_workspace_p)
    estist_free(w->estist_workspace_p);

  if (w->estist_calc_workspace_p)
    estist_calc_free(w->estist_calc_workspace_p);

  if (w->ace_workspace_p)
    ace_free(w->ace_workspace_p);

  free(w);
} /* pomme_free() */

/*
pomme_calc()
  Call the POMME model

Inputs: theta - geographic colatitude (radians)
        phi   - geographic longitude (radians)
        t     - timestamp (UTC)
        alt   - altitude above sea level in km
        B     - (output) where to store magnetic field
        w     - workspace

Return: success (0) or error

Notes: B is a 4 element array; on output:

B[0] = X component of field (northward, -B_theta)
B[1] = Y component of field (eastward, B_phi)
B[2] = Z component of field (vertically down, -B_r)
B[3] = F component (field intensity)

Magnetic field units in output are Tesla
*/

int
pomme_calc(double theta, double phi, time_t t,
           double alt, double *B, pomme_workspace *w)
{
  int s = 0;
  double E_st = 0.0, I_st = 0.0;
  double IMF_By = 0.0, Em = 0.5;
  double f107 = 120.0;

  s += pomme_get_indices(1, t, &E_st, &I_st, &IMF_By, &Em, &f107, w);

  s += pomme_call(theta, phi, t, alt, E_st, I_st, IMF_By, f107,
                  Em, B, w);

  return s;
} /* pomme_calc() */

/*
pomme_calc_int()
  Call the POMME model for given parameters for the internal field

Inputs: theta  - geographic colatitude (radians)
        phi    - geographic longitude (radians)
        t      - timestamp
        alt    - altitude in km
        B      - (output) where to store B field
        w      - workspace

Notes: B is a 4 element array; on output:

B[0] = X component of field
B[1] = Y component of field
B[2] = Z component of field
B[3] = F component (field intensity)

Output units are Tesla
*/

int
pomme_calc_int(double theta, double phi, time_t t, double alt,
               double B[4], pomme_workspace *w)
{
  int s = 0;
  double fday = time2fday(t);
  double fyear = 2000.0 + fday / 365.25;
  double r = w->r_earth + alt;
  double mx, my, mz;
  double cx, cy, cz;

  pomme_internal_geocentric(w->control[POS_0], fyear, r, theta, phi,
                            &mx, &my, &mz, &cx, &cy, &cz);

  B[0] = (mx + cx) * 1.0e-9;
  B[1] = (my + cy) * 1.0e-9;
  B[2] = (mz + cz) * 1.0e-9;

  B[3] = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

  return s;
} /* pomme_calc_int() */

/*
pomme_calc_ext()
  Call the POMME model for given parameters for the external field

Inputs: theta  - geographic colatitude (radians)
        phi    - geographic longitude (radians)
        t      - timestamp
        alt    - altitude in km
        B      - (output) where to store B field
        w      - workspace

Notes: B is a 4 element array; on output:

B[0] = X component of field
B[1] = Y component of field
B[2] = Z component of field
B[3] = F component (field intensity)

Output units are Tesla
*/

int
pomme_calc_ext(double theta, double phi, time_t t, double alt,
               double B[4], pomme_workspace *w)
{
  int s = 0;
  double E_st = 0.0, I_st = 0.0;
  double IMF_By = 0.0, Em = 0.5;
  double f107 = 120.0;

  s += pomme_get_indices(1, t, &E_st, &I_st, &IMF_By, &Em, &f107, w);

  s += pomme_calc_ext_indices(theta, phi, t, alt, E_st, I_st,
                              IMF_By, Em, f107, B, w);

  return s;
} /* pomme_calc_ext() */

/*
pomme_calc_ext2()
  Call the POMME model for given parameters for the external field,
except E_st and I_st are passed as parameters and not computed
internally. Since these can be time consuming to compute from D_st,
the user has the option to pass these in to avoid the computation.

Inputs: theta  - geographic colatitude (radians)
        phi    - geographic longitude (radians)
        t      - timestamp
        alt    - altitude in km
        E_st   - E_st index (nT)
        I_st   - I_st index (nT)
        B      - (output) where to store B field
        w      - workspace

Notes: B is a 4 element array; on output:

B[0] = X component of field
B[1] = Y component of field
B[2] = Z component of field
B[3] = F component (field intensity)

Output units are nT
*/

int
pomme_calc_ext2(const double theta, const double phi, const time_t t, const double alt,
                const double E_st, const double I_st, double B[4], pomme_workspace *w)
{
  int s = 0;
  double dummy;
  double IMF_By = 0.0, Em = 0.5;
  double f107 = 120.0;
  size_t i;

  s += pomme_get_indices(0, t, &dummy, &dummy, &IMF_By, &Em, &f107, w);

  s += pomme_calc_ext_indices(theta, phi, t, alt, E_st, I_st,
                              IMF_By, Em, f107, B, w);

  for (i = 0; i < 4; ++i)
    B[i] *= 1.0e9;

  return s;
}

/*
pomme_calc_ext_indices()
  Call the POMME model for given parameters for the external field

Inputs: theta  - geographic colatitude (radians)
        phi    - geographic longitude (radians)
        t      - timestamp
        alt    - altitude in km
        E_st   - Est index
        I_st   - Ist index
        IMF_By - IMF By
        Em     - merging electric field
        f107   - F10.7 20 months prior
        B      - (output) where to store B field
        w      - workspace

Notes: B is a 4 element array; on output:

B[0] = X component of field
B[1] = Y component of field
B[2] = Z component of field
B[3] = F component (field intensity)

Output units are Tesla
*/

int
pomme_calc_ext_indices(double theta, double phi, time_t t, double alt,
                       double E_st, double I_st, double IMF_By, double Em,
                       double f107, double B[4], pomme_workspace *w)
{
  int s = 0;
  double fday = time2fday(t);
  double rrel = (w->r_earth + alt) / w->r_earth;
  double x, y, z;

  pomme_external_geocentric(w->control, fday, rrel, theta, phi,
                            E_st, I_st, IMF_By, f107, Em,
                            &x, &y, &z);

  /* convert to T */
  B[0] = x * 1.0e-9;
  B[1] = y * 1.0e-9;
  B[2] = z * 1.0e-9;

  B[3] = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

  return s;
} /* pomme_calc_ext_indices() */

/*
pomme_calc_sph()
  Call the POMME model for the magnetic field and return results in
spherical coordinates

Inputs: theta - geographic colatitude (radians)
        phi   - geographic longitude (radians)
        t     - timestamp (UTC)
        alt   - altitude above sea level in km
        B     - (output) where to store magnetic field
        w     - workspace

Return: success (0) or error

Notes: B is a 4 element array; on output:

B[0] = r component of field
B[1] = theta component of field
B[2] = phi component of field
B[3] = F component (field intensity)

Magnetic field units in output are Tesla
*/

int
pomme_calc_sph(double theta, double phi, time_t t,
               double alt, double B[4], pomme_workspace *w)
{
  int s = 0;
  double E_st = 0.0, I_st = 0.0;
  double IMF_By = 0.0, Em = 0.5;
  double f107 = 120.0;
  double V[4];

  s += pomme_get_indices(1, t, &E_st, &I_st, &IMF_By, &Em, &f107, w);

  s += pomme_call(theta, phi, t, alt, E_st, I_st, IMF_By, f107,
                  Em, V, w);

  B[0] = -V[2];
  B[1] = -V[0];
  B[2] = V[1];
  B[3] = V[3];

  return s;
} /* pomme_calc_sph() */

/*
pomme_call()
  Call the POMME model for given parameters

Inputs: theta  - geographic colatitude (radians)
        phi    - geographic longitude (radians)
        t      - timestamp
        alt    - altitude in km
        E_st   - Est index
        I_st   - Ist index
        IMF_By - IMF By value
        f107   - F10.7 index
        Em     - merging electric field (mV/m)
        B      - (output) where to store B field
        w      - workspace

Notes: B is a 4 element array; on output:

B[0] = X component of field
B[1] = Y component of field
B[2] = Z component of field
B[3] = F component (field intensity)

Output units are Tesla
*/

int
pomme_call(double theta, double phi, time_t t, double alt,
           double E_st, double I_st, double IMF_By, double f107,
           double Em, double B[4], pomme_workspace *w)
{
  int s = 0;
  double fday = time2fday(t);
  double x, y, z;
  double lat_deg = 90.0 - theta * 180.0 / M_PI;
  double lon_deg = phi * 180.0 / M_PI;
  double r = w->r_earth + alt;

  pomme_geocentric(w->control, fday, r, lat_deg, lon_deg,
                   E_st, I_st, IMF_By, f107, Em,
                   &x, &y, &z);

  B[0] = x * 1.0e-9;
  B[1] = y * 1.0e-9;
  B[2] = z * 1.0e-9;

  B[3] = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

  return s;
} /* pomme_call() */

/*
pomme_get_indices()
  Look up indices needed for POMME external field model

Inputs: calc_estist - compute/retrieve E_st/Ist indices?
        t           - timestamp
        E_st        - (output) E_st index
        I_st        - (output) I_st index
        IMF_By      - (output) IMF B_y
        Em          - (output) merging electric field Em
        f107        - (output) F10.7 index
        w           - workspace

Notes:
1) If F10.7, IMF or Em data is not available, default values are used
2) If Est/Ist is not available, an error is returned, since this can
dramatically affect the POMME external field calculation
*/

int
pomme_get_indices(const int calc_estist, time_t t, double *E_st, double *I_st, double *IMF_By,
                  double *Em, double *f107, pomme_workspace *w)
{
  int s = 0;
  int status = 0;
  double IMF_B[4]; /* nT */
  double V;  /* solar wind velocity (km/s) */
  double BT; /* component of B perpendicular to Bx */
  double beta;

  /* calculate F10.7A 20 months prior */
  s = f107a_get(t - 20.0*365.25/12.0*86400.0, f107, w->f107_workspace_p);
  if (s)
    *f107 = 120.0; /* default value */

  if (calc_estist)
    {
      if (w->estist_calc_workspace_p)
        s = estist_calc_get(t, E_st, I_st, w->estist_calc_workspace_p);
      else
        s = estist_get(t, E_st, I_st, w->estist_workspace_p);

    if (s)
      {
        /* raise error since these are critical parameters */
        *E_st = 0.0;
        *I_st = 0.0;
        status += s;
      }
    }

  if (w->ace_workspace_p)
    {
      s = ace_get(t, &IMF_B[0], &IMF_B[1], &IMF_B[2], &V, w->ace_workspace_p);
      if (s)
        {
          *Em = 0.5;
          *IMF_By = 0.0;
        }
      else
        {
          /* compute sqrt(By^2 + Bz^2) */
          BT = sqrt(IMF_B[1]*IMF_B[1] + IMF_B[2]*IMF_B[2]);

          if (BT == 0.0)
            beta = 1.0;
          else
            beta = acos(IMF_B[2] / BT);

          *Em = V * BT * sin(beta / 2.0) * sin(beta / 2.0) / 1000.0;
          *IMF_By = IMF_B[1];
        }
    }
  else
    {
      *Em = 0.5;
      *IMF_By = 0.0;
    }

  return status;
} /* pomme_get_indices() */
