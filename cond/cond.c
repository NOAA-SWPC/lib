/*
 * cond.c
 *
 * This file provides a C interface to the IRI-2007 and MSIS
 * models (written in fortran).
 *
 * 28 Feb 2009: fixed cond_direct and cond_pederson to use
 *              IRI ion densities. cond_hall does not give
 *              correct output with IRI ion densities.
 *
 * 24 March 2012: fixed cond_hall to use IRI ion densities
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>

#include <common/common.h>
#include <msynth/msynth.h>

#include "cond.h"
#include "iri.h"
#include "msis.h"

static int cond_error_scaling(iri_result *iri, msis_result *msis,
                              cond_workspace *w);
static double cond_direct(cond_result *result, cond_workspace *w);
static double cond_pederson(cond_result *result, cond_workspace *w);
static double cond_hall(cond_result *result, cond_workspace *w);

cond_workspace *
cond_alloc(size_t nalt, const char *f107_datafile)
{
  cond_workspace *w;

  w = (cond_workspace *) calloc(1, sizeof(cond_workspace));
  if (!w)
    {
      fprintf(stderr, "cond_alloc: calloc failed: %s\n",
              strerror(errno));
      return 0;
    }

  w->cond_results = malloc(nalt * sizeof(cond_result));
  if (!w->cond_results)
    {
      fprintf(stderr, "cond_alloc: malloc failed: %s\n", strerror(errno));
      cond_free(w);
      return 0;
    }

  w->msynth_workspace_p = msynth_igrf_read(MSYNTH_IGRF_FILE);
  w->iri_workspace_p = iri_alloc(nalt, f107_datafile);
  w->msis_workspace_p = msis_alloc(nalt, f107_datafile);

  w->e_sq = pow(GSL_CONST_MKSA_ELECTRON_CHARGE, 2.0);

  /* initialize factor-4 correction */
  w->alpha = 4.0;

  w->iri_n_scale = 1.0;
  w->iri_T_scale = 1.0;

  return (w);
} /* cond_alloc() */

void
cond_free(cond_workspace *w)
{
  if (w->cond_results)
    free(w->cond_results);

  if (w->msynth_workspace_p)
    msynth_free(w->msynth_workspace_p);

  if (w->iri_workspace_p)
    iri_free(w->iri_workspace_p);

  if (w->msis_workspace_p)
    msis_free(w->msis_workspace_p);

  free(w);
} /* cond_free() */

void
cond_f107_override(double f107, double f107a, cond_workspace *w)
{
  iri_f107_override(f107, w->iri_workspace_p);
  msis_f107_override(f107, f107a, w->msis_workspace_p);
} /* cond_f107_override() */

/*
cond_set_alpha()
  Set factor which will multiply the electron collision frequency
to account for the gradient drift instability
*/

void
cond_set_alpha(double alpha, cond_workspace *w)
{
  w->alpha = alpha;
} /* cond_set_alpha() */

/*
cond_set_error_scale
  Set scaling factors for error analysis in IRI/MSIS

Inputs: iri_n_scale - scale IRI densities (both electrons and ions) by
                      this factor prior to computing conductivities
        iri_T_scale - scale IRI temperatures by this factor prior to
                      computing conductivities
        w           - workspace
*/

int
cond_set_error_scale(double iri_n_scale, double iri_T_scale,
                     cond_workspace *w)
{
  int s = 0;

  w->iri_n_scale = iri_n_scale;
  w->iri_T_scale = iri_T_scale;

  return s;
} /* cond_set_error_scale() */

/*
cond_error_scaling()
  Perform scaling if requested on IRI/MSIS parameters for error analysis
*/

static int
cond_error_scaling(iri_result *iri, msis_result *msis, cond_workspace *w)
{
  int s = 0;

  iri->Ne *= w->iri_n_scale;
  iri->n_Op *= w->iri_n_scale;
  iri->n_Hp *= w->iri_n_scale;
  iri->n_HEp *= w->iri_n_scale;
  iri->n_O2p *= w->iri_n_scale;
  iri->n_NOp *= w->iri_n_scale;
  iri->n_Np *= w->iri_n_scale;
  iri->n_clus *= w->iri_n_scale;

  iri->Tn *= w->iri_T_scale;
  iri->Ti *= w->iri_T_scale;
  iri->Te *= w->iri_T_scale;

  return s;
} /* cond_error_scaling() */

/*
cond_calc()
  Call the IRI and MSIS models and compute relevant conductivity
parameters.

Inputs: theta             - geographic colatitude (radians)
        phi               - geographic longitude (radians)
        t                 - timestamp
        altmin            - minimum altitude in km
        altstp            - altitude step size in km
        nalt              - number of altitude steps
        w                 - workspace

Return: success or error

Notes:

The units of the sigma outputs are SI (A / V / m)
*/

int
cond_calc(double theta, double phi, time_t t, double altmin,
          double altstp, size_t nalt, cond_workspace *w)
{
  int s = 0;
  const double R = R_EARTH_KM;
  const double tyear = get_year(t);
  size_t i, j;
  cond_result result;

  /* compute IRI and MSIS values */
  s += iri_calc(theta, phi, t, altmin, altstp, nalt, w->iri_workspace_p);
  s += msis_calc(theta, phi, t, altmin, altstp, nalt, w->msis_workspace_p);

  if (s)
    {
      fprintf(stderr, "cond_calc: error computing IRI/MSIS quantities\n");
      return s;
    }

  for (i = 0; i < nalt; ++i)
    {
      double alt = altmin + i * altstp;
      double r = R + alt;
      iri_result *result_iri = iri_get_result(i, w->iri_workspace_p);
      msis_result *result_msis = msis_get_result(i, w->msis_workspace_p);
      double A;     /* mean molecular weight */
      double B[4];  /* magnetic field */

      /* perform scaling for error analysis if necessary */
      cond_error_scaling(result_iri, result_msis, w);

      /*result.n_n = result_msis->n_n;*/
      result.n_n = result_msis->n_O +
                   result_msis->n_O2 +
                   result_msis->n_N2;

      A = (16.0*result_msis->n_O + 32.0*result_msis->n_O2 +
           28.02*result_msis->n_N2) / result.n_n;

      /*
       * convert densities to cm^{-3} since Kelley's formulas use that
       * assumption
       */

      result.n_n *= 1.0e-6;
      result.n_e = result_iri->Ne * 1.0e-6;

      result.n_Op = result_iri->n_Op * 1.0e-6;
      result.n_Hp = result_iri->n_Hp * 1.0e-6;
      result.n_HEp = result_iri->n_HEp * 1.0e-6;
      result.n_O2p = result_iri->n_O2p * 1.0e-6;
      result.n_NOp = result_iri->n_NOp * 1.0e-6;
      result.n_Np = result_iri->n_Np * 1.0e-6;

      /*
       * sometimes IRI calculates densities where:
       * Sum [ n_i ] != n_e
       * but the conductivity formulas require this condition to work
       * (particularly the Hall conductivity)
       */
      result.n_e = result.n_Op + result.n_Hp + result.n_HEp +
                   result.n_O2p + result.n_NOp + result.n_Np;

      result.n_He = result_msis->n_He * 1.0e-6;
      result.n_O = result_msis->n_O * 1.0e-6;
      result.n_O2 = result_msis->n_O2 * 1.0e-6;
      result.n_N2 = result_msis->n_N2 * 1.0e-6;
      result.n_Ar = result_msis->n_Ar * 1.0e-6;
      result.n_H = result_msis->n_H * 1.0e-6;
      result.n_N = result_msis->n_N * 1.0e-6;

      /* store temperatures */

      result.T_e = result_iri->Te;
      result.T_i = result_iri->Ti;
      result.T_n_iri = result_iri->Tn;
      result.T_n_msis = result_msis->T_n;

      /*
       * Below 120km, IRI does not provide temperatures, but assumes
       * thermal equilibrium: T_e = T_i = T_n, so use the neutral
       * temperature from MSIS for this region
       */
      if (result.T_i < 0.0)
        result.T_i = result_msis->T_n;
      if (result.T_e < 0.0)
        result.T_e = result_msis->T_n;
      if (result.T_n_iri < 0.0)
        result.T_n_iri = result_msis->T_n;

      result.v_en = 5.4e-10 * result.n_n * sqrt(result.T_e);

      /*
       * For altitudes below 80km, IRI does not produce an electron
       * density since it is so small. In this region, set v_ei = 0
       */
      if (result.n_e < 1.0e-5)
        result.v_ei = 0.0;
      else
        result.v_ei = (34.0 + 4.18*log10(result.T_e * result.T_e * result.T_e / result.n_e)) *
                     result.n_e * pow(result.T_e, -1.5);

      result.v_e = result.v_en + result.v_ei;

      /* effective collision frequency (factor-4 correction) */
      result.v_e *= w->alpha;

      result.v_e_sq = result.v_e * result.v_e;

      /* N2+ is not given by IRI so estimate from MSIS */
      result.n_N2p = result.n_N2 / result.n_n * result.n_e;

      result.v_Op = 2.6e-9 * (result.n_Op + result.n_n) * sqrt(1.0 / A);
      result.v_Op_sq = result.v_Op * result.v_Op;
      result.v_Hp = 2.6e-9 * (result.n_Hp + result.n_n) * sqrt(1.0 / A);
      result.v_Hp_sq = result.v_Hp * result.v_Hp;
      result.v_HEp = 2.6e-9 * (result.n_HEp + result.n_n) * sqrt(1.0 / A);
      result.v_HEp_sq = result.v_HEp * result.v_HEp;
      result.v_O2p = 2.6e-9 * (result.n_O2p + result.n_n) * sqrt(1.0 / A);
      result.v_O2p_sq = result.v_O2p * result.v_O2p;
      result.v_NOp = 2.6e-9 * (result.n_NOp + result.n_n) * sqrt(1.0 / A);
      result.v_NOp_sq = result.v_NOp * result.v_NOp;
      result.v_Np = 2.6e-9 * (result.n_Np + result.n_n) * sqrt(1.0 / A);
      result.v_Np_sq = result.v_Np * result.v_Np;

      result.v_O = 2.6e-9 * (result.n_O + result.n_n) * sqrt(1.0 / A);
      result.v_O_sq = result.v_O * result.v_O;
      result.v_O2 = 2.6e-9 * (result.n_O2 + result.n_n) * sqrt(1.0 / A);
      result.v_O2_sq = result.v_O2 * result.v_O2;
      result.v_N2 = 2.6e-9 * (result.n_N2 + result.n_n) * sqrt(1.0 / A);
      result.v_N2_sq = result.v_N2 * result.v_N2;

      msynth_eval(tyear, r, theta, phi, B, w->msynth_workspace_p);
      
      /* convert to T */
      for (j = 0; j < 4; ++j)
        B[j] *= 1.0e-9;

      /* compute gyro-frequencies */

      result.w_e = GSL_CONST_MKSA_ELECTRON_CHARGE * B[3] /
                   GSL_CONST_MKSA_MASS_ELECTRON;
      result.w_e_sq = result.w_e * result.w_e;

      result.w_Op = GSL_CONST_MKSA_ELECTRON_CHARGE * B[3] / COND_MASS_O;
      result.w_Op_sq = result.w_Op * result.w_Op;
      result.w_Hp = GSL_CONST_MKSA_ELECTRON_CHARGE * B[3] / COND_MASS_H;
      result.w_Hp_sq = result.w_Hp * result.w_Hp;
      result.w_HEp = GSL_CONST_MKSA_ELECTRON_CHARGE * B[3] / COND_MASS_HE;
      result.w_HEp_sq = result.w_HEp * result.w_HEp;
      result.w_O2p = GSL_CONST_MKSA_ELECTRON_CHARGE * B[3] / COND_MASS_O2;
      result.w_O2p_sq = result.w_O2p * result.w_O2p;
      result.w_NOp = GSL_CONST_MKSA_ELECTRON_CHARGE * B[3] / COND_MASS_NO;
      result.w_NOp_sq = result.w_NOp * result.w_NOp;
      result.w_Np = GSL_CONST_MKSA_ELECTRON_CHARGE * B[3] / COND_MASS_N;
      result.w_Np_sq = result.w_Np * result.w_Np;

      result.w_O = GSL_CONST_MKSA_ELECTRON_CHARGE * B[3] / COND_MASS_O;
      result.w_O_sq = result.w_O * result.w_O;
      result.w_O2 = GSL_CONST_MKSA_ELECTRON_CHARGE * B[3] / COND_MASS_O2;
      result.w_O2_sq = result.w_O2 * result.w_O2;
      result.w_N2 = GSL_CONST_MKSA_ELECTRON_CHARGE * B[3] / COND_MASS_N2;
      result.w_N2_sq = result.w_N2 * result.w_N2;

      result.B = B[3];
      result.height = alt;

      /* compute and store conductivities */
      result.sigma_0 = cond_direct(&result, w);
      result.sigma_p = cond_pederson(&result, w);
      result.sigma_h = cond_hall(&result, w);

      w->cond_results[i] = result;

#if 0
      printf("%f %e %e %e %e %e %e %e %e %e %e %e\n",
             alt,
             result.n_n,
             result.n_e,
             result.T_e,
             result.T_i,
             result.v_e,
             result.v_Op,
             result.v_O2p,
             result.v_NOp,
             result.v_Np,
             result.v_HEp,
             result.v_Hp);
#endif
    }

  return s;
} /* cond_calc() */

/*
cond_get_result()
*/

cond_result *
cond_get_result(size_t idx, cond_workspace *w)
{
  return (&(w->cond_results[idx]));
} /* cond_get_result() */

/****************************************
 * INTERNAL ROUTINES                    *
 ****************************************/

static double
cond_direct(cond_result *result, cond_workspace *w)
{
  double sigma =
    1.0e6 * w->e_sq *
    (result->n_e / GSL_CONST_MKSA_MASS_ELECTRON / result->v_e +
     result->n_Op / COND_MASS_O / result->v_Op +
     result->n_O2p / COND_MASS_O2 / result->v_O2p +
     result->n_NOp / COND_MASS_NO / result->v_NOp +
     result->n_Hp / COND_MASS_H / result->v_Hp +
     result->n_HEp / COND_MASS_HE / result->v_HEp +
     result->n_Np / COND_MASS_N / result->v_Np);

  return sigma;
} /* cond_direct() */

static double
cond_pederson(cond_result *result, cond_workspace *w)
{
  double sigma =
    1.0e6 * w->e_sq *
    ((result->n_e * result->v_e / GSL_CONST_MKSA_MASS_ELECTRON) /
     (result->v_e_sq + result->w_e_sq) +
     (result->n_Op * result->v_Op / COND_MASS_O) /
     (result->v_Op_sq + result->w_Op_sq) +
     (result->n_O2p * result->v_O2p / COND_MASS_O2) /
     (result->v_O2p_sq + result->w_O2p_sq) +
     (result->n_NOp * result->v_NOp / COND_MASS_NO) /
     (result->v_NOp_sq + result->w_NOp_sq) +
     (result->n_Hp * result->v_Hp / COND_MASS_H) /
     (result->v_Hp_sq + result->w_Hp_sq) +
     (result->n_HEp * result->v_HEp / COND_MASS_HE) /
     (result->v_HEp_sq + result->w_HEp_sq) +
     (result->n_Np * result->v_Np / COND_MASS_N) /
     (result->v_Np_sq + result->w_Np_sq));

  return sigma;
} /* cond_pederson() */

static double
cond_hall(cond_result *result, cond_workspace *w)
{
  double sigma =
    1.0e6 * w->e_sq *
    ((result->n_e * result->w_e / GSL_CONST_MKSA_MASS_ELECTRON) /
     (result->v_e_sq + result->w_e_sq) -
     (result->n_Op * result->w_Op / COND_MASS_O) /
     (result->v_Op_sq + result->w_Op_sq) -
     (result->n_O2p * result->w_O2p / COND_MASS_O2) /
     (result->v_O2p_sq + result->w_O2p_sq) -
     (result->n_Hp * result->w_Hp / COND_MASS_H) /
     (result->v_Hp_sq + result->w_Hp_sq) -
     (result->n_HEp * result->w_HEp / COND_MASS_HE) /
     (result->v_HEp_sq + result->w_HEp_sq) -
     (result->n_Np * result->w_Np / COND_MASS_N) /
     (result->v_Np_sq + result->w_Np_sq) -
     (result->n_NOp * result->w_NOp / COND_MASS_NO) /
     (result->v_NOp_sq + result->w_NOp_sq));

  return sigma;
} /* cond_hall() */
