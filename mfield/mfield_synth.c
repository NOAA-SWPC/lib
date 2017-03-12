/*
 * mfield_synth.c
 *
 * Routines for handling synthetic test case
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>

#include "euler.h"

#include "mfield.h"
#include "mfield_synth.h"

/* n m gnm dgnm ddgnm */
static mfield_synth_coeff test_gnm[] = {
  /*
   * need to start somewhat close to IGRF for main field, otherwise
   * Euler angles won't converge
   */
  { 1, 0, -30000.0, -25.0, 1.2 },
  { 1, -1, 5000.0, 30.0, -3.2 },
  { 2, 1, 3000.0, -4.0, 4.2 },
  { 3, -2, 250.0, -3.0, -6.2 },

  { 0, 0, 0.0, 0.0, 0.0 }
};

static int mfield_synth_calc(const double t, const double r, const double theta, const double phi,
                             const gsl_vector * g, double B[3], mfield_workspace *w);

/* fill in internal coefficient vector of synthetic gauss coefficients */
int
mfield_synth_g(gsl_vector * g, mfield_workspace * w)
{
  mfield_synth_coeff *gptr;

  gsl_vector_set_zero(g);

  /* initialize MF/SV/SA coefficients g from test_gnm[] */
  for (gptr = &test_gnm[0]; gptr->n != 0; ++gptr)
    {
      size_t n = gptr->n;
      int m = gptr->m;
      size_t cidx = mfield_coeff_nmidx(n, m);

      mfield_set_mf(g, cidx, gptr->gnm, w);
      mfield_set_sv(g, cidx, gptr->dgnm, w);
      mfield_set_sa(g, cidx, gptr->ddgnm, w);
    }

  /* fill in crustal field part of g with MF7 values */
  if (w->nmax_mf >= 16)
    {
      msynth_workspace *crust_p = msynth_mf7_read(MSYNTH_MF7_FILE);
      size_t n;

      for (n = 16; n <= w->nmax_mf; ++n)
        {
          int M = (int) n;
          int m;

          for (m = -M; m <= M; ++m)
            {
              size_t cidx = mfield_coeff_nmidx(n, m);
              size_t didx = msynth_nmidx(n, m, crust_p);
              mfield_set_mf(g, cidx, crust_p->c[didx], w);
            }
        }

      msynth_free(crust_p);
    }

  return 0;
}

/* replace with synthetic data for testing */
int
mfield_synth_replace(mfield_workspace *w)
{
  size_t i, j;
  gsl_vector *g = gsl_vector_alloc(w->p_int);

#if MFIELD_FIT_EULER
  /* Euler angles */
  const double alpha = -13.1 * M_PI / 180.0 * 1.0;
  const double beta = -5.2 * M_PI / 180.0 * 1.0;
  const double gamma = 3.4 * M_PI / 180.0 * 1.0;
#endif

  /* initialize synthetic gauss coefficients */
  mfield_synth_g(g, w);

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

      for (j = 0; j < mptr->n; ++j)
        {
          double t = satdata_epoch2year(mptr->t[j]);
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          double B[3];

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* synthesize magnetic field vector */
          mfield_synth_calc(t, r, theta, phi, g, B, w);

          mptr->Bx_nec[j] = B[0];
          mptr->By_nec[j] = B[1];
          mptr->Bz_nec[j] = B[2];
          mptr->F[j] = gsl_hypot3(B[0], B[1], B[2]);

          /* no crustal/external field for synthetic data */
          mptr->Bx_model[j] = 0.0;
          mptr->By_model[j] = 0.0;
          mptr->Bz_model[j] = 0.0;

#if MFIELD_FIT_EULER
          /* rotate NEC vector to VFM frame */
          {
            double *q = &(mptr->q[4*j]);
            double B_vfm[3];

            euler_nec2vfm(EULER_FLG_ZYX, alpha, beta, gamma, q, B, B_vfm);

            mptr->Bx_vfm[j] = B_vfm[0];
            mptr->By_vfm[j] = B_vfm[1];
            mptr->Bz_vfm[j] = B_vfm[2];
          }
#endif

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
            {
              t = satdata_epoch2year(mptr->t_ns[j]);
              mfield_synth_calc(t, mptr->r_ns[j], mptr->theta_ns[j], mptr->phi_ns[j], g, B, w);

              mptr->Bx_nec_ns[j] = B[0];
              mptr->By_nec_ns[j] = B[1];
              mptr->Bz_nec_ns[j] = B[2];
              mptr->F_ns[j] = gsl_hypot3(B[0], B[1], B[2]);

              mptr->Bx_model_ns[j] = 0.0;
              mptr->By_model_ns[j] = 0.0;
              mptr->Bz_model_ns[j] = 0.0;

#if MFIELD_FIT_EULER
              /* rotate NEC vector to VFM frame */
              {
                double *q = &(mptr->q_ns[4*j]);
                double B_vfm[3];

                euler_nec2vfm(EULER_FLG_ZYX, alpha, beta, gamma, q, B, B_vfm);

                mptr->Bx_vfm_ns[j] = B_vfm[0];
                mptr->By_vfm_ns[j] = B_vfm[1];
                mptr->Bz_vfm_ns[j] = B_vfm[2];
              }
#endif
            }

        }
    }

  gsl_vector_free(g);

  return 0;
}

/* compute B(r,theta,phi) using Gauss coefficients 'g' */
static int
mfield_synth_calc(const double t, const double r, const double theta, const double phi, const gsl_vector * g,
                  double B[3], mfield_workspace *w)
{
  int s = 0;
  const double t0 = w->epoch;
  const double t1 = t - t0;
  const double t2 = 0.5 * t1 * t1;
  const size_t nmax = w->nmax_mf;
  const double ratio = w->R / r;
  const double sint = sin(theta);
  const double cost = cos(theta);
  double rterm = ratio * ratio;
  double *Plm = w->Plm;
  double *dPlm = w->dPlm;
  size_t n;

  B[0] = 0.0;
  B[1] = 0.0;
  B[2] = 0.0;

  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax, cost, Plm, dPlm);

  for (n = 0; n <= nmax; ++n)
    {
      w->cosmphi[n] = cos(n * phi);
      w->sinmphi[n] = sin(n * phi);
    }

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;
      int m;

      /* (R/r)^{n+2} */
      rterm *= ratio;

      for (m = 0; m <= ni; ++m)
        {
          double c = w->cosmphi[m];
          double s = w->sinmphi[m];
          size_t pidx = gsl_sf_legendre_array_index(n, m);
          size_t cidx = mfield_coeff_nmidx(n, m);
          double gnm, hnm = 0.0;

          gnm = mfield_get_mf(g, cidx, w) +
                mfield_get_sv(g, cidx, w) * t1 +
                mfield_get_sa(g, cidx, w) * t2;

          if (m > 0)
            {
              cidx = mfield_coeff_nmidx(n, -m);
              hnm = mfield_get_mf(g, cidx, w) +
                    mfield_get_sv(g, cidx, w) * t1 +
                    mfield_get_sa(g, cidx, w) * t2;
            }

          B[0] += rterm * (gnm * c + hnm * s) * dPlm[pidx];
          B[1] += rterm / sint * m * (gnm * s - hnm * c) * Plm[pidx];
          B[2] -= (n + 1.0) * rterm * (gnm * c + hnm * s) * Plm[pidx];
        }
    }

  return s;
}
