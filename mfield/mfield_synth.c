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

/* fill in internal coefficient vector of synthetic gauss coefficients */
int
mfield_synth_g(gsl_vector * g, mfield_workspace * w)
{
  mfield_synth_coeff *gptr;

  gsl_vector_set_zero(g);

  /* initialize coefficients g from test_gnm[] */
  for (gptr = &test_gnm[0]; gptr->n != 0; ++gptr)
    {
      size_t n = gptr->n;
      int m = gptr->m;
      size_t cidx = mfield_coeff_nmidx(n, m);

      mfield_set_mf(g, cidx, gptr->gnm, w);
      mfield_set_sv(g, cidx, gptr->dgnm, w);
      mfield_set_sa(g, cidx, gptr->ddgnm, w);
    }

  return 0;
}

/* replace with synthetic data for testing */
int
mfield_synth_replace(mfield_workspace *w)
{
  const size_t nmax = w->nmax_mf;
  const double t0 = w->epoch;
  size_t i, j;
  size_t plm_size = gsl_sf_legendre_array_n(nmax);
  double *Plm = malloc(plm_size * sizeof(double));
  double *dPlm = malloc(plm_size * sizeof(double));
  size_t c_size = w->p;
  gsl_vector *g = gsl_vector_alloc(c_size);

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
          double t1 = t - t0;
          double t2 = 0.5 * t1 * t1;
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          double sint = sin(theta);
          double cost = cos(theta);
          double X = 0.0, Y = 0.0, Z = 0.0;
          double ratio = w->R / r;
          size_t n;
          int m;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                          cost, Plm, dPlm);

          for (n = 1; n <= nmax; ++n)
            {
              int ni = (int) n;
              double rterm = pow(ratio, n + 2.0);

              for (m = 0; m <= ni; ++m)
                {
                  double c = cos(m * phi);
                  double s = sin(m * phi);
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

                  X += rterm * (gnm * c + hnm * s) * dPlm[pidx];
                  Y += rterm / sint * m * (gnm * s - hnm * c) * Plm[pidx];
                  Z -= (n + 1.0) * rterm *
                       (gnm * c + hnm * s) * Plm[pidx];
                }
            }

          mptr->Bx_nec[j] = X;
          mptr->By_nec[j] = Y;
          mptr->Bz_nec[j] = Z;
          mptr->F[j] = gsl_hypot3(X, Y, Z);

          /* no crustal/external field for synthetic data */
          mptr->Bx_model[j] = 0.0;
          mptr->By_model[j] = 0.0;
          mptr->Bz_model[j] = 0.0;

#if MFIELD_FIT_EULER
          /* rotate NEC vector to VFM frame */
          {
            double *q = &(mptr->q[4*j]);
            double B_nec[3], B_vfm[3];

            B_nec[0] = X;
            B_nec[1] = Y;
            B_nec[2] = Z;

            euler_nec2vfm(EULER_FLG_ZYX, alpha, beta, gamma, q, B_nec, B_vfm);

            mptr->Bx_vfm[j] = B_vfm[0];
            mptr->By_vfm[j] = B_vfm[1];
            mptr->Bz_vfm[j] = B_vfm[2];
          }
#endif
        }
    }

  free(Plm);
  free(dPlm);
  gsl_vector_free(g);

  return 0;
}
