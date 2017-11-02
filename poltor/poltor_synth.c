/*
 * poltor_synth.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_test.h>

#include "magdata.h"
#include "poltor.h"
#include "poltor_test.h"

/* internal potential coefficients: n m gnm */
static poltor_test_coeff test_gnm[] = {
  { 1, 0, -3.0 },
  { 1, -1, 5.8 },
  { 2, 1, 3.2 },
  { 3, -2, 20.1 },
  { 14, -5, 0.25 },

  { 0, 0, 0.0 }
};

/* external potential coefficients: n m knm */
static poltor_test_coeff test_knm[] = {
  { 1, 0, 15.6 },
  { 1, 1, 8.1 },
  { 2, 0, -2.2 },
  { 2, 1, 12.5 },
  { 2, 2, -2.7 },

  { 0, 0, 0.0 }
};

/* degree 0 shell poloidal coefficients: n m qnm */
static poltor_test_coeff test_qnm[] = {
  { 1, 0, -152.3 },
  { 1, -1, -18.1 },
  { 2, 0, -24.2 },
  { 2, 1, 125.5 },
  { 2, 2, -32.7 },

  { 0, 0, 0.0 }
};

/* toroidal coefficients: n m phinm */
static poltor_test_coeff test_phinm[] = {
  { 1, 0, 0.1 },
  { 1, 1, -0.1 },
  { 12, 0, -4.0 },
  { 13, -1, 12.0 },
  { 13, 2, -80.7 },
  { 25, -6, -1.1 },

  { 0, 0, 0.0 }
};

/* initialize coefficient vector with synthetic coefficients */
static int
poltor_synth_init(gsl_vector_complex *g, poltor_workspace *w)
{
  int s = 0;
  poltor_test_coeff *cptr;

  /* initialize internal coefficients from test_gnm[] */
  for (cptr = &test_gnm[0]; cptr->n != 0; ++cptr)
    {
      size_t n = cptr->n;
      int m = cptr->m;
      size_t cidx;
      gsl_complex val;

      if (n > w->nmax_int || abs(m) > (int) w->mmax_int)
        continue;

      cidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
      GSL_SET_COMPLEX(&val, cptr->gnm, 0.0);
      gsl_vector_complex_set(g, cidx, val);

      /*
       * to get a real-valued magnetic field, we require
       * g_{n,-m} = conj(g_{nm})
       */
      cidx = poltor_nmidx(POLTOR_IDX_PINT, n, -m, w);
      GSL_SET_COMPLEX(&val, cptr->gnm, 0.0);
      gsl_vector_complex_set(g, cidx, val);
    }

  /* initialize external coefficients from test_knm[] */
  for (cptr = &test_knm[0]; cptr->n != 0; ++cptr)
    {
      size_t n = cptr->n;
      int m = cptr->m;
      size_t cidx;
      gsl_complex val;

      if (n > w->nmax_ext || abs(m) > (int) w->mmax_ext)
        continue;

      cidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
      GSL_SET_COMPLEX(&val, cptr->gnm, 0.0);
      gsl_vector_complex_set(g, cidx, val);

      /*
       * to get a real-valued magnetic field, we require
       * k_{n,-m} = conj(k_{nm})
       */
      cidx = poltor_nmidx(POLTOR_IDX_PEXT, n, -m, w);
      GSL_SET_COMPLEX(&val, cptr->gnm, 0.0);
      gsl_vector_complex_set(g, cidx, val);
    }

  /* initialize shell poloidal coefficients from test_qnm[] */
  {
    size_t j;

    for (j = 0; j <= w->shell_J; ++j)
      {
        double fac = j + 1.0;

        for (cptr = &test_qnm[0]; cptr->n != 0; ++cptr)
          {
            size_t n = cptr->n;
            int m = cptr->m;
            size_t cidx;
            gsl_complex val;

            if (n > w->nmax_sh || abs(m) > (int) w->mmax_sh)
              continue;

            cidx = poltor_jnmidx(j, n, m, w);
            GSL_SET_COMPLEX(&val, fac * cptr->gnm, 0.0);
            gsl_vector_complex_set(g, cidx, val);

            /*
             * to get a real-valued magnetic field, we require
             * q_{n,-m} = conj(q_{nm})
             */
            cidx = poltor_jnmidx(j, n, -m, w);
            GSL_SET_COMPLEX(&val, fac * cptr->gnm, 0.0);
            gsl_vector_complex_set(g, cidx, val);
          }
      }
  }

  /* initialize toroidal coefficients from test_phinm[] */
  for (cptr = &test_phinm[0]; cptr->n != 0; ++cptr)
    {
      size_t n = cptr->n;
      int m = cptr->m;
      size_t cidx;
      gsl_complex val;

      if (n > w->nmax_tor || abs(m) > (int) w->mmax_tor)
        continue;

      cidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
      GSL_SET_COMPLEX(&val, cptr->gnm, 0.0);
      gsl_vector_complex_set(g, cidx, val);

      /*
       * to get a real-valued magnetic field, we require
       * phi_{n,-m} = conj(phi_{nm})
       */
      cidx = poltor_nmidx(POLTOR_IDX_TOR, n, -m, w);
      GSL_SET_COMPLEX(&val, cptr->gnm, 0.0);
      gsl_vector_complex_set(g, cidx, val);
    }

  return s;
} /* poltor_synth_init() */

/* calculate synthetic field value for single point */
static int
poltor_synth_calc(const double r, const double theta, const double phi,
                  gsl_vector_complex *g, double B[3], poltor_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax_max;
  const size_t mmax = w->mmax_max;
  const double invsint = 1.0 / sin(theta);
  const double cost = cos(theta);
  complex double X = 0.0, Y = 0.0, Z = 0.0;
  double ratio;
  size_t n, j;
  int m;

  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                  cost, w->Pnm, w->dPnm);

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= (int) mmax; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);

      for (n = GSL_MAX(m, 1); n <= nmax; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;
        }
    }

  /* internal potential field */
  ratio = w->R / r;
  for (n = 1; n <= w->nmax_int; ++n)
    {
      int ni = (int) GSL_MIN(n, w->mmax_int);
      double rterm = pow(ratio, n + 2.0);

      for (m = -ni; m <= ni; ++m)
        {
          size_t mabs = (size_t) abs(m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          size_t cidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          gsl_complex val = gsl_vector_complex_get(g, cidx);
          complex double gnm = GSL_REAL(val) + I * GSL_IMAG(val);
          complex double Ynm, dYnm;

          if (m >= 0)
            {
              Ynm = w->Ynm[pidx];
              dYnm = w->dYnm[pidx];
            }
          else
            {
              Ynm = conj(w->Ynm[pidx]);
              dYnm = conj(w->dYnm[pidx]);
            }

          X += rterm * gnm * dYnm;
          Y -= rterm * I * m * invsint * gnm * Ynm;
          Z -= (n + 1.0) * rterm * gnm * Ynm;
        }
    }

  /* external potential field */
  ratio = r / w->R;
  for (n = 1; n <= w->nmax_ext; ++n)
    {
      int ni = (int) GSL_MIN(n, w->mmax_ext);
      double rterm = pow(ratio, n - 1.0);

      for (m = -ni; m <= ni; ++m)
        {
          size_t mabs = (size_t) abs(m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          size_t cidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
          gsl_complex val = gsl_vector_complex_get(g, cidx);
          complex double knm = GSL_REAL(val) + I * GSL_IMAG(val);
          complex double Ynm, dYnm;

          if (m >= 0)
            {
              Ynm = w->Ynm[pidx];
              dYnm = w->dYnm[pidx];
            }
          else
            {
              Ynm = conj(w->Ynm[pidx]);
              dYnm = conj(w->dYnm[pidx]);
            }

          X += rterm * knm * dYnm;
          Y -= rterm * I * m * invsint * knm * Ynm;
          Z += n * rterm * knm * Ynm;
        }
    }

  /* shell poloidal field */
  ratio = w->R / r;
  for (j = 0; j <= w->shell_J; ++j)
    {
      for (n = 1; n <= w->nmax_sh; ++n)
        {
          int ni = (int) GSL_MIN(n, w->mmax_sh);
          double An, Bn, ABsum, dABsum;

          poltor_shell_An(n, j, r, &An, w);
          poltor_shell_Bn(n, j, r, &Bn, w);
          ABsum = An + Bn;
          dABsum = -(double)n * An + (n + 1.0) * Bn;

          for (m = -ni; m <= ni; ++m)
            {
              size_t mabs = (size_t) abs(m);
              size_t pidx = gsl_sf_legendre_array_index(n, mabs);
              size_t cidx = poltor_jnmidx(j, n, m, w);
              gsl_complex val = gsl_vector_complex_get(g, cidx);
              complex double gnm = GSL_REAL(val) + I * GSL_IMAG(val);
              complex double Ynm, dYnm;

              if (m >= 0)
                {
                  Ynm = w->Ynm[pidx];
                  dYnm = w->dYnm[pidx];
                }
              else
                {
                  Ynm = conj(w->Ynm[pidx]);
                  dYnm = conj(w->dYnm[pidx]);
                }

              X += ratio * gnm * dABsum * dYnm;
              Y -= ratio * I * m * invsint * gnm * dABsum * Ynm;
              Z += ratio * n * (n + 1.0) * gnm * ABsum * Ynm;
            }
        }
    }

  /* shell toroidal field */
  for (n = 1; n <= w->nmax_tor; ++n)
    {
      int ni = (int) GSL_MIN(n, w->mmax_tor);

      for (m = -ni; m <= ni; ++m)
        {
          size_t mabs = (size_t) abs(m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          size_t cidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
          gsl_complex val = gsl_vector_complex_get(g, cidx);
          complex double phinm = GSL_REAL(val) + I * GSL_IMAG(val);
          complex double Ynm, dYnm;

          if (m >= 0)
            {
              Ynm = w->Ynm[pidx];
              dYnm = w->dYnm[pidx];
            }
          else
            {
              Ynm = conj(w->Ynm[pidx]);
              dYnm = conj(w->dYnm[pidx]);
            }

          X += I * m * invsint * phinm * Ynm;
          Y += phinm * dYnm;
        }
    }

  B[0] = creal(X);
  B[1] = creal(Y);
  B[2] = creal(Z);

  return s;
} /* poltor_synth_calc() */

/* replace with synthetic data for testing */
int
poltor_synth(poltor_workspace *w)
{
  int s = 0;
  size_t i, j;
  gsl_vector_complex *g = gsl_vector_complex_calloc(w->p);
  magdata_list *list = w->data;

  s = poltor_synth_init(g, w);

  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);

      for (j = 0; j < mptr->n; ++j)
        {
          double B[3];

          s += poltor_synth_calc(mptr->r[j], mptr->theta[j], mptr->phi[j], g, B, w);

          mptr->Bx_nec[j] = B[0];
          mptr->By_nec[j] = B[1];
          mptr->Bz_nec[j] = B[2];
          mptr->Bx_model[j] = 0.0;
          mptr->By_model[j] = 0.0;
          mptr->Bz_model[j] = 0.0;

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS))
            {
              double B_ns[3];

              s += poltor_synth_calc(mptr->r_ns[j], mptr->theta_ns[j], mptr->phi_ns[j], g, B_ns, w);

              mptr->Bx_nec_ns[j] = B_ns[0];
              mptr->By_nec_ns[j] = B_ns[1];
              mptr->Bz_nec_ns[j] = B_ns[2];
              mptr->Bx_model_ns[j] = 0.0;
              mptr->By_model_ns[j] = 0.0;
              mptr->Bz_model_ns[j] = 0.0;
            }
        }
    }

  gsl_vector_complex_free(g);

  return 0;
} /* poltor_synth() */

int
poltor_synth_print(poltor_workspace *w)
{
  int s = 0;
  const double relerr = 1.0e-6;
  poltor_test_coeff *cptr;

  /* print internal poloidal coefficients */
  fprintf(stderr, "Internal poloidal coefficients:\n");
  for (cptr = &test_gnm[0]; cptr->n != 0 && cptr->n <= w->nmax_int; ++cptr)
    {
      size_t n = cptr->n;
      int m = cptr->m;
      size_t cidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
      gsl_complex coef = poltor_get(cidx, w);
      double gnm = GSL_REAL(coef);

      fprintf(stderr, "g(%2zu,%2d) = %12g (%5g) relerr = %15.6e\n",
              n,
              m,
              gnm,
              cptr->gnm,
              fabs(gnm - cptr->gnm) / cptr->gnm);

      gsl_test_rel(gnm, cptr->gnm, relerr, "g(%zu,%d)", n, m);
    }

  /* print external poloidal coefficients */
  fprintf(stderr, "External poloidal coefficients:\n");
  for (cptr = &test_knm[0]; cptr->n != 0 && cptr->n <= w->nmax_ext; ++cptr)
    {
      size_t n = cptr->n;
      int m = cptr->m;
      size_t cidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
      gsl_complex coef = poltor_get(cidx, w);
      double knm = GSL_REAL(coef);

      fprintf(stderr, "k(%2zu,%2d) = %12g (%5g) relerr = %15.6e\n",
              n,
              m,
              knm,
              cptr->gnm,
              fabs(knm - cptr->gnm) / cptr->gnm);

      gsl_test_rel(knm, cptr->gnm, relerr, "k(%zu,%d)", n, m);
    }

  /* print shell poloidal coefficients */
  fprintf(stderr, "Shell poloidal coefficients:\n");
  {
    size_t j;

    for (j = 0; j <= w->shell_J; ++j)
      {
        double fac = j + 1.0;

        for (cptr = &test_qnm[0]; cptr->n != 0 && cptr->n <= w->nmax_sh; ++cptr)
          {
            size_t n = cptr->n;
            int m = cptr->m;
            size_t cidx = poltor_jnmidx(j, n, m, w);
            gsl_complex coef = poltor_get(cidx, w);
            double qnm = GSL_REAL(coef);

            /* only output j=0 terms */
            if (j == 0)
              {
                fprintf(stderr, "q(%2zu,%2d) = %12g (%5g) relerr = %15.6e\n",
                        n,
                        m,
                        qnm,
                        cptr->gnm,
                        fabs(qnm - cptr->gnm) / cptr->gnm);
              }

            gsl_test_rel(qnm, fac * cptr->gnm, relerr, "j=%zu q(%zu,%d)", j, n, m);
          }
      }
  }

  /* print shell toroidal coefficients */
  fprintf(stderr, "Shell toroidal coefficients:\n");
  for (cptr = &test_phinm[0]; cptr->n != 0 && cptr->n <= w->nmax_tor; ++cptr)
    {
      size_t n = cptr->n;
      int m = cptr->m;
      size_t cidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
      gsl_complex coef = poltor_get(cidx, w);
      double phinm = GSL_REAL(coef);

      fprintf(stderr, "phi(%2zu,%2d) = %10g (%5g) relerr = %15.6e\n",
              n,
              m,
              phinm,
              cptr->gnm,
              fabs(phinm - cptr->gnm) / cptr->gnm);

      gsl_test_rel(phinm, cptr->gnm, relerr, "phi(%zu,%d)", n, m);
    }

  return s;
}
