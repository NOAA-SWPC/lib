/*
 * poltor.c
 *
 * Decompose satellite magnetic field measurements to toroidal
 * currents mapped onto a thin shell in the ionosphere
 *
 * Calling sequence:
 *
 * 1. magdata_alloc     - allocate data structure for internal storage
 *                        of satellite magnetic measurements
 * 2. magdata_add       - add single datum to internal data structure
 * 3. magdata_calc      - calculate spatial weights for all data
 * 4. poltor_alloc      - allocate a workspace, specify thin current shell
 *                        radii, maximum spherical harmonic degree,
 *                        and previously constructed data structure
 * 5. poltor_calc       - calculate LS matrix, rhs vector, and store both
 *                        to a file for future quick manipulation (LLS.dat)
 * 6. poltor_solve      - regularize and solve LS system for coefficients
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <complex.h>
#include <sys/time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <common/common.h>
#include <common/oct.h>

#include "lls.h"
#include "magdata.h"
#include "poltor.h"

static int poltor_row(const double r, const double theta, const double phi,
                      const double r_ns, const double theta_ns, const double phi_ns,
                      gsl_vector_complex *x, gsl_vector_complex *y, gsl_vector_complex *z,
                      gsl_vector_complex *dx, gsl_vector_complex *dy, gsl_vector_complex *dz,
                      poltor_workspace *w);
static int poltor_row_pint(const double r, const double theta,
                           const double r_ns, const double theta_ns,
                           gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
                           gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
                           poltor_workspace *w);
static int poltor_row_pext(const double r, const double theta,
                           const double r_ns, const double theta_ns,
                           gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
                           gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
                           poltor_workspace *w);
static int poltor_row_psh(const double r, const double theta,
                          const double r_ns, const double theta_ns,
                          gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
                          gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
                          poltor_workspace *w);
static int poltor_row_tor(const double r, const double theta,
                          const double r_ns, const double theta_ns,
                          gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
                          gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
                          poltor_workspace *w);
static int poltor_eval_chi(const double r, const double theta, const double phi,
                           const size_t idx, double *chi, poltor_workspace *w);

/*
poltor_alloc()
  Allocate a poltor workspace

Inputs: params - parameters
          R        - reference radius (km)
          nmax_int - maximum internal spherical harmonic degree
          mmax_int - maximum internal spherical harmonic order
          nmax_ext - maximum external spherical harmonic degree
          mmax_ext - maximum external spherical harmonic order
          nmax_tor - maximum spherical harmonic degree for B_tor
          mmax_tor - maximum spherical harmonic order for B_tor
          data     - satellite data to invert
*/

poltor_workspace *
poltor_alloc(const poltor_parameters *params)
{
  poltor_workspace *w;
  size_t psize; /* size of Legendre arrays */
  size_t n;

  w = calloc(1, sizeof(poltor_workspace));
  if (!w)
    return 0;

  w->params = *params;

  w->R = params->R;
  w->b = params->b;
  w->d = params->d;
  w->rmin = params->rmin;
  w->rmax = params->rmax;
  w->nmax_int = params->nmax_int;
  w->mmax_int = params->mmax_int;
  w->nmax_ext = params->nmax_ext;
  w->mmax_ext = params->mmax_ext;
  w->nmax_sh = params->nmax_sh;
  w->mmax_sh = params->mmax_sh;
  w->nmax_tor = params->nmax_tor;
  w->mmax_tor = params->mmax_tor;
  w->shell_J = params->shell_J;
  w->data = params->data;
  w->alpha_int = params->alpha_int;
  w->alpha_sh = params->alpha_sh;
  w->alpha_tor = params->alpha_tor;
  w->flags = params->flags;

  w->rmid = 0.5 * (w->rmin + w->rmax);

  w->nmax_max = GSL_MAX(w->nmax_int,
                  GSL_MAX(w->nmax_ext,
                    GSL_MAX(w->nmax_sh, w->nmax_tor)));
  w->mmax_max = GSL_MAX(w->mmax_int,
                  GSL_MAX(w->mmax_ext,
                    GSL_MAX(w->mmax_sh, w->mmax_tor)));

  if (w->data)
    {
      w->n = w->data->nres;
    }
  else
    {
      /* minimal workspace for poltor_eval_B() */
      w->n = 6;
    }

  w->base_int = calloc(1, (w->nmax_int + 1) * sizeof(size_t));
  w->base_ext = calloc(1, (w->nmax_ext + 1) * sizeof(size_t));
  w->base_sh = calloc(1, (w->nmax_sh + 1) * sizeof(size_t));
  w->base_tor = calloc(1, (w->nmax_tor + 1) * sizeof(size_t));

  /*
   * precompute base offsets for internal (n,m) indexing, and count total
   * internal coefficients
   */
  w->p_pint = 0;
  for (n = 1; n <= w->nmax_int; ++n)
    {
      int ns = (int) GSL_MIN(n, w->mmax_int);

      w->base_int[n] = w->p_pint;
      w->p_pint += 2*ns + 1;
    }

  /*
   * precompute base offsets for external (n,m) indexing, and count total
   * external coefficients
   */
  w->p_pext = 0;
  for (n = 1; n <= w->nmax_ext; ++n)
    {
      int ns = (int) GSL_MIN(n, w->mmax_ext);

      w->base_ext[n] = w->p_pext;
      w->p_pext += 2*ns + 1;
    }

  /*
   * precompute base offsets for poloidal shell (n,m) indexing, and count total
   * internal coefficients
   */
  w->nnm_sh = 0;
  for (n = 1; n <= w->nmax_sh; ++n)
    {
      int ns = (int) GSL_MIN(n, w->mmax_sh);

      w->base_sh[n] = w->nnm_sh;
      w->nnm_sh += 2*ns + 1;
    }

  /* compute total poloidal shell coefficients accouting for Taylor expansion */
  w->p_psh = w->nnm_sh * (w->shell_J + 1);

  /*
   * precompute base offsets for B_tor (n,m) indexing, and count total
   * coefficients
   */
  w->p_tor = 0;
  for (n = 1; n <= w->nmax_tor; ++n)
    {
      int ns = (int) GSL_MIN(n, w->mmax_tor);

      w->base_tor[n] = w->p_tor;
      w->p_tor += 2*ns + 1;
    }

  w->p = w->p_pint + w->p_pext + w->p_psh + w->p_tor;

  w->pint_offset = 0;
  w->pext_offset = w->pint_offset + w->p_pint;
  w->psh_offset = w->pext_offset + w->p_pext;
  w->tor_offset = w->psh_offset + w->p_psh;

  if (w->n > 0)
    {
      /*
       * Need to fold in measurements in blocks, otherwise 'A'
       * matrix will be too large to fit in memory for some cases;
       * multiply by 6 for X,Y,Z,DX,DY,DZ
       */
      w->nblock = GSL_MIN(w->n, 6 * POLTOR_BLOCK_SIZE);

      w->A = gsl_matrix_complex_alloc(w->nblock, w->p);
      if (!w->A)
        {
          fprintf(stderr, "poltor_alloc: cannot allocate design matrix: %s\n",
                  strerror(errno));
          poltor_free(w);
          return 0;
        }

      w->rhs = gsl_vector_complex_alloc(w->nblock);
      if (!w->rhs)
        {
          fprintf(stderr, "poltor_alloc: cannot allocate rhs vector: %s\n",
                  strerror(errno));
          poltor_free(w);
          return 0;
        }

      w->weights = gsl_vector_alloc(w->nblock);
      if (!w->weights)
        {
          fprintf(stderr, "poltor_alloc: cannot allocate weights vector: %s\n",
                  strerror(errno));
          poltor_free(w);
          return 0;
        }

      w->lls_workspace_p = lls_complex_alloc(w->nblock, w->p);
      if (!w->lls_workspace_p)
        {
          fprintf(stderr, "poltor_alloc: cannot allocate lls workspace: %s\n",
                  strerror(errno));
          poltor_free(w);
          return 0;
        }

      w->dof = (int) w->n - (int) w->p;
    }

  w->c = gsl_vector_complex_alloc(w->p);
  if (!w->c)
    {
      fprintf(stderr, "poltor_alloc: cannot allocate coefficient vector: %s\n",
              strerror(errno));
      poltor_free(w);
      return 0;
    }

  w->L = gsl_vector_alloc(w->p);
  if (!w->L)
    {
      fprintf(stderr, "poltor_alloc: cannot allocate Tikhonov vector: %s\n",
              strerror(errno));
      poltor_free(w);
      return 0;
    }

  /* initialize to identity matrix */
  gsl_vector_set_all(w->L, 1.0);

  w->residuals = gsl_vector_complex_alloc(w->n);
  if (!w->residuals)
    {
      fprintf(stderr, "poltor_alloc: cannot allocate residual vector: %s\n",
              strerror(errno));
      poltor_free(w);
      return 0;
    }

  psize = gsl_sf_legendre_array_n(w->nmax_max);

  w->Pnm = malloc(psize * sizeof(double));
  w->dPnm = malloc(psize * sizeof(double));
  w->Pnm2 = malloc(psize * sizeof(double));
  w->dPnm2 = malloc(psize * sizeof(double));
  w->Ynm = malloc(psize * sizeof(complex double));
  w->dYnm = malloc(psize * sizeof(complex double));
  w->Ynm2 = malloc(psize * sizeof(complex double));
  w->dYnm2 = malloc(psize * sizeof(complex double));
  if (!w->Pnm || !w->dPnm || !w->Pnm2 || !w->dPnm2 ||
      !w->Ynm || !w->dYnm || !w->Ynm2 || !w->dYnm2)
    {
      fprintf(stderr, "poltor_alloc: cannot allocate legendre arrays: %s\n",
              strerror(errno));
      poltor_free(w);
      return 0;
    }

  w->cquad_workspace_p = gsl_integration_cquad_workspace_alloc(1000);
  if (!w->cquad_workspace_p)
    {
      fprintf(stderr, "poltor_alloc: cannot allocate integration workspace: %s\n",
              strerror(errno));
      poltor_free(w);
      return 0;
    }

  w->nreg = 200;
  w->reg_param = gsl_vector_alloc(w->nreg);
  w->rho = gsl_vector_alloc(w->nreg);
  w->eta = gsl_vector_alloc(w->nreg);

  return w;
} /* poltor_alloc() */

void
poltor_free(poltor_workspace *w)
{
  if (w->base_int)
    free(w->base_int);

  if (w->base_ext)
    free(w->base_ext);

  if (w->base_sh)
    free(w->base_sh);

  if (w->base_tor)
    free(w->base_tor);

  if (w->A)
    gsl_matrix_complex_free(w->A);

  if (w->rhs)
    gsl_vector_complex_free(w->rhs);

  if (w->weights)
    gsl_vector_free(w->weights);

  if (w->c)
    gsl_vector_complex_free(w->c);

  if (w->L)
    gsl_vector_free(w->L);

  if (w->residuals)
    gsl_vector_complex_free(w->residuals);

  if (w->Pnm)
    free(w->Pnm);

  if (w->dPnm)
    free(w->dPnm);

  if (w->Pnm2)
    free(w->Pnm2);

  if (w->dPnm2)
    free(w->dPnm2);

  if (w->Ynm)
    free(w->Ynm);

  if (w->dYnm)
    free(w->dYnm);

  if (w->Ynm2)
    free(w->Ynm2);

  if (w->dYnm2)
    free(w->dYnm2);

  if (w->lls_workspace_p)
    lls_complex_free(w->lls_workspace_p);

  if (w->cquad_workspace_p)
    gsl_integration_cquad_workspace_free(w->cquad_workspace_p);

  if (w->reg_param)
    gsl_vector_free(w->reg_param);

  if (w->rho)
    gsl_vector_free(w->rho);

  if (w->eta)
    gsl_vector_free(w->eta);

  free(w);
}

/*
poltor_calc()
  Build LS system: A^H W A and A^H W b

Notes:
1) On output, w->L contains diag(L) the regularization matrix
*/

int
poltor_calc(poltor_workspace *w)
{
  int s = 0;
  const char *lls_file = "LLS.dat";
  magdata *data = w->data;

  fprintf(stderr, "poltor_calc: nvec   = %zu\n", data->nvec);
  fprintf(stderr, "poltor_calc: ngrad  = %zu\n", data->ngrad);
  fprintf(stderr, "poltor_calc: nres   = %zu\n", data->nres);
  fprintf(stderr, "poltor_calc: ncoeff = %zu\n", w->p);
  fprintf(stderr, "poltor_calc: nres / ncoeff = %.1f\n",
          (double) data->nres / (double) w->p);

  /* initialize lls module */
  lls_complex_reset(w->lls_workspace_p);

  fprintf(stderr, "poltor_calc: computing regularization matrix...");

#if !POLTOR_SYNTH_DATA

  if (w->params.flags & POLTOR_FLG_REGULARIZE)
    {
      /* construct L = diag(L) */
      s = poltor_regularize(w->L, w);
      if (s)
        return s;
    }

#endif

  fprintf(stderr, "done\n");

  fprintf(stderr, "poltor_calc: building LS system...");

  s += poltor_build_ls(1, w);
  if (s)
    return s;

  fprintf(stderr, "done\n");

  fprintf(stderr, "poltor_calc: writing LS system to %s...",
          lls_file);
  s = lls_complex_save(lls_file, w->lls_workspace_p);
  fprintf(stderr, "done\n");

  return s;
} /* poltor_calc() */

/*
poltor_solve()
  Solve LS system previously constructed by poltor_calc()

Inputs: w - workspace

Notes:
1) On output, w->c contains coefficient vector

2) On output, w->chisq contains the residual sum of squares r^T r
*/

int
poltor_solve(poltor_workspace *w)
{
  int s = 0;
  struct timeval tv0, tv1;
  const double dof = (double) w->data->nres - (double) w->p;
  double lambda = 0.0; /* regularization parameter */

#if !POLTOR_SYNTH_DATA

  if (w->params.flags & POLTOR_FLG_REGULARIZE)
    {
      const char *lcurve_file = "lcurve.dat";

      /* construct L = diag(L); this must be done again here
       * in case user loads LLS system from a file
       */
      s = poltor_regularize(w->L, w);
      if (s)
        return s;

      /* calculate L-curve */
      fprintf(stderr, "poltor_solve: calculating L-curve...");
      lls_complex_lcurve(w->reg_param, w->rho, w->eta, w->lls_workspace_p);
      fprintf(stderr, "done\n");

      /* find L-corve corner */
      fprintf(stderr, "poltor_solve: calculating L-curve corner...");
      gsl_multifit_linear_lcorner(w->rho, w->eta, &(w->reg_idx));
      lambda = gsl_vector_get(w->reg_param, w->reg_idx);
      fprintf(stderr, "done (lambda = %f)\n", lambda);

      fprintf(stderr, "poltor_solve: writing L-curve to %s...", lcurve_file);
      poltor_print_lcurve(lcurve_file, w);
      fprintf(stderr, "done\n");
    }

#endif

  fprintf(stderr, "poltor_solve: solving LS system...");
  gettimeofday(&tv0, NULL);

  /* solve A^H W A c = A^H W b */
  s = lls_complex_solve(lambda, w->c, w->lls_workspace_p);
  if (s)
    {
      fprintf(stderr, "poltor_solve: error solving system: %d\n", s);
      return s;
    }

  /* backtransform standard form system to get original solution vector */
  lls_complex_btransform(w->L, w->c, w->lls_workspace_p);

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, residual = %.12e, chisq/dof = %.12e)\n",
          time_diff(tv0, tv1),
          w->lls_workspace_p->residual,
          w->lls_workspace_p->chisq / dof);

  fprintf(stderr, "poltor_solve: matrix condition number = %f\n",
          w->lls_workspace_p->cond);

  w->chisq = w->lls_workspace_p->chisq;

  return s;
} /* poltor_solve() */

/*
poltor_eval_J_shell()
  Calculate shell toroidal current vector at a given point

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        J     - (output) shell current density vector in units of A/m^2
                (J_r,J_theta,J_phi) - J_r = 0 for toroidal currents
        w     - workspace

Notes:
1) J should be initialized prior to calling this function, since
on output, J <- J + J_tor
*/

int
poltor_eval_J_shell(const double r, const double theta, const double phi,
                    double J[3], poltor_workspace *w)
{
  int s = 0;
  size_t n, j;
  int m;
  const double invsint = 1.0 / sin(theta);
  const double ratio = (r - w->rmid) / w->R;
  double rterm = 1.0; /* [(r - r_0) / R]^j */
  complex double Jt = 0.0, Jp = 0.0;

  /* compute legendre functions */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, w->nmax_sh,
                                  cos(theta), w->Pnm, w->dPnm);

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= w->mmax_sh; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);

      for (n = GSL_MAX(m, 1); n <= w->nmax_sh; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;
        }
    }

  for (j = 0; j <= w->shell_J; ++j)
    {
      for (n = 1; n <= w->nmax_sh; ++n)
        {
          int M = (int) GSL_MIN(n, w->mmax_sh);

          for (m = -M; m <= M; ++m)
            {
              int mabs = abs(m);
              size_t nmidx = poltor_jnmidx(j, n, m, w);
              size_t pidx = gsl_sf_legendre_array_index(n, mabs);
              gsl_complex cnm = gsl_vector_complex_get(w->c, nmidx);
              complex double qnm = GSL_REAL(cnm) + I * GSL_IMAG(cnm);
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

              /* compute (theta,phi) components of J */
              Jt -= rterm * I * m * invsint * qnm * Ynm;
              Jp += rterm * qnm * dYnm;
            }
        }

      rterm *= ratio;
    }

  /* convert to A/m^2 */
  J[1] += creal(Jt) / w->R / POLTOR_MU_0 * 1.0e-3;
  J[2] += creal(Jp) / w->R / POLTOR_MU_0 * 1.0e-3;

  return s;
} /* poltor_eval_J_shell */

/*
poltor_eval_K_tor()
  Calculate toroidal sheet current vector at a given point
on the spherical shell r = b

Inputs: b     - radius of spherical shell (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        K     - (output) sheet current vector in units of kA/km
                (K_r,K_theta,K_phi) - K_r = 0 for toroidal currents
        w     - workspace

Notes:
1) K should be initialized prior to calling this function, since
on output, K <- K + K_tor
*/

int
poltor_eval_K_tor(const double b, const double theta, const double phi,
                  double J[3], poltor_workspace *w)
{
  int s = 0;
  size_t n, m;
  const double invsint = 1.0 / sin(theta);
  complex double Jt = 0.0, Jp = 0.0;

  /* compute legendre functions */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, w->nmax_int,
                                  cos(theta), w->Pnm, w->dPnm);

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= w->mmax_int; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);

      for (n = GSL_MAX(m, 1); n <= w->nmax_int; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;
        }
    }

  for (n = 1; n <= w->nmax_int; ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax_int);
      int m;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          gsl_complex cnm = poltor_get(nmidx, w);
          double gnm = GSL_REAL(cnm);
          double qnm = -(2.0*n + 1.0) / (double)n * gnm *
                       pow(w->R / b, n + 3.0);
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

          /* compute (theta,phi) components of J */
          Jt -= I * m * invsint * qnm * Ynm;
          Jp += qnm * dYnm;
        }
    }

  /* convert to kA/km */
  J[1] += (b / w->R) * creal(Jt) / POLTOR_MU_0;
  J[2] += (b / w->R) * creal(Jp) / POLTOR_MU_0;

  return s;
} /* poltor_eval_K_tor */

/*
poltor_eval_K_ext()
  Calculate toroidal sheet current vector at a given point
on the spherical shell r = b due to external sources (ie:
use external coefficients in calculation)

Inputs: r     - radius of shell b (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        K     - (output) sheet current vector in units of kA/km
                (K_r,K_theta,K_phi) - K_r = 0 for toroidal currents
        w     - workspace

Notes:
1) K should be initialized prior to calling this function, since
on output, K <- K + K_ext
*/

int
poltor_eval_K_ext(const double r, const double theta, const double phi,
                  double J[3], poltor_workspace *w)
{
  int s = 0;
  size_t n, m;
  const size_t nmax = w->nmax_ext;
  const size_t mmax = w->mmax_ext;
  const double invsint = 1.0 / sin(theta);
  complex double Jt = 0.0, Jp = 0.0;
  const double ratio = r / w->R;
  double rterm = pow(ratio, -2.0);

  /* compute legendre functions */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                  cos(theta), w->Pnm, w->dPnm);

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= mmax; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);

      for (n = GSL_MAX(m, 1); n <= nmax; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;
        }
    }

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;
      double nfac = (2.0 * n + 1.0) / (n + 1.0);

      rterm *= ratio;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          gsl_complex cnm = gsl_vector_complex_get(w->c, nmidx);
          complex double qnm;
          complex double Ynm, dYnm;

          qnm = nfac * rterm * (GSL_REAL(cnm) + I * GSL_IMAG(cnm));

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

          /* compute (theta,phi) components of J */
          Jt -= I * m * invsint * qnm * Ynm;
          Jp += qnm * dYnm;
        }
    }

  /* convert to kA/km */
  J[1] += ratio * creal(Jt) / POLTOR_MU_0;
  J[2] += ratio * creal(Jp) / POLTOR_MU_0;

  return s;
} /* poltor_eval_K_ext */

/*
poltor_eval_chi_int()
  Calculate toroidal current scalar function chi at given point
on the spherical shell r = b

chi = -Q b / mu_0 

Inputs: b     - radius of spherical shell (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        chi   - (output) current scalar in units of kA
        w     - workspace

Return: success/error
*/

int
poltor_eval_chi_int(const double b, const double theta, const double phi,
                    double *chi, poltor_workspace *w)
{
  int s = 0;
  size_t n, m;
  complex double sum = 0.0;
  const size_t nmax = w->nmax_int;
  const size_t mmax = w->mmax_int;
  const double ratio = w->R / b;
  double rterm = pow(ratio, 2.0);

  /* compute legendre functions */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                  cos(theta), w->Pnm, w->dPnm);

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= mmax; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);

      for (n = GSL_MAX(m, 1); n <= nmax; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;
        }
    }

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;
      double nfac = -(2.0 * n + 1.0) / (double) n;

      /* (R/b)^(n+2) */
      rterm *= ratio;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          gsl_complex cnm = gsl_vector_complex_get(w->c, nmidx);
          complex double qnm;
          complex double Ynm;

          qnm = nfac * rterm * (GSL_REAL(cnm) + I * GSL_IMAG(cnm));

          if (m >= 0)
            Ynm = w->Ynm[pidx];
          else
            Ynm = conj(w->Ynm[pidx]);

          sum += qnm * Ynm;
        }
    }

  /* current function = -Qr/mu0 (units of kA) */
  sum = -sum / POLTOR_MU_0 * b;
  *chi = creal(sum);

  return s;
} /* poltor_eval_chi_int() */

/*
poltor_eval_chi_sh()
  Calculate toroidal current scalar function chi at given point
on the spherical shell r = d

chi = -Q d / mu_0 

Inputs: r     - radius (km) (not currently used)
        theta - colatitude (radians)
        phi   - longitude (radians)
        chi   - (output) current scalar in units of kA
        w     - workspace

Return: success/error
*/

int
poltor_eval_chi_sh(const double r, const double theta, const double phi,
                   double *chi, poltor_workspace *w)
{
  int s = 0;

  s = poltor_eval_chi(w->d, theta, phi, POLTOR_IDX_PSH, chi, w);

  return s;
} /* poltor_eval_chi_sh() */

/*
poltor_eval_chi_ext()
  Calculate toroidal current scalar function chi at given point
on the spherical shell r = b

chi = -Q r / mu_0 

Inputs: b          - radius of current shell (km)
        theta      - colatitude (radians)
        phi        - longitude (radians)
        chi        - (output) current scalar in units of kA
        w          - workspace

Return: success/error
*/

int
poltor_eval_chi_ext(const double b, const double theta, const double phi,
                    double *chi, poltor_workspace *w)
{
  int s = 0;
  size_t n, m;
  complex double sum = 0.0;
  const size_t nmax = w->nmax_ext;
  const size_t mmax = w->mmax_ext;
  const double ratio = b / w->R;
  double rterm = pow(ratio, -2.0);

  /* compute legendre functions */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                  cos(theta), w->Pnm, w->dPnm);

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= mmax; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);

      for (n = GSL_MAX(m, 1); n <= nmax; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;
        }
    }

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;
      double nfac = (2.0 * n + 1.0) / (n + 1.0);

      /* (b/R)^(n-1) */
      rterm *= ratio;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          gsl_complex cnm = gsl_vector_complex_get(w->c, nmidx);
          complex double qnm;
          complex double Ynm;

          qnm = nfac * rterm * (GSL_REAL(cnm) + I * GSL_IMAG(cnm));

          if (m >= 0)
            Ynm = w->Ynm[pidx];
          else
            Ynm = conj(w->Ynm[pidx]);

          sum += qnm * Ynm;
        }
    }

  /* current function = -Qr/mu0 (units of kA) */
  sum = -sum / POLTOR_MU_0 * b;
  *chi = creal(sum);

  return s;
} /* poltor_eval_chi_ext() */

/*
poltor_eval_B()
  Calculate total magnetic field model at a given point

Inputs: r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - longitude (radians)
        B     - (output) NEC magnetic field model in units of nT
                (B_x,B_y,B_z,|B|)
        w     - workspace

Return: success/error
*/

int
poltor_eval_B(const double r, const double theta, const double phi,
              double B[4], poltor_workspace *w)
{
  int s = 0;
  gsl_vector_complex_view vx, vy, vz;
  gsl_complex valx, valy, valz;

  /* views of 3 rows of the design matrix for X,Y,Z */
  vx = gsl_matrix_complex_row(w->A, 0);
  vy = gsl_matrix_complex_row(w->A, 1);
  vz = gsl_matrix_complex_row(w->A, 2);

  /* build these 3 rows of design matrix at this point */
  poltor_row(r, theta, phi, 0.0, 0.0, 0.0,
             &vx.vector, &vy.vector, &vz.vector,
             NULL, NULL, NULL, w);

  gsl_blas_zdotu(&vx.vector, w->c, &valx);
  gsl_blas_zdotu(&vy.vector, w->c, &valy);
  gsl_blas_zdotu(&vz.vector, w->c, &valz);

  B[0] = GSL_REAL(valx);
  B[1] = GSL_REAL(valy);
  B[2] = GSL_REAL(valz);
  B[3] = gsl_hypot3(B[0], B[1], B[2]);

  return s;
} /* poltor_eval_B() */

/*
poltor_eval_B_all()
  Calculate magnetic field model at a given point

Inputs: r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - longitude (radians)
        B_int - (output) internal poloidal NEC magnetic field model in units of nT
                (B_x,B_y,B_z,|B|)
        B_ext - (output) external poloidal NEC magnetic field model in units of nT
                (B_x,B_y,B_z,|B|)
        B_sh  - (output) shell poloidal NEC magnetic field model in units of nT
                (B_x,B_y,B_z,|B|)
        B_tor - (output) shell toroidal NEC magnetic field model in units of nT
                (B_x,B_y,B_z,|B|)
        w     - workspace

Return: success/error
*/

int
poltor_eval_B_all(const double r, const double theta, const double phi,
                  double B_int[4], double B_ext[4], double B_sh[4], double B_tor[4],
                  poltor_workspace *w)
{
  int s = 0;
  gsl_vector_complex_view vx, vy, vz;
  gsl_vector_complex_view c, sx, sy, sz;
  gsl_complex valx, valy, valz;
  size_t i;

  /* initialize to 0 */
  for (i = 0; i < 3; ++i)
    {
      B_int[i] = 0.0;
      B_ext[i] = 0.0;
      B_sh[i] = 0.0;
      B_tor[i] = 0.0;
    }

  /* views of 3 rows of the design matrix for X,Y,Z */
  vx = gsl_matrix_complex_row(w->A, 0);
  vy = gsl_matrix_complex_row(w->A, 1);
  vz = gsl_matrix_complex_row(w->A, 2);

  /* build these 3 rows of design matrix at this point */
  poltor_row(r, theta, phi, 0.0, 0.0, 0.0,
             &vx.vector, &vy.vector, &vz.vector,
             NULL, NULL, NULL, w);

  if (w->p_pint > 0)
    {
      c = gsl_vector_complex_subvector(w->c, w->pint_offset, w->p_pint);
      sx = gsl_vector_complex_subvector(&vx.vector, w->pint_offset, w->p_pint);
      sy = gsl_vector_complex_subvector(&vy.vector, w->pint_offset, w->p_pint);
      sz = gsl_vector_complex_subvector(&vz.vector, w->pint_offset, w->p_pint);

      gsl_blas_zdotu(&sx.vector, &c.vector, &valx);
      gsl_blas_zdotu(&sy.vector, &c.vector, &valy);
      gsl_blas_zdotu(&sz.vector, &c.vector, &valz);

      B_int[0] = GSL_REAL(valx);
      B_int[1] = GSL_REAL(valy);
      B_int[2] = GSL_REAL(valz);
    }

  if (w->p_pext > 0)
    {
      c = gsl_vector_complex_subvector(w->c, w->pext_offset, w->p_pext);
      sx = gsl_vector_complex_subvector(&vx.vector, w->pext_offset, w->p_pext);
      sy = gsl_vector_complex_subvector(&vy.vector, w->pext_offset, w->p_pext);
      sz = gsl_vector_complex_subvector(&vz.vector, w->pext_offset, w->p_pext);

      gsl_blas_zdotu(&sx.vector, &c.vector, &valx);
      gsl_blas_zdotu(&sy.vector, &c.vector, &valy);
      gsl_blas_zdotu(&sz.vector, &c.vector, &valz);

      B_ext[0] = GSL_REAL(valx);
      B_ext[1] = GSL_REAL(valy);
      B_ext[2] = GSL_REAL(valz);
    }

  if (w->p_psh > 0)
    {
      c = gsl_vector_complex_subvector(w->c, w->psh_offset, w->p_psh);
      sx = gsl_vector_complex_subvector(&vx.vector, w->psh_offset, w->p_psh);
      sy = gsl_vector_complex_subvector(&vy.vector, w->psh_offset, w->p_psh);
      sz = gsl_vector_complex_subvector(&vz.vector, w->psh_offset, w->p_psh);

#if 0
      /* XXX */
      {
        /* set degrees [nmin,nmax] to zero */
        const size_t nmin = 1;
        const size_t nmax = 6;
        size_t j, n;
        for (j = 0; j <= w->shell_J; ++j)
          {
            for (n = nmin; n <= nmax; ++n)
              {
                int M = (int) GSL_MIN(n, w->mmax_sh);
                int m;

                for (m = -M; m <= M; ++m)
                  {
                    size_t nmidx = poltor_jnmidx(j, n, m, w);
                    gsl_vector_complex_set(w->c, nmidx, GSL_COMPLEX_ZERO);
                  }
              }
          }
      }
#endif

      gsl_blas_zdotu(&sx.vector, &c.vector, &valx);
      gsl_blas_zdotu(&sy.vector, &c.vector, &valy);
      gsl_blas_zdotu(&sz.vector, &c.vector, &valz);

      B_sh[0] = GSL_REAL(valx);
      B_sh[1] = GSL_REAL(valy);
      B_sh[2] = GSL_REAL(valz);
    }

  if (w->p_tor > 0)
    {
      c = gsl_vector_complex_subvector(w->c, w->tor_offset, w->p_tor);
      sx = gsl_vector_complex_subvector(&vx.vector, w->tor_offset, w->p_tor);
      sy = gsl_vector_complex_subvector(&vy.vector, w->tor_offset, w->p_tor);
      sz = gsl_vector_complex_subvector(&vz.vector, w->tor_offset, w->p_tor);

      gsl_blas_zdotu(&sx.vector, &c.vector, &valx);
      gsl_blas_zdotu(&sy.vector, &c.vector, &valy);
      gsl_blas_zdotu(&sz.vector, &c.vector, &valz);

      B_tor[0] = GSL_REAL(valx);
      B_tor[1] = GSL_REAL(valy);
      B_tor[2] = GSL_REAL(valz);
    }

  /* compute |B| */
  B_int[3] = gsl_hypot3(B_int[0], B_int[1], B_int[2]);
  B_ext[3] = gsl_hypot3(B_ext[0], B_ext[1], B_ext[2]);
  B_sh[3] = gsl_hypot3(B_sh[0], B_sh[1], B_sh[2]);
  B_tor[3] = gsl_hypot3(B_tor[0], B_tor[1], B_tor[2]);

  return s;
} /* poltor_eval_B_all() */

int
poltor_write(const char *filename, poltor_workspace *w)
{
  int s = 0;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "poltor_write: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  fwrite(&(w->params), sizeof(poltor_parameters), 1, fp);
  gsl_vector_complex_fwrite(fp, w->c);

  fclose(fp);

  return s;
} /* poltor_write() */

poltor_workspace *
poltor_read(const char *filename)
{
  poltor_workspace *w;
  poltor_parameters params;
  FILE *fp;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "poltor_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  fread(&params, sizeof(poltor_parameters), 1, fp);

  params.data = NULL;

  w = poltor_alloc(&params);

  gsl_vector_complex_fread(fp, w->c);

  fclose(fp);

  return w;
} /* poltor_read() */

/*
poltor_get()
  Return coefficient corresponding to a given index.
The index should be obtained from poltor_nmidx() and so
it will already account for internal/external etc

Inputs: cidx - coefficient index from poltor_nmidx()
        w    - workspace
*/

gsl_complex
poltor_get(const size_t cidx, const poltor_workspace *w)
{
  gsl_complex coeff;
  
  if (cidx >= w->c->size)
    GSL_SET_COMPLEX(&coeff, 0.0, 0.0);
  else
    coeff = gsl_vector_complex_get(w->c, cidx);

  return coeff;
}

double
poltor_spectrum_int(const size_t n, poltor_workspace *w)
{
  int ni = (int) GSL_MIN(n, w->mmax_int);
  int m;
  double sum = 0.0;

  if (n > w->nmax_int)
    return 0.0;

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
      gsl_complex coef = poltor_get(cidx, w);
      double gnm = GSL_REAL(coef);

      sum += gnm * gnm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
}

double
poltor_spectrum_ext(const size_t n, poltor_workspace *w)
{
  int ni = (int) GSL_MIN(n, w->mmax_ext);
  int m;
  double sum = 0.0;

  if (n > w->nmax_ext)
    return 0.0;

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
      gsl_complex coef = poltor_get(cidx, w);
      double knm = GSL_REAL(coef);

      sum += knm * knm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
} /* poltor_spectrum_ext() */

double
poltor_spectrum_sh(const size_t n, poltor_workspace *w)
{
  int ni = (int) GSL_MIN(n, w->mmax_sh);
  int m;
  double sum = 0.0;

  if (n > w->nmax_sh)
    return 0.0;

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = poltor_nmidx(POLTOR_IDX_PSH, n, m, w);
      gsl_complex coef = poltor_get(cidx, w);
      double qnm = GSL_REAL(coef);

      sum += qnm * qnm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
} /* poltor_spectrum_sh() */

double
poltor_spectrum_tor(const size_t n, poltor_workspace *w)
{
  int ni = (int) GSL_MIN(n, w->mmax_tor);
  int m;
  double sum = 0.0;

  if (n > w->nmax_tor)
    return 0.0;

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
      gsl_complex coef = poltor_get(cidx, w);
      double phinm = GSL_REAL(coef);

      sum += phinm * phinm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
} /* poltor_spectrum_tor() */

int
poltor_print_spectrum(const char *filename, poltor_workspace *w)
{
  const size_t nmax = w->nmax_max;
  size_t n;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "poltor_print_spectrum: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  n = 1;
  fprintf(fp, "# Field %zu: spherical harmonic degree\n", n++);
  fprintf(fp, "# Field %zu: internal spectrum g_n (nT^2)\n", n++);
  fprintf(fp, "# Field %zu: external spectrum k_n (nT^2)\n", n++);
  fprintf(fp, "# Field %zu: shell poloidal spectrum q_n (nT^2)\n", n++);
  fprintf(fp, "# Field %zu: shell toroidal spectrum phi_n (nT^2)\n", n++);

  for (n = 1; n <= nmax; ++n)
    {
      double gn = poltor_spectrum_int(n, w);
      double kn = poltor_spectrum_ext(n, w);
      double qn = poltor_spectrum_sh(n, w);
      double phin = poltor_spectrum_tor(n, w);

      fprintf(fp, "%zu %.12e %.12e %.12e %.12e\n",
              n,
              gn,
              kn,
              qn,
              phin);
    }

  fclose(fp);

  return 0;
} /* poltor_print_spectrum() */

int
poltor_print_lcurve(const char *filename, const poltor_workspace *w)
{
  FILE *fp;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "poltor_print_lcurve: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: regularization parameter lambda\n", i++);
  fprintf(fp, "# Field %zu: residual norm ||b - A x||\n", i++);
  fprintf(fp, "# Field %zu: solution norm ||L x||\n", i++);

  for (i = 0; i < w->nreg; ++i)
    {
      fprintf(fp, "%.12e %.12e %.12e\n",
              gsl_vector_get(w->reg_param, i),
              gsl_vector_get(w->rho, i),
              gsl_vector_get(w->eta, i));
    }

  fprintf(fp, "\n\n");

  fprintf(fp, "%.12e %.12e %.12e\n",
          gsl_vector_get(w->reg_param, w->reg_idx),
          gsl_vector_get(w->rho, w->reg_idx),
          gsl_vector_get(w->eta, w->reg_idx));

  fclose(fp);

  return 0;
}

/*
poltor_theta()
  Determine colatitude (either geocentric or QD)

Inputs: idx - index of data point
        w   - workspace

Return: colatitude in radians
*/

double
poltor_theta(const size_t idx, const poltor_workspace *w)
{
  if (w->flags & POLTOR_FLG_QD_HARMONICS)
    return (M_PI / 2.0 - w->data->qdlat[idx] * M_PI / 180.0);

  return w->data->theta[idx];
} /* poltor_theta() */

/*
poltor_theta_ns()
  Determine colatitude (either geocentric or QD)

Inputs: idx - index of data point
        w   - workspace

Return: colatitude in radians
*/

double
poltor_theta_ns(const size_t idx, const poltor_workspace *w)
{
  if (w->flags & POLTOR_FLG_QD_HARMONICS)
    return (M_PI / 2.0 - w->data->qdlat_ns[idx] * M_PI / 180.0);

  return w->data->theta_ns[idx];
} /* poltor_theta_ns() */

/***********************************************
 * INTERNAL ROUTINES                           *
 ***********************************************/

/*
poltor_build_ls()
  Build LS design matrix

Inputs: fold - if set to 1, LS matrix and rhs vector are folded
               into A^H W A and A^H W b
               if set to 0, it is assumed we already solved the
               LS system, and now wish to compute residual norm
               ||b - A x|| for L-curve analysis
        w    - workspace

Notes:
1) If fold = 1, on output, w->lls_workspace_p is updated with A^H W A and A^H W b,
ready for solving

2) If fold = 0, on output, w->chisq contains ||sqrt(W) (b - A x)||^2

3) If fold = 0, on output, w->residuals contains sqrt(W) * (b - A x)

4) Flagged outliers are not added to matrix
*/

int
poltor_build_ls(const int fold, poltor_workspace *w)
{
  int s = 0;
  magdata *data = w->data;
  size_t i;
  struct timeval tv0, tv1;
  size_t data_idx = 0; /* index into 'data' */
  size_t res_idx = 0;  /* index into 'residuals' */
  int nleft;           /* number of vector data left to fold in */

  fprintf(stderr, "\n");

  w->chisq = 0.0;

  nleft = (int) data->n;
  while (nleft > 0)
    {
      gsl_matrix_complex_view A;
      gsl_vector_complex_view rhs;
      gsl_vector_view wts;
      size_t ridx = 0; /* index of current row in A */

      /* number of vector measurements to process */
      size_t ndata = (size_t) GSL_MIN(POLTOR_BLOCK_SIZE, nleft);

      fprintf(stderr, "\t adding measurements [%zu,%zu] of %zu into matrix...",
              data_idx, data_idx + ndata - 1, data->n);
      gettimeofday(&tv0, NULL);

      /* initialize matrix to 0 */
      gsl_matrix_complex_set_zero(w->A);

      for (i = 0; i < ndata; ++i, ++data_idx)
        {
          double theta = poltor_theta(data_idx, w);
          double theta_ns = poltor_theta_ns(data_idx, w);
          gsl_vector_complex *x = NULL,
                             *y = NULL,
                             *z = NULL,
                             *dx = NULL,
                             *dy = NULL,
                             *dz = NULL;
          gsl_vector_complex_view vx, vy, vz, vdx, vdy, vdz;
          gsl_complex val;
          double B[4], dB[4];

          /* compute residual: B_i - B_main - B_crust - B_ext */
          magdata_residual(data_idx, B, data);

          if (!gsl_finite(B[0]) || !gsl_finite(B[1]) || !gsl_finite(B[2]))
            {
              fprintf(stderr, "ACK1\n");
            }

          if (data->flags[data_idx] & (MAGDATA_FLG_DX_NS|MAGDATA_FLG_DY_NS|MAGDATA_FLG_DZ_NS))
            {
              /* compute dB residual: B_i_ns - B_i */
              magdata_residual_dB_ns(data_idx, dB, data);
              if (!gsl_finite(dB[0]) || !gsl_finite(dB[1]) || !gsl_finite(dB[2]))
                {
                  fprintf(stderr, "ACK2\n");
                }
            }

          if (data->flags[data_idx] & MAGDATA_FLG_X)
            {
              /* weights vector */
              gsl_vector_set(w->weights, ridx, POLTOR_WEIGHT_X * data->weights[data_idx]);

              /* rhs vector */
              GSL_SET_COMPLEX(&val, B[0], 0.0);
              gsl_vector_complex_set(w->rhs, ridx, val);

              vx = gsl_matrix_complex_row(w->A, ridx++);
              x = &vx.vector;
            }

          if (data->flags[data_idx] & MAGDATA_FLG_DX_NS)
            {
              gsl_vector_set(w->weights, ridx, POLTOR_WEIGHT_DX * data->weights[data_idx]);

              GSL_SET_COMPLEX(&val, dB[0], 0.0);
              gsl_vector_complex_set(w->rhs, ridx, val);

              vdx = gsl_matrix_complex_row(w->A, ridx++);
              dx = &vdx.vector;
            }

          if (data->flags[data_idx] & MAGDATA_FLG_Y)
            {
              /* weights vector */
              gsl_vector_set(w->weights, ridx, POLTOR_WEIGHT_Y * data->weights[data_idx]);

              /* rhs vector */
              GSL_SET_COMPLEX(&val, B[1], 0.0);
              gsl_vector_complex_set(w->rhs, ridx, val);

              vy = gsl_matrix_complex_row(w->A, ridx++);
              y = &vy.vector;
            }

          if (data->flags[data_idx] & MAGDATA_FLG_DY_NS)
            {
              gsl_vector_set(w->weights, ridx, POLTOR_WEIGHT_DY * data->weights[data_idx]);

              GSL_SET_COMPLEX(&val, dB[1], 0.0);
              gsl_vector_complex_set(w->rhs, ridx, val);

              vdy = gsl_matrix_complex_row(w->A, ridx++);
              dy = &vdy.vector;
            }

          if (data->flags[data_idx] & MAGDATA_FLG_Z)
            {
              /* weights vector */
              gsl_vector_set(w->weights, ridx, POLTOR_WEIGHT_Z * data->weights[data_idx]);

              /* rhs vector */
              GSL_SET_COMPLEX(&val, B[2], 0.0);
              gsl_vector_complex_set(w->rhs, ridx, val);

              vz = gsl_matrix_complex_row(w->A, ridx++);
              z = &vz.vector;
            }

          if (data->flags[data_idx] & MAGDATA_FLG_DZ_NS)
            {
              gsl_vector_set(w->weights, ridx, POLTOR_WEIGHT_DZ * data->weights[data_idx]);

              GSL_SET_COMPLEX(&val, dB[2], 0.0);
              gsl_vector_complex_set(w->rhs, ridx, val);

              vdz = gsl_matrix_complex_row(w->A, ridx++);
              dz = &vdz.vector;
            }

          /* build these 6 rows of design matrix */
          poltor_row(data->r[data_idx], theta, data->phi[data_idx],
                     data->r_ns[data_idx], theta_ns, data->phi_ns[data_idx],
                     x, y, z, dx, dy, dz, w);
        }

      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      /* add design matrix to least squares system */

      A = gsl_matrix_complex_submatrix(w->A, 0, 0, ridx, w->p);
      rhs = gsl_vector_complex_subvector(w->rhs, 0, ridx);
      wts = gsl_vector_subvector(w->weights, 0, ridx);

      if (fold)
        {
          /* convert LS system to standard form */
          fprintf(stderr, "\t converting (A,b) to standard form...");
          gettimeofday(&tv0, NULL);
          lls_complex_stdform(&A.matrix, &rhs.vector, &wts.vector, w->L, w->lls_workspace_p);
          gettimeofday(&tv1, NULL);
          fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

          /* compute A^H A and A^H b */
          fprintf(stderr, "\t folding into A^H A and A^H b...");
          gettimeofday(&tv0, NULL);
          lls_complex_fold(&A.matrix, &rhs.vector, w->lls_workspace_p);
          gettimeofday(&tv1, NULL);
          fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
        }
      else
        {
          gsl_vector_complex_view resv = gsl_vector_complex_subvector(w->residuals, res_idx, ridx);
          gsl_complex one, mone;
          double rnorm;

          GSL_SET_COMPLEX(&one, 1.0, 0.0);
          GSL_SET_COMPLEX(&mone, -1.0, 0.0);

          fprintf(stderr, "\t computing ||sqrt(W) (b - Ac)||^2...");

          /* calculate r = b - A c with previously computed coefficients */
          gsl_vector_complex_memcpy(&resv.vector, &rhs.vector);
          gsl_blas_zgemv(CblasNoTrans, mone, &A.matrix, w->c, one, &resv.vector);

          /* compute weighted residual term in chi^2: r_i' = sqrt(w_i) r_i */
          for (i = 0; i < ridx; ++i)
            {
              double wi = gsl_vector_get(&wts.vector, i);
              gsl_complex ri = gsl_vector_complex_get(&resv.vector, i);

              gsl_vector_complex_set(&resv.vector, i, gsl_complex_mul_real(ri, sqrt(wi)));
            }

          /* compute || sqrt(W) (b - A c) || */
          rnorm = gsl_blas_dznrm2(&resv.vector);

          w->chisq += rnorm * rnorm;

          fprintf(stderr, "done (chi^2 = %f)\n", w->chisq);

          res_idx += ridx;
        }

      nleft -= (int) ndata;
    }

  return s;
} /* poltor_build_ls() */

/*
poltor_row()
  Build 6 rows of LS design matrix for X, Y, Z, DX, DY, DZ components

Inputs: r      - geocentric radius (km)
        theta  - geocentric colatitude (radians)
        phi    - longitude (radians)
        r_ns   - geocentric radius for along-track point (km)
        theta_ns - geocentric colatitude for along-track point (radians)
        phi_ns - longitude for along-track point (radians)
        x      - (output) row corresponding to X (or NULL)
        y      - (output) row corresponding to Y (or NULL)
        z      - (output) row corresponding to Z (or NULL)
        dx     - (output) row corresponding to DX (or NULL)
        dy     - (output) row corresponding to DY (or NULL)
        dz     - (output) row corresponding to DZ (or NULL)
        w      - workspace

Return: success/error

Notes: row vectors should be initialized to 0 prior to calling function
*/

static int
poltor_row(const double r, const double theta, const double phi,
           const double r_ns, const double theta_ns, const double phi_ns,
           gsl_vector_complex *x, gsl_vector_complex *y, gsl_vector_complex *z,
           gsl_vector_complex *dx, gsl_vector_complex *dy, gsl_vector_complex *dz,
           poltor_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax_max;
  const size_t mmax = w->mmax_max;
  const int grad = dx || dy || dz;
  size_t n, m;

  /* compute legendre functions */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                  cos(theta), w->Pnm, w->dPnm);

  if (grad)
    {
      gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                      cos(theta_ns), w->Pnm2, w->dPnm2);
    }

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= mmax; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);
      complex double expimphi_ns = cos(m * phi_ns) + I * sin(m * phi_ns);

      for (n = GSL_MAX(m, 1); n <= nmax; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;

          if (grad)
            {
              w->Ynm2[pidx] = w->Pnm2[pidx] * expimphi_ns;
              w->dYnm2[pidx] = w->dPnm2[pidx] * expimphi_ns;
            }
        }
    }

  /* build row for B_pol^i coefficients */
  if (w->p_pint > 0)
    poltor_row_pint(r, theta, r_ns, theta_ns, x, y, z, dx, dy, dz, w);

  /* build row for B_pol^e coefficients */
  if (w->p_pext > 0)
    poltor_row_pext(r, theta, r_ns, theta_ns, x, y, z, dx, dy, dz, w);

  /* build row for B_pol^{sh} coefficients */
  if (w->p_psh > 0)
    {
#if 0
      double phi_deg = wrap180(phi * 180.0 / M_PI);
      double qdlat = 90.0 - theta * 180.0 / M_PI;

      /* exclude EEJ region from shell poloidal calculation */
      if (!(fabs(qdlat) <= 10.0 && (phi_deg >= -105.0 && phi_deg <= 105.0)))
#endif
        poltor_row_psh(r, theta, r_ns, theta_ns, x, y, z, dx, dy, dz, w);
    }

  /* build row for B_tor coefficients */
  if (w->p_tor > 0)
    poltor_row_tor(r, theta, r_ns, theta_ns, x, y, z, dx, dy, dz, w);

  return s;
} /* poltor_row() */

/*
poltor_row_pint()
  Build a row of the LS design matrix corresponding to poloidal
internal magnetic field B_pol^i

Inputs: r     - radius (km)
        theta - colatitude (radians)
        vx    - (output) where to store matrix row X component (or NULL)
        vy    - (output) where to store matrix row Y component (or NULL)
        vz    - (output) where to store matrix row Z component (or NULL)
        vdx   - (output) where to store matrix row DX component (or NULL)
        vdy   - (output) where to store matrix row DY component (or NULL)
        vdz   - (output) where to store matrix row DZ component (or NULL)
        w     - workspace

Return: success/error

Notes:
1) w->Ynm, w->dYnm must be initialized prior to calling this function
*/

static int
poltor_row_pint(const double r, const double theta,
                const double r_ns, const double theta_ns,
                gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
                gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
                poltor_workspace *w)
{
  int s = 0;
  size_t n;
  gsl_complex val;
  const complex double invisint = I / sin(theta);
  const complex double invisint2 = I / sin(theta_ns);
  const double ratio = w->R / r;
  const double ratio2 = w->R / r_ns;
  double rterm = ratio * ratio;
  double rterm2 = ratio2 * ratio2;

  /*
   * construct 3 rows of the design matrix corresponding to the
   * X,Y,Z data observation
   */
  for (n = 1; n <= w->nmax_int; ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax_int);
      int m;

      /* (R/r)^(n+2) */
      rterm *= ratio;

      /* (R/r_ns)^(n+2) */
      rterm2 *= ratio2;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          complex double X, Y, Z, Ynm, dYnm, Ynm2, dYnm2;

          if (m >= 0)
            {
              Ynm = w->Ynm[pidx];
              dYnm = w->dYnm[pidx];
              Ynm2 = w->Ynm2[pidx];
              dYnm2 = w->dYnm2[pidx];
            }
          else
            {
              Ynm = conj(w->Ynm[pidx]);
              dYnm = conj(w->dYnm[pidx]);
              Ynm2 = conj(w->Ynm2[pidx]);
              dYnm2 = conj(w->dYnm2[pidx]);
            }

          /* compute (X,Y,Z) components of B_{lm} */
          if (vx)
            {
              X = rterm * dYnm;
              GSL_SET_COMPLEX(&val, creal(X), cimag(X));
              gsl_vector_complex_set(vx, nmidx, val);
            }

          if (vdx)
            {
              X = rterm2 * dYnm2 - rterm * dYnm;
              GSL_SET_COMPLEX(&val, creal(X), cimag(X));
              gsl_vector_complex_set(vdx, nmidx, val);
            }

          if (vy)
            {
              Y = -rterm * m * invisint * Ynm;
              GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
              gsl_vector_complex_set(vy, nmidx, val);
            }

          if (vdy)
            {
              Y = -m * (rterm2 * invisint2 * Ynm2 -
                        rterm * invisint * Ynm);
              GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
              gsl_vector_complex_set(vdy, nmidx, val);
            }

          if (vz)
            {
              Z = -rterm * (n + 1.0) * Ynm;
              GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
              gsl_vector_complex_set(vz, nmidx, val);
            }

          if (vdz)
            {
              Z = -(n + 1.0) * (rterm2 * Ynm2 - rterm * Ynm);
              GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
              gsl_vector_complex_set(vdz, nmidx, val);
            }
        }
    }

  return s;
} /* poltor_row_pint() */

/*
poltor_row_pext()
  Build a row of the LS design matrix corresponding to poloidal
external magnetic field B_pol^e

Inputs: r     - radius (km)
        theta - colatitude (radians)
        vx    - (output) where to store matrix row X component (or NULL)
        vy    - (output) where to store matrix row Y component (or NULL)
        vz    - (output) where to store matrix row Z component (or NULL)
        w     - workspace

Return: success/error

Notes:
1) w->Ynm, and w->dYnm must be initialized prior to calling this function
*/

static int
poltor_row_pext(const double r, const double theta,
                const double r_ns, const double theta_ns,
                gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
                gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
                poltor_workspace *w)
{
  int s = 0;
  double ratio = r / w->R;
  double rterm = 1.0; /* (r / R)^(n-1) */
  double ratio2 = r_ns / w->R;
  double rterm2 = 1.0; /* (r_ns / R)^(n-1) */
  const complex double invisint = I / sin(theta);
  const complex double invisint2 = I / sin(theta_ns);
  size_t n;
  gsl_complex val;

  /*
   * construct 3 rows of the design matrix corresponding to the
   * X,Y,Z data observation
   */
  for (n = 1; n <= w->nmax_ext; ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax_ext);
      int m;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          complex double X, Y, Z, Ynm, dYnm, Ynm2, dYnm2;

          if (m >= 0)
            {
              Ynm = w->Ynm[pidx];
              dYnm = w->dYnm[pidx];
              Ynm2 = w->Ynm2[pidx];
              dYnm2 = w->dYnm2[pidx];
            }
          else
            {
              Ynm = conj(w->Ynm[pidx]);
              dYnm = conj(w->dYnm[pidx]);
              Ynm2 = conj(w->Ynm2[pidx]);
              dYnm2 = conj(w->dYnm2[pidx]);
            }

          if (vx)
            {
              X = rterm * dYnm;
              GSL_SET_COMPLEX(&val, creal(X), cimag(X));
              gsl_vector_complex_set(vx, nmidx, val);
            }

          if (vdx)
            {
              X = rterm2 * dYnm2 - rterm * dYnm;
              GSL_SET_COMPLEX(&val, creal(X), cimag(X));
              gsl_vector_complex_set(vdx, nmidx, val);
            }

          if (vy)
            {
              Y = -rterm * m * invisint * Ynm;
              GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
              gsl_vector_complex_set(vy, nmidx, val);
            }

          if (vdy)
            {
              Y = -m * (rterm2 * invisint2 * Ynm2 -
                        rterm * invisint * Ynm);
              GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
              gsl_vector_complex_set(vdy, nmidx, val);
            }

          if (vz)
            {
              Z = rterm * n * Ynm;
              GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
              gsl_vector_complex_set(vz, nmidx, val);
            }

          if (vdz)
            {
              Z = n * (rterm2 * Ynm2 - rterm * Ynm);
              GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
              gsl_vector_complex_set(vdz, nmidx, val);
            }
        }

      /* (r/R)^{n-1} */
      rterm *= ratio;

      /* (r_ns/R)^{n-1} */
      rterm2 *= ratio2;
    }

  return s;
} /* poltor_row_pext() */

/*
poltor_row_psh()
  Build a row of the LS design matrix corresponding to poloidal
magnetic field B_pol^{sh}

Inputs: r     - radius (km)
        theta - colatitude (radians)
        vx    - (output) where to store matrix row X component
        vy    - (output) where to store matrix row Y component
        vz    - (output) where to store matrix row Z component
        w     - workspace

Return: success/error

Notes:
1) w->Ynm, and w->dYnm must be initialized prior to calling this function
*/

static int
poltor_row_psh(const double r, const double theta,
               const double r_ns, const double theta_ns,
               gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
               gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
               poltor_workspace *w)
{
  int s = 0;
  size_t n, j;
  const double ratio = w->R / r;
  const double ratio2 = w->R / r_ns;
  const complex double invisint = I / sin(theta);
  const complex double invisint2 = I / sin(theta_ns);
  gsl_complex val;

  /* outer loop over Taylor series */
  for (j = 0; j <= w->shell_J; ++j)
    {
      for (n = 1; n <= w->nmax_sh; ++n)
        {
          int M = (int) GSL_MIN(n, w->mmax_sh);
          int m;
          double An, Bn, ABsum, dABsum, ABsum2, dABsum2;

          /* compute Taylor series integrals */
          s = poltor_shell_An(n, j, r, &An, w);
          s += poltor_shell_Bn(n, j, r, &Bn, w);
          if (s)
            return s;

          ABsum = An + Bn;
          dABsum = -(double)n * An + (n + 1.0) * Bn;

          s = poltor_shell_An(n, j, r_ns, &An, w);
          s += poltor_shell_Bn(n, j, r_ns, &Bn, w);
          if (s)
            return s;

          ABsum2 = An + Bn;
          dABsum2 = -(double)n * An + (n + 1.0) * Bn;

          for (m = -M; m <= M; ++m)
            {
              int mabs = abs(m);
              size_t nmidx = poltor_jnmidx(j, n, m, w);
              size_t pidx = gsl_sf_legendre_array_index(n, mabs);
              complex double X, Y, Z, Ynm, dYnm, Ynm2, dYnm2;

              if (m >= 0)
                {
                  Ynm = w->Ynm[pidx];
                  dYnm = w->dYnm[pidx];
                  Ynm2 = w->Ynm2[pidx];
                  dYnm2 = w->dYnm2[pidx];
                }
              else
                {
                  Ynm = conj(w->Ynm[pidx]);
                  dYnm = conj(w->dYnm[pidx]);
                  Ynm2 = conj(w->Ynm2[pidx]);
                  dYnm2 = conj(w->dYnm2[pidx]);
                }

              /* compute (X,Y,Z) components of B_{lm} */

              if (vx)
                {
                  X = ratio * dABsum * dYnm;
                  GSL_SET_COMPLEX(&val, creal(X), cimag(X));
                  gsl_vector_complex_set(vx, nmidx, val);
                }

              if (vdx)
                {
                  X = ratio2 * dABsum2 * dYnm2 - ratio * dABsum * dYnm;
                  GSL_SET_COMPLEX(&val, creal(X), cimag(X));
                  gsl_vector_complex_set(vdx, nmidx, val);
                }

              if (vy)
                {
                  Y = -ratio * m * invisint * dABsum * Ynm;
                  GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
                  gsl_vector_complex_set(vy, nmidx, val);
                }

              if (vdy)
                {
                  Y = -m * (ratio2 * invisint2 * dABsum2 * Ynm2 -
                            ratio * invisint * dABsum * Ynm);
                  GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
                  gsl_vector_complex_set(vdy, nmidx, val);
                }

              if (vz)
                {
                  Z = ratio * n * (n + 1.0) * ABsum * Ynm;
                  GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
                  gsl_vector_complex_set(vz, nmidx, val);
                }

              if (vdz)
                {
                  Z = n * (n + 1.0) * (ratio2 * ABsum2 * Ynm2 - ratio * ABsum * Ynm);
                  GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
                  gsl_vector_complex_set(vdz, nmidx, val);
                }
            }
        }
    }

  return s;
} /* poltor_row_psh() */

/*
poltor_row_tor()
  Build a row of the LS design matrix corresponding to toroidal
magnetic field B_tor

Inputs: r     - radius (km)
        theta - colatitude (radians)
        vx    - (output) where to store matrix row X component (or NULL)
        vy    - (output) where to store matrix row Y component (or NULL)
        vz    - (output) where to store matrix row Z component (or NULL)
        w     - workspace

Return: success/error

Notes:
1) w->Ynm, and w->dYnm must be initialized prior to calling this function
*/

static int
poltor_row_tor(const double r, const double theta,
               const double r_ns, const double theta_ns,
               gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
               gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
               poltor_workspace *w)
{
  int s = 0;
  const complex double invisint = I / sin(theta);
  const complex double invisint2 = I / sin(theta_ns);
  size_t n;
  gsl_complex val;

  /*
   * construct 3 rows of the design matrix corresponding to the
   * X,Y,Z data observation
   */
  for (n = 1; n <= w->nmax_tor; ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax_tor);
      int m;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          complex double X, Y, Z, Ynm, dYnm, Ynm2, dYnm2;

          if (m >= 0)
            {
              Ynm = w->Ynm[pidx];
              dYnm = w->dYnm[pidx];
              Ynm2 = w->Ynm2[pidx];
              dYnm2 = w->dYnm2[pidx];
            }
          else
            {
              Ynm = conj(w->Ynm[pidx]);
              dYnm = conj(w->dYnm[pidx]);
              Ynm2 = conj(w->Ynm2[pidx]);
              dYnm2 = conj(w->dYnm2[pidx]);
            }

          /* compute (r,theta,phi) components of B_{lm} */

          if (vx)
            {
              X = m * invisint * Ynm;
              GSL_SET_COMPLEX(&val, creal(X), cimag(X));
              gsl_vector_complex_set(vx, nmidx, val);
            }

          if (vdx)
            {
              X = m * (invisint2 * Ynm2 - invisint * Ynm);
              GSL_SET_COMPLEX(&val, creal(X), cimag(X));
              gsl_vector_complex_set(vdx, nmidx, val);
            }

          if (vy)
            {
              Y = dYnm;
              GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
              gsl_vector_complex_set(vy, nmidx, val);
            }

          if (vdy)
            {
              Y = dYnm2 - dYnm;
              GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
              gsl_vector_complex_set(vdy, nmidx, val);
            }

          if (vz)
            {
              Z = 0.0;
              GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
              gsl_vector_complex_set(vz, nmidx, val);
            }

          if (vdz)
            {
              Z = 0.0;
              GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
              gsl_vector_complex_set(vdz, nmidx, val);
            }
        }
    }

  return s;
} /* poltor_row_tor() */

/*
poltor_nmidx()
  This function returns a unique index in [0,p-1] corresponding
to a given (n,m) pair. The array will look like:

[(1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) ...]

(the (0,0) coefficient is not solved for)

Things get a little more tricky when mmax != nmax, so the
base indices of each 'n' block are precomputed in _alloc and
stored for easy reference here. The offset of (n,m) in a given
'n' block is:

offset = m + min(n,mmax)

which defaults to the standard m + n for the case mmax = nmax

Inputs: type - POLTOR_IDX_xxx
        n    - SH degree (> 0)
        m    - SH order (-n <= m <= n)

Return: index in [0,nnm-1]
*/

size_t
poltor_nmidx(const size_t type, const size_t n, const int m, poltor_workspace *w)
{
  int ns;
  size_t mmax, *baseptr;
  size_t base; /* index of block for this n */
  size_t type_offset; /* offset for coefficient type (internal/external) */
  int offset;  /* offset within block for this m */
  size_t nmidx;

  if (type == POLTOR_IDX_PINT)
    {
      mmax = w->mmax_int;
      baseptr = w->base_int;
      type_offset = w->pint_offset;
    }
  else if (type == POLTOR_IDX_PEXT)
    {
      mmax = w->mmax_ext;
      baseptr = w->base_ext;
      type_offset = w->pext_offset;
    }
  else if (type == POLTOR_IDX_PSH)
    {
      mmax = w->mmax_sh;
      baseptr = w->base_sh;
      type_offset = w->psh_offset;
    }
  else if (type == POLTOR_IDX_TOR)
    {
      mmax = w->mmax_tor;
      baseptr = w->base_tor;
      type_offset = w->tor_offset;
    }
  else
    {
      fprintf(stderr, "poltor_nmidx: error: unknown coefficient type\n");
      return 0;
    }

  ns = (int) GSL_MIN(n, mmax);

  if (n == 0)
    {
      fprintf(stderr, "poltor_nmidx: error: n = 0\n");
      return 0;
    }
  else if (abs(m) > (int) mmax)
    {
      fprintf(stderr, "poltor_nmidx: error: m = %d\n", m);
      return 0;
    }

  base = baseptr[n]; /* precomputed */
  offset = m + ns;

  nmidx = base + offset + type_offset;

  return nmidx;
} /* poltor_nmidx() */

/*
poltor_jnmidx()
  Return index of shell poloidal coefficient corresponding
to (j,n,m) where j is the term in the Taylor series expansion
for q_{nm}(r)
*/

size_t
poltor_jnmidx(const size_t j, const size_t n, const int m, poltor_workspace *w)
{
  size_t nmidx = poltor_nmidx(POLTOR_IDX_PSH, n, m, w);

  return nmidx + j * w->nnm_sh;
} /* poltor_jnmidx() */

/*
poltor_regularize()
  Regularize least-squares system with Tikhonov matrix

Inputs: L         - (output) diag(L) damping matrix
        w         - workspace

Notes:
1) Following Olsen et al, 1997, L = diag(alpha * n^3) so that
alpha^2 n^6 is added to the diagonal of the LS matrix
*/

int
poltor_regularize(gsl_vector *L, poltor_workspace *w)
{
  int s = 0;
  const double alpha_int = w->alpha_int;
  const double alpha_sh = w->alpha_sh;
  const double alpha_tor = w->alpha_tor;
  const double beta = 1.0; /* power of n */
  const double gamma_int = 2000.0;
  const double gamma_sh = 100.0;
  size_t n, j;

  gsl_vector_set_zero(L);

  for (n = 1; n <= w->nmax_int; ++n)
    {
#if 1
      double nfac = pow((double) n, beta);
      double val = alpha_int * nfac;
#else
      double val = alpha_int * pow((double) n, beta) + gamma_int;
#endif
      int M = (int) GSL_MIN(n, w->mmax_int);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          gsl_vector_set(L, nmidx, val);
        }
    }

  {
    size_t nmax = w->nmax_sh;

    for (j = 0; j <= w->shell_J; ++j)
      {
        for (n = 1; n <= nmax; ++n)
          {
#if 1
            double nfac = pow((double) n, beta);
            double val = alpha_sh * nfac;
#else
            double val = alpha_sh * pow((double) n, beta) + gamma_sh;
#endif
            int M = (int) GSL_MIN(n, w->mmax_sh);
            int m;

            for (m = -M; m <= M; ++m)
              {
                size_t nmidx = poltor_jnmidx(j, n, m, w);
                gsl_vector_set(L, nmidx, val);
              }
          }
      }
  }

  for (n = 1; n <= w->nmax_tor; ++n)
    {
      double nfac = pow((double) n, beta);
      double val = alpha_tor * nfac;
      int M = (int) GSL_MIN(n, w->mmax_tor);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t nmidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
          gsl_vector_set(L, nmidx, val);
        }
    }

  return s;
} /* poltor_regularize() */

/*
poltor_eval_chi()
  Calculate toroidal current scalar function chi at given point
on the spherical shell r

chi = -Q r / mu_0 

Inputs: r          - radius of current shell (km)
        theta      - colatitude (radians)
        phi        - longitude (radians)
        idx        - POLTOR_IDX_xxx
        chi        - (output) current scalar in units of kA
        w          - workspace

Return: success/error
*/

static int
poltor_eval_chi(const double r, const double theta, const double phi,
                const size_t idx, double *chi, poltor_workspace *w)
{
  int s = 0;
  size_t n, m;
  complex double sum = 0.0;
  size_t nmax, mmax;

  if (idx == POLTOR_IDX_PINT)
    {
      nmax = w->nmax_int;
      mmax = w->mmax_int;
    }
  else if (idx == POLTOR_IDX_PSH)
    {
      nmax = w->nmax_sh;
      mmax = w->mmax_sh;
    }

  /* compute legendre functions */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                  cos(theta), w->Pnm, w->dPnm);

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= mmax; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);

      for (n = GSL_MAX(m, 1); n <= nmax; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;
        }
    }

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(idx, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          gsl_complex cnm = gsl_vector_complex_get(w->c, nmidx);
          complex double qnm = GSL_REAL(cnm) + I * GSL_IMAG(cnm);
          complex double Ynm;

          if (m >= 0)
            Ynm = w->Ynm[pidx];
          else
            Ynm = conj(w->Ynm[pidx]);

          sum += qnm * Ynm;
        }
    }

  /* current function = -Qr/mu0 (units of kA) */
  sum = -sum / POLTOR_MU_0 * r;
  *chi = creal(sum);

  return s;
} /* poltor_eval_chi() */
