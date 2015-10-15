/*
 * euler_calc.c
 *
 * This module contains routines for computing Euler angles
 * to rotate to a satellite-fixed basis defined by
 *
 * s1 = v / |v|
 * s2 = s3 x s1
 * s3 = -rhat
 *
 * This module has 2 modes, defined by EULER_CALC_CRF2S
 * If undefined, it calculates the Euler angles going from
 * VFM to s1,s2,s3.
 *
 * Alternatively, if defined, Euler angles are computed
 * to go from CRF to s1,s2,s3, by computing angles such that
 * R_s R(alpha,beta,gamma) =~ R_q
 * This bypasses the magnetometer completely, and no magnetic
 * data, nor any magnetic main field model is used in computing
 * the angles.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

#include "common.h"
#include "eci.h"
#include "ellipsoid.h"
#include "levmar.h"
#include "magdata.h"

#include "euler.h"
#include "euler_calc.h"

static int euler_calc_f(const gsl_vector *m, void *params, gsl_vector *f);
static int euler_calc_df(const gsl_vector *m, void *params, gsl_matrix *J);
static int euler_calc_vfm2sat(const size_t flags, const gsl_vector *m,
                              const double B_vfm[3], double B_sat[3],
                              const magdata *data);
static int euler_calc_vfm2nec(const size_t flags, const size_t idx, const gsl_vector *m,
                              double B_nec[3], const magdata *data);
static int euler_calc_vel(const size_t idx, double V[3], const magdata *data);
static int euler_calc_sph_basis(const size_t idx, double rhat[3], double that[3],
                                double phat[3], const magdata *data);
static int euler_calc_basis(const size_t idx, double s1_hat[3], double s2_hat[3],
                            double s3_hat[3], const magdata *data);
static int euler_calc_Rs(const size_t idx, gsl_matrix *Rs,
                         const magdata *data);
static void euler_calc_print_state (size_t iter, gsl_multifit_fdfsolver * s);

/*
 * Define to compute angles alpha,beta,gamma such that:
 *
 * R_s R(alpha,beta,gamma) =~ R_q
 *
 * so that R(alpha,beta,gamma) rotates from CRF to s1,s2,s3
 * system
 */
#define EULER_CALC_CRF2S         0

static double V_CRF[3] = { -1.0, 1.0, 2.0 };

typedef struct
{
  size_t n;
  const magdata *data;
  size_t flags;       /* EULER_FLG_xxx */
} euler_calc_params;

/*
euler_calc_proc()
  Compute Euler angles alpha,beta,gamma such that

R_q R_3(alpha,beta,gamma) B_VFM - B_main

is minimized in a least-squares sense.

Inputs: flags - EULER_FLG_xxx
        m     - (input/output)
                On input, initial guess of 3 angles in radians
                On output, vector of 3 Euler angles in radians
        data  - satellite data

Return: success/error
*/

int
euler_calc_proc(const size_t flags, gsl_vector *m, const magdata *data)
{
  int s = 0;
  const size_t p = 3;
  const size_t n = 3 * magdata_neuler(data); /* 3 field components */
  gsl_multifit_function_fdf f;
  euler_calc_params params;
  double fnorm0, fnorm;
  size_t niter;

  if (n < 500)
    {
      fprintf(stderr, "euler_calc_proc: insufficient number of data to process: %zu\n",
              n);
      return -1;
    }

  f.f = &euler_calc_f;
  f.df = &euler_calc_df;
  f.n = n;
  f.p = p;
  f.params = &params;

  params.n = n;
  params.data = data;
  params.flags = flags;

  fprintf(stderr, "euler_calc_proc: number of data = %zu\n", n);

#if 1

  /* use GSL */
  {
    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    const double xtol = 1.0e-8;
    const double gtol = 1.0e-8;
    gsl_multifit_fdfsolver *fdf_s = gsl_multifit_fdfsolver_alloc(T, n, p);
    int info;

    fprintf(stderr, "euler_calc_proc: initializing fdfsolver...");
    gsl_multifit_fdfsolver_set (fdf_s, &f, m);
    fprintf(stderr, "done\n");

    fnorm0 = gsl_blas_dnrm2(fdf_s->f);

    fprintf(stderr, "euler_calc_proc: computing Euler angles...");
    s = gsl_multifit_fdfsolver_driver (fdf_s, 100, xtol, gtol, 0.0, &info);
    if (s != GSL_SUCCESS)
      fprintf(stderr, "euler_calc_proc: error in fdfsolver: %d\n", s);
    fprintf(stderr, "done\n");

    niter = gsl_multifit_fdfsolver_niter(fdf_s),
    fnorm = gsl_blas_dnrm2(fdf_s->f);

    /* save angles */
    gsl_vector_memcpy(m, fdf_s->x);

    gsl_multifit_fdfsolver_free(fdf_s);
  }

#else

  /* use levmar */
  {
    levmar_workspace *lm_workspace_p = levmar_alloc(n, p);

    levmar_proc(m, &f, lm_workspace_p);

    fnorm0 = sqrt(levmar_chisq0(lm_workspace_p));
    fnorm = sqrt(levmar_chisq(lm_workspace_p));
    niter = levmar_niter(lm_workspace_p),

    levmar_free(lm_workspace_p);
  }

#endif

  fprintf(stderr,
          "iter: %3zu\n"
          "E = %15.8f %15.8f %15.8f [deg]\n"
          "initial |f(x)| = %g\n"
          "final   |f(x)| = %g\n",
          niter,
          gsl_vector_get (m, EULER_IDX_ALPHA) * 180.0 / M_PI,
          gsl_vector_get (m, EULER_IDX_BETA) * 180.0 / M_PI,
          gsl_vector_get (m, EULER_IDX_GAMMA) * 180.0 / M_PI,
          fnorm0,
          fnorm);

  return s;
} /* euler_calc_proc() */

static int
euler_calc_f(const gsl_vector *m, void *params, gsl_vector *f)
{
  int s = GSL_SUCCESS;
  euler_calc_params *p = ((euler_calc_params *) params);
  const magdata *data = p->data;
  size_t i;
  size_t idx = 0;

  gsl_vector_set_zero(f);

  for (i = 0; i < data->n; ++i)
    {
      double rhs[3], B_nec[3];

      /* skip flagged data */
      if (!MAGDATA_FitEuler(data->flags[i]))
        continue;

      s = euler_calc_vfm2nec(p->flags, i, m, B_nec, data);
      if (s)
        return s;

      rhs[0] = data->Bx_model[i];
      rhs[1] = data->By_model[i];
      rhs[2] = data->Bz_model[i];

      gsl_vector_set(f, idx++, B_nec[0] - rhs[0]);
      gsl_vector_set(f, idx++, B_nec[1] - rhs[1]);
      gsl_vector_set(f, idx++, B_nec[2] - rhs[2]);
    }

  assert(idx == p->n);

  return s;
} /* euler_calc_f() */

static int
euler_calc_df(const gsl_vector *m, void *params, gsl_matrix *J)
{
  euler_calc_params *p = ((euler_calc_params *) params);
  const magdata *data = p->data;
  size_t i;
  size_t idx = 0;

  gsl_matrix_set_zero(J);

  for (i = 0; i < data->n; ++i)
    {
      double B_nec[3];

      /* skip flagged data */
      if (!MAGDATA_FitEuler(data->flags[i]))
        continue;

      /* compute d/dalpha f_i */
      p->flags |= EULER_FLG_DERIV_ALPHA;
      euler_calc_vfm2nec(p->flags, i, m, B_nec, data);
      p->flags &= ~EULER_FLG_DERIV_ALPHA;

      gsl_matrix_set(J, idx, EULER_IDX_ALPHA, B_nec[0]);
      gsl_matrix_set(J, idx + 1, EULER_IDX_ALPHA, B_nec[1]);
      gsl_matrix_set(J, idx + 2, EULER_IDX_ALPHA, B_nec[2]);

      /* compute d/dbeta f_i */
      p->flags |= EULER_FLG_DERIV_BETA;
      euler_calc_vfm2nec(p->flags, i, m, B_nec, data);
      p->flags &= ~EULER_FLG_DERIV_BETA;

      gsl_matrix_set(J, idx, EULER_IDX_BETA, B_nec[0]);
      gsl_matrix_set(J, idx + 1, EULER_IDX_BETA, B_nec[1]);
      gsl_matrix_set(J, idx + 2, EULER_IDX_BETA, B_nec[2]);

      /* compute d/dgamma f_i */
      p->flags |= EULER_FLG_DERIV_GAMMA;
      euler_calc_vfm2nec(p->flags, i, m, B_nec, data);
      p->flags &= ~EULER_FLG_DERIV_GAMMA;

      gsl_matrix_set(J, idx, EULER_IDX_GAMMA, B_nec[0]);
      gsl_matrix_set(J, idx + 1, EULER_IDX_GAMMA, B_nec[1]);
      gsl_matrix_set(J, idx + 2, EULER_IDX_GAMMA, B_nec[2]);

      idx += 3;
    }

  assert(idx == p->n);

  return GSL_SUCCESS;
} /* euler_df() */

#if 0

/*
euler_calc_vfm2sat()
  Convert from VFM axes to satellite fixed frame.
These two coordinate axes are both fixed with respect to the satellite,
and are therefore not expected to change over the orbit. Therefore,
a standard Euler rotation with 3 Euler angles should bring them into
alignment.

The transformation applied is:

B_sat = R_x(alpha) R_y(beta) R_z(gamma) B_vfm

where alpha, beta, gamma are the Euler angles

Inputs: m     - euler parameters (Euler angles)
        B_vfm - input vector in VFM frame
        B_sat - (output) rotated vector in satellite frame
        w     - workspace

Notes: w->flags can control the output of this function. If
the EULER_FLG_DERIV_xxx flags are set, we take derivatives
of the corresponding rotation matrices when converting from
the magnetometer frame to satellite frame
*/

static int
euler_calc_vfm2sat(const size_t flags, const gsl_vector *m,
                 const double B_vfm[3], double B_sat[3],
                 const magdata *data)
{
  int s = 0;
  double alpha = gsl_vector_get(m, EULER_IDX_ALPHA),
         beta = gsl_vector_get(m, EULER_IDX_BETA),
         gamma = gsl_vector_get(m, EULER_IDX_GAMMA);

  /* apply 3D Euler rotation */
  euler_apply_deriv(flags, alpha, beta, gamma, B_vfm, B_sat);

  return s;
} /* euler_calc_vfm2sat() */

#endif

/*
euler_calc_vfm2nec()
  Convert a B vector measurement in the VFM frame to NEC

Inputs: flags - EULER_FLG_xxx, for derivative and Euler convention
                info
        idx   - measurement index to transform
        m     - euler parameters (for Euler angle rotation)
        B_nec - (output) geocentric magnetic field vector in
                north-east-down coordinates
                B_nec[0] = B_x (nT) (north)
                B_nec[1] = B_y (nT) (east)
                B_nec[2] = B_z (nT) (down)
        w     - euler workspace

Return: success or error

Notes: w->flags can control the output of this function. If
the EULER_FLG_DERIV_xxx flags are set, we take derivatives
of the corresponding rotation matrices when converting from
the magnetometer frame to satellite frame
*/

static int
euler_calc_vfm2nec(const size_t flags, const size_t idx, const gsl_vector *m,
                   double B_nec[3], const magdata *data)
{
  int s = 0;
  double B_vfm[3],  /* measured B in VFM frame */
         B_sat[3],  /* measured B in satellite-fixed frame */
         B_eci[3];  /* measured B in Cartesian ECI frame */
  double s1_hat[3], s2_hat[3], s3_hat[3];
  double rhat[3], that[3], phat[3];
  size_t j;

  B_vfm[0] = data->Bx_vfm[idx];
  B_vfm[1] = data->By_vfm[idx];
  B_vfm[2] = data->Bz_vfm[idx];

  if (flags & EULER_FLG_ROTSC)
    {
      /* use star camera quaternions for VFM->NEC transformation */

      double alpha = gsl_vector_get(m, EULER_IDX_ALPHA);
      double beta = gsl_vector_get(m, EULER_IDX_BETA);
      double gamma = gsl_vector_get(m, EULER_IDX_GAMMA);
      double *q = &(data->q[4*idx]);

      s = euler_vfm2nec(flags, alpha, beta, gamma, q, B_vfm, B_nec);
      return s;
    }

  fprintf(stderr, "UH OH\n");
  exit(1);

#if 0
  /* convert from VFM axes to satellite-fixed frame */
  euler_calc_vfm2sat(flags, m, B_vfm, B_sat, data);


  {
    double Rs_data[9];
    gsl_matrix_view Rs = gsl_matrix_view_array(Rs_data, 3, 3);
    gsl_vector_view sat = gsl_vector_view_array(B_sat, 3);
    gsl_vector_view nec = gsl_vector_view_array(B_nec, 3);

    /* construct R_s matrix explicitly */
    euler_calc_Rs(idx, &Rs.matrix, data);

    /* rotate: B_nec = R_s B_sat */
    gsl_blas_dgemv(CblasNoTrans, 1.0, &Rs.matrix, &sat.vector, 0.0, &nec.vector);
  }

#else

  /*
   * compute Cartesian ECI components of satellite frame basis vectors
   * at time t
   */
  euler_calc_basis(idx, s1_hat, s2_hat, s3_hat, data);

  /* convert B_sat to ECI frame */
  for (j = 0; j < 3; ++j)
    {
      B_eci[j] = B_sat[0] * s1_hat[j] +
                 B_sat[1] * s2_hat[j] +
                 B_sat[2] * s3_hat[j];
    }

  /* compute ECI components of spherical basis vectors */
  euler_calc_sph_basis(idx, rhat, that, phat, data);

  /* convert to local north-east-down coordinates */
  B_nec[0] = -vec_dot(B_eci, that); /* B_x = -B_t */
  B_nec[1] = vec_dot(B_eci, phat);  /* B_y = B_p */
  B_nec[2] = -vec_dot(B_eci, rhat); /* B_z = -B_r */

#endif

  return s;
} /* euler_calc_vfm2nec() */

/* compute ECI velocity using finite differences */
static int
euler_calc_vel(const size_t idx, double V[3], const magdata *data)
{
  int s = 0;
  double r0 = data->r[idx];
  double theta0 = data->theta[idx];
  double phi0 = data->phi[idx];
  double r1, theta1, phi1;
  double ecef0[3], ecef1[3];
  double eci0[3], eci1[3];
  time_t t0, t1, dt;
  size_t i;

  assert(idx < data->n - 1);

  t0 = satdata_epoch2timet(data->t[idx]);
  t1 = satdata_epoch2timet(data->t[idx + 1]);
  dt = t1 - t0;

  if (dt > 5)
    {
      fprintf(stderr, "euler_calc_vel: warning: dt = %ld [s]\n", dt);
    }

  r1 = data->r[idx + 1];
  theta1 = data->theta[idx + 1];
  phi1 = data->phi[idx + 1];

  ecef0[0] = r0 * sin(theta0) * cos(phi0);
  ecef0[1] = r0 * sin(theta0) * sin(phi0);
  ecef0[2] = r0 * cos(theta0);

  ecef1[0] = r1 * sin(theta1) * cos(phi1);
  ecef1[1] = r1 * sin(theta1) * sin(phi1);
  ecef1[2] = r1 * cos(theta1);

  ecef2eci_pos(t0, ecef0, eci0);
  ecef2eci_pos(t1, ecef1, eci1);

  for (i = 0; i < 3; ++i)
    V[i] = (eci1[i] - eci0[i]) / (double)dt;

  return s;
} /* euler_calc_vel() */

/*
euler_calc_sph_basis()
  Compute ECI components of spherical basis vectors for a given
time t
*/

static int
euler_calc_sph_basis(const size_t idx, double rhat[3], double that[3],
                   double phat[3], const magdata *data)
{
  int s = 0;
  double theta = data->theta[idx];
  double phi = data->phi[idx];
  double rhat_ecef[3], that_ecef[3], phat_ecef[3]; /* ECEF components */

  /* compute ECEF components of spherical unit vectors */
  sph_basis(theta, phi, rhat_ecef, that_ecef, phat_ecef);

  {
    time_t unix_time = satdata_epoch2timet(data->t[idx]);

    /* transform to ECI components */
    ecef2eci_pos(unix_time, rhat_ecef, rhat);
    ecef2eci_pos(unix_time, that_ecef, that);
    ecef2eci_pos(unix_time, phat_ecef, phat);
  }

  return s;
} /* euler_calc_sph_basis() */

/*
euler_calc_basis()
  Compute satellite-frame basis vectors in Cartesian (ECI)
coordinates, defined as:

s1 = v_t / |v_t|
s2 = s3 x s1
s3 = -rhat

where

v_t = v - (v . s3) s3

Inputs: idx    - index specifying satellite position
        s1_hat - (output) Cartesian ECI components of s1_hat
        s2_hat - (output) Cartesian ECI components of s2_hat
        s3_hat - (output) Cartesian ECI components of s3_hat
        data   - data
*/

static int
euler_calc_basis(const size_t idx, double s1_hat[3], double s2_hat[3],
               double s3_hat[3], const magdata *data)
{
  int s = 0;
  double V[3],     /* velocity vector in ECI */
         Vt[3];    /* tangent velocity vector in ECI */
  double V_3;      /* V . s3_hat */
  size_t j;
  time_t unix_time = satdata_epoch2timet(data->t[idx]);
  double theta = data->theta[idx];
  double phi = data->phi[idx];

  /* get ECI velocity */
  euler_calc_vel(idx, V, data);

#if 1

  /* s3 = -e_mu */
  {
    double r = data->r[idx];
    double e_mu[3], e_nu[3], e_phi[3]; /* ECEF components of ellipsoidal basis vectors */
    double e_mu_eci[3];                /* ECI components of e_mu */
    double r_ECEF[3];

    /* compute ECEF position */
    r_ECEF[0] = r * sin(theta) * cos(phi);
    r_ECEF[1] = r * sin(theta) * sin(phi);
    r_ECEF[2] = r * cos(theta);

    /* compute ellipsoidal basis vectors at time t */
    ellipsoid_basis_mu(r_ECEF, WGS84_MU, e_mu, e_nu, e_phi);

    /* transform e_mu from ECEF to ECI */
    ecef2eci_pos(unix_time, e_mu, e_mu_eci);

    /* compute s3 = -e_mu */
    for (j = 0; j < 3; ++j)
      s3_hat[j] = -e_mu_eci[j];
  }

#else

  /* s3 = -rhat */
  {
    double rhat_ecef[3], that_ecef[3], phat_ecef[3];
    double rhat_eci[3];

    /* compute spherical basis vectors in ECEF */
    sph_basis(theta, phi, rhat_ecef, that_ecef, phat_ecef);

    /* convert rhat to ECI */
    ecef2eci_pos(unix_time, rhat_ecef, rhat_eci);

    /* compute s3 = -rhat */
    for (j = 0; j < 3; ++j)
      s3_hat[j] = -rhat_eci[j];
  }

#endif

  /* compute V_3 component of V */
  V_3 = vec_dot(V, s3_hat);

  /* compute V_t = V - V_3 s3_hat */
  for (j = 0; j < 3; ++j)
    Vt[j] = V[j] - V_3 * s3_hat[j];

  /* s1 = V_t / |V_t| */
  for (j = 0; j < 3; ++j)
    s1_hat[j] = Vt[j] / vec_norm(Vt);

  /* s2 = s3 x s1 */
  sphcross(s3_hat, s1_hat, s2_hat);

  return s;
} /* euler_calc_basis() */

/*
euler_calc_Rs()
  Explicitely construct R_s rotation matrix
to take vector from s1,s2,s3 basis to NEC

Inputs: idx  - index specifying satellite position
        Rs   - (output) rotation matrix R_s
        data - satellite data
*/

static int
euler_calc_Rs(const size_t idx, gsl_matrix *Rs,
            const magdata *data)
{
  double s1_hat[3], s2_hat[3], s3_hat[3];
  double rhat_ecef[3], that_ecef[3], phat_ecef[3];
  double rhat_eci[3], that_eci[3], phat_eci[3];
  time_t unix_time = satdata_epoch2timet(data->t[idx]);
  double theta = data->theta[idx];
  double phi = data->phi[idx];

  /* compute s1,s2,s3 in ECI */
  euler_calc_basis(idx, s1_hat, s2_hat, s3_hat, data);

  /* compute spherical basis vectors in ECEF */
  sph_basis(theta, phi, rhat_ecef, that_ecef, phat_ecef);

  /* convert spherical basis vectors to ECI */
  ecef2eci_pos(unix_time, rhat_ecef, rhat_eci);
  ecef2eci_pos(unix_time, that_ecef, that_eci);
  ecef2eci_pos(unix_time, phat_ecef, phat_eci);

  gsl_matrix_set(Rs, 0, 0, -vec_dot(s1_hat, that_eci));
  gsl_matrix_set(Rs, 0, 1, -vec_dot(s2_hat, that_eci));
  gsl_matrix_set(Rs, 0, 2, -vec_dot(s3_hat, that_eci));

  gsl_matrix_set(Rs, 1, 0, vec_dot(s1_hat, phat_eci));
  gsl_matrix_set(Rs, 1, 1, vec_dot(s2_hat, phat_eci));
  gsl_matrix_set(Rs, 1, 2, vec_dot(s3_hat, phat_eci));

  gsl_matrix_set(Rs, 2, 0, -vec_dot(s1_hat, rhat_eci));
  gsl_matrix_set(Rs, 2, 1, -vec_dot(s2_hat, rhat_eci));
  gsl_matrix_set(Rs, 2, 2, -vec_dot(s3_hat, rhat_eci));

  return 0;
} /* euler_calc_Rs() */

static void
euler_calc_print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
}

int
euler_calc_print_residuals(char *filename, const size_t flags,
                           const magdata *data, const euler_workspace *w)
{
  int s = 0;
  FILE *fp;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "euler_calc_print_residuals: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: time (years)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: longitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: X NEC residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y NEC residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z NEC residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC X measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Y measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Z measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC X main field (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Y main field (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Z main field (nT)\n", i++);
  fprintf(fp, "# Field %zu: satellite direction (+/- 1)\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      double B_vfm[3], B_nec[3], B_main[3];
      double *q = &(data->q[4*i]);

      if (!MAGDATA_FitEuler(data->flags[i]))
        continue;

      B_vfm[0] = data->Bx_vfm[i];
      B_vfm[1] = data->By_vfm[i];
      B_vfm[2] = data->Bz_vfm[i];

      euler_vfm2nec_t(data->t[i], q, B_vfm, B_nec, w);

      B_main[0] = data->Bx_model[i];
      B_main[1] = data->By_model[i];
      B_main[2] = data->Bz_model[i];

      fprintf(fp, "%f %f %f %f %f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %d\n",
              satdata_epoch2year(data->t[i]),
              data->r[i] - R_EARTH_KM,
              90.0 - data->theta[i] * 180.0 / M_PI,
              data->phi[i] * 180.0 / M_PI,
              data->qdlat[i],
              B_nec[0] - B_main[0],
              B_nec[1] - B_main[1],
              B_nec[2] - B_main[2],
              B_nec[0],
              B_nec[1],
              B_nec[2],
              B_main[0],
              B_main[1],
              B_main[2],
              data->satdir[i]);
    }

  fclose(fp);

  return s;
}
