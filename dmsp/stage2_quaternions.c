/*
 * stage2_quaternions.c
 *
 * This module computes a correction to the quaternions to obtain
 * a spacecraft-fixed coordinate system, since it was observed with
 * F-15 that the satellite slowly rotations wrt the flight direction
 * during the orbit. The steps are:
 *
 * 1. For each point along the orbit, find rotation angle gamma_i which
 *    aligns B_VFM and a main field model in the NEC frame:
 *
 *      min || B_main - R_q R_3(gamma_i) B_VFM ||
 *
 *    This step uses an a priori model to define a NEC system to achieve
 *    the rotation from the VFM frame; this seems unavoidable since the
 *    satellite velocity vector is not fixed wrt to VFM axes. R_3 is
 *    a rotation about the VFM z axis (i.e. geodetic vertical):
 *
 *    R_3(gamma) = [ cos(gamma) -sin(gamma) 0 ]
 *                 [ sin(gamma)  cos(gamma) 0 ]
 *                 [     0           0      1 ]
 *
 * 2. Estimate a model gamma(theta) separately for ascending/descending tracks,
 *    where theta is colatitude:
 *
 *    gamma(theta) = c1 + c2*theta + c3*sin(theta)
 *
 * 3. Update quaternions from S/C frame to NEC using model:
 *
 *    R_q := R_q * R_3(gamma(theta))
 */

#include <gsl/gsl_test.h>

#include "euler.h"
#include "oct.h"
#include "quat.h"
#include "ellipsoid.h"

#if 1
struct quat_data
{
  const double *B;
  const double *E;
  const double *q;
};

static int
quat_f(const gsl_vector * x, void * params, gsl_vector * f)
{
  struct quat_data *qdata = (struct quat_data *) params;
  double gamma = gsl_vector_get(x, 0);
  double sg = sin(gamma);
  double cg = cos(gamma);
  const double *E = qdata->E;
  const double *B = qdata->B;
  double v[3];

  /* v = R_3(gamma) E */
  v[0] = E[0]*cg - E[1]*sg;
  v[1] = E[0]*sg + E[1]*cg;
  v[2] = E[2];

  /* v = R_q R_3(gamma) E */
  quat_apply(qdata->q, v, v);

  gsl_vector_set(f, 0, B[0] - v[0]);
  gsl_vector_set(f, 1, B[1] - v[1]);
  gsl_vector_set(f, 2, B[2] - v[2]);

  return GSL_SUCCESS;
}

static int
quat_df(const gsl_vector * x, void * params, gsl_matrix * J)
{
  struct quat_data *qdata = (struct quat_data *) params;
  double gamma = gsl_vector_get(x, 0);
  double sg = sin(gamma);
  double cg = cos(gamma);
  const double *E = qdata->E;
  double v[3];

  /* v = d/dgamma R_3(gamma) E */
  v[0] = -E[0]*sg - E[1]*cg;
  v[1] = E[0]*cg - E[1]*sg;
  v[2] = 0.0;

  /* v = R_q d/dgamma R_3(gamma) E */
  quat_apply(qdata->q, v, v);

  gsl_matrix_set(J, 0, 0, -v[0]);
  gsl_matrix_set(J, 1, 0, -v[1]);
  gsl_matrix_set(J, 2, 0, -v[2]);

  return GSL_SUCCESS;
}

static void
callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
#if 0
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  double rcond;

  gsl_multifit_nlinear_rcond(&rcond, w);

  fprintf(stderr,
          "iter %zu:\n"
          "\t |x|    = %.12e\n"
          "\t |f(x)| = %.12e\n"
          "\t cond(J) = %.12e\n",
          iter,
          gsl_blas_dnrm2(x),
          gsl_blas_dnrm2(f),
          1.0 / rcond);
#endif
}

/*
find_gamma2()

Solve:

min || B - R(gamma) E ||

where

R(gamma) = [ cos(gamma) -sin(gamma) ]
           [ sin(gamma)  cos(gamma) ]
*/

static double
find_gamma2(const double *q, const double B[3], const double E[3])
{
  const size_t n = 3;
  const size_t p = 1;
  const double xtol = 1e-12;
  const double gtol = 1e-12;
  const double ftol = 0.0;
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();
  gsl_vector *c = gsl_vector_calloc(p);
  double gamma;
  struct quat_data params;
  int info;
  gsl_vector *f;
  double chisq, chisq0;

  params.B = B;
  params.E = E;
  params.q = q;

  /* define the function to be minimized */
  fdf.f = quat_f;
  fdf.df = quat_df;
  fdf.fvv = NULL;     /* not using geodesic acceleration */
  fdf.n = n;
  fdf.p = p;
  fdf.params = &params;

  /* allocate workspace with default parameters */
  fdf_params.trs = gsl_multifit_nlinear_trs_lm;
  fdf_params.solver = gsl_multifit_nlinear_solver_svd;
  fdf_params.scale = gsl_multifit_nlinear_scale_levenberg;
  w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

  /* initialize solver with starting point and weights */
  gsl_multifit_nlinear_init (c, &fdf, w);

  f = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(f, f, &chisq0);

  /* solve the system with a maximum of 20 iterations */
  gsl_multifit_nlinear_driver(200, xtol, gtol, ftol, callback, &params, &info, w);

  gsl_blas_ddot(f, f, &chisq);
  gsl_vector_memcpy(c, w->x);

  gamma = gsl_vector_get(c, 0);

#if 0
  fprintf(stderr, "summary from method '%s/%s'\n",
          gsl_multifit_nlinear_name(w),
          gsl_multifit_nlinear_trs_name(w));
  fprintf(stderr, "number of iterations: %zu\n",
          gsl_multifit_nlinear_niter(w));
  fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
  fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
  fprintf(stderr, "reason for stopping: %s\n",
          (info == 1) ? "small step size" : "small gradient");
  fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
  fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));
#endif

  gsl_vector_free(c);
  gsl_multifit_nlinear_free(w);

  return gamma;
}

#endif

/*
find_gamma()
Find a rotation angle gamma which solves:

  min || H_main - R_q R(gamma) H_VFM ||

where:

  R(gamma) = [ cos(gamma) -sin(gamma) ]
             [ sin(gamma)  cos(gamma) ]

and the quaternion matrix R_q is assumed to be
a rotation about z, so R_q = R_q(1:2,1:2)

The solution is:

  min || H_main - M x ||

where:

  M = R_q * [ X_VFM -Y_VFM ]
            [ Y_VFM  X_VFM ]

  x = [ cos(gamma) ]
      [ sin(gamma) ]

Inputs: q      - quaternions
        H_VFM  - horizontal components of B_VFM
        H_main - horizontal components of B_main
*/

static double
find_gamma(const double *q, const double H_VFM[2], const double H_main[2])
{
  double Rq_data[9], Ainv_data[4], b_data[3], x_data[2];
  gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
  gsl_matrix_view Ainv = gsl_matrix_view_array(Ainv_data, 2, 2);
  gsl_vector_view b = gsl_vector_view_array(b_data, 2);
  gsl_vector_view x = gsl_vector_view_array(x_data, 2);
  double gamma;

  b_data[0] = H_main[0];
  b_data[1] = H_main[1];
  b_data[2] = 0.0;

  /* compute R_q */
  quat_q2R(q, &Rq.matrix);

  /*
   * compute: A^{-1} = 1 / (X_VFM^2 + Y_VFM^2) [  X_VFM Y_VFM ]
   *                                           [ -Y_VFM X_VFM ]
   */
  gsl_matrix_set(&Ainv.matrix, 0, 0, H_VFM[0]);
  gsl_matrix_set(&Ainv.matrix, 0, 1, H_VFM[1]);
  gsl_matrix_set(&Ainv.matrix, 1, 0, -H_VFM[1]);
  gsl_matrix_set(&Ainv.matrix, 1, 1, H_VFM[0]);
  gsl_matrix_scale(&Ainv.matrix, 1.0 / (H_VFM[0]*H_VFM[0] + H_VFM[1]*H_VFM[1]));

  /* compute b := R_q^T b */
  quat_apply_inverse(q, b_data, b_data);

  /* x = A^{-1} R_q^T b */
  gsl_blas_dgemv(CblasNoTrans, 1.0, &Ainv.matrix, &b.vector, 0.0, &x.vector);

  gamma = atan2(x_data[1], x_data[0]);

  return gamma;
}

static int
stage2_correct_quaternions(const char *filename, satdata_mag *data, track_workspace *track_p)
{
  size_t i, j;
  const size_t n = 10000; /* max data in a track */
  const size_t p = 3;     /* number of model parameters */
  gsl_vector *c = gsl_vector_alloc(p);
  gsl_matrix *cov = gsl_matrix_alloc(p, p);
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_rng *rng_p = gsl_rng_alloc(gsl_rng_default);
  gsl_multifit_linear_workspace *multifit_p = gsl_multifit_linear_alloc(n, p);
  double *c1 = gsl_vector_ptr(c, 0);
  double *c2 = gsl_vector_ptr(c, 1);
  double *c3 = gsl_vector_ptr(c, 2);
  FILE *fp = fopen(filename, "w");

  i = 1;
  fprintf(fp, "# Field %zu: time (decimal years)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: geocentric latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: gamma angle (degrees)\n", i++);
  fprintf(fp, "# Field %zu: modeled gamma angle (degrees)\n", i++);
  fprintf(fp, "# Field %zu: satellite flight direction (+1 north, -1 south)\n", i++);

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;
      size_t idx = 0;

      for (j = start_idx; j <= end_idx; ++j)
        {
          double theta = M_PI / 2.0 - data->latitude[j] * M_PI / 180.0;
          double B_model[4];
          double gamma;
          double *q = &(data->q[4*j]);

          /* calculate field model */
          satdata_mag_model(j, B_model, data);

#if 0
          gamma = find_gamma(q, &(data->B_VFM[3*j]), B_model);
#else
          gamma = find_gamma2(q, B_model, &(data->B_VFM[3*j]));
#endif

          if (fabs(data->qdlat[j]) <= 40.0)
            {
              gsl_matrix_set(X, idx, 0, 1.0);
              gsl_matrix_set(X, idx, 1, theta);
              gsl_matrix_set(X, idx, 2, sin(theta));

              gsl_vector_set(y, idx, gamma);

              ++idx;
            }
        }

      {
        gsl_matrix_view A = gsl_matrix_submatrix(X, 0, 0, idx, p);
        gsl_vector_view b = gsl_vector_subvector(y, 0, idx);
        double chisq;

        /* solve LS system */
        gsl_multifit_linear(&A.matrix, &b.vector, c, cov, &chisq, multifit_p);
      }

      /* now loop again through track and correct quaternions */

      for (j = start_idx; j <= end_idx; ++j)
        {
          double Rq_data[9], Rq_new_data[9], Rz_data[9];
          gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
          gsl_matrix_view Rq_new = gsl_matrix_view_array(Rq_new_data, 3, 3);
          gsl_matrix_view Rz = gsl_matrix_view_array(Rz_data, 3, 3);
          double theta = M_PI / 2.0 - data->latitude[j] * M_PI / 180.0;
          double *q = &(data->q[4*j]);
          double gamma_model = *c1 + *c2 * theta + *c3 * sin(theta);
          double cg = cos(gamma_model);
          double sg = sin(gamma_model);

          quat_q2R(q, &Rq.matrix);

          gsl_matrix_set(&Rz.matrix, 0, 0, cg);
          gsl_matrix_set(&Rz.matrix, 0, 1, -sg);
          gsl_matrix_set(&Rz.matrix, 0, 2, 0.0);

          gsl_matrix_set(&Rz.matrix, 1, 0, sg);
          gsl_matrix_set(&Rz.matrix, 1, 1, cg);
          gsl_matrix_set(&Rz.matrix, 1, 2, 0.0);

          gsl_matrix_set(&Rz.matrix, 2, 0, 0.0);
          gsl_matrix_set(&Rz.matrix, 2, 1, 0.0);
          gsl_matrix_set(&Rz.matrix, 2, 2, 1.0);

          /* R_q_new = R_q * R_z(model) */
          gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Rq.matrix, &Rz.matrix, 0.0, &Rq_new.matrix);

          if (j % 60 == 0)
            {
              double gamma;
              double B_model[4];

              satdata_mag_model(j, B_model, data);
              gamma = find_gamma(q, &(data->B_VFM[3*j]), B_model);

              fprintf(fp, "%f %f %f %f %f %d\n",
                      satdata_epoch2year(data->t[j]),
                      data->qdlat[j],
                      data->latitude[j],
                      gamma * 180.0 / M_PI,
                      gamma_model * 180.0 / M_PI,
                      tptr->satdir);
            }

          /* compute new quaternions */
          quat_R2q(&Rq_new.matrix, q);
        }
    }

  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_rng_free(rng_p);
  gsl_multifit_linear_free(multifit_p);

  fclose(fp);

  return GSL_SUCCESS;
}
