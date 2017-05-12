/*
 * mfield_multifit.c
 *
 * Contains routines for fitting magnetic field module using gsl_multifit_nlinear framework
 */

static int mfield_calc_f(const gsl_vector *x, void *params, gsl_vector *f);
static int mfield_calc_df(const gsl_vector *x, void *params, gsl_matrix *J);
static int mfield_nonlinear_model(const int res_flag, const gsl_vector * x, const magdata * mptr, const size_t idx,
                                  const size_t thread_id, double B_model[3], mfield_workspace *w);
static double mfield_nonlinear_model_int(const double t, const gsl_vector *v,
                                         const gsl_vector *g, const mfield_workspace *w);
static inline int jacobian_row_int(const double t, const gsl_vector * dB, gsl_vector * J, mfield_workspace * w);
static inline int jacobian_row_int_grad(const double t, const double t2, const gsl_vector * dB, const gsl_vector * dB2,
                                        gsl_vector * J, mfield_workspace * w);

/*
mfield_calc_f()
  Construct residual vector f(x) using OpenMP parallel
processing to compute all the Green's functions quickly.

Inputs: x      - model coefficients
        params - parameters
        f      - (output) residual vector
                 if set to NULL, the residual histograms
                 w->hf and w->hz are updated with residuals

Notes:
1) For the histograms, w->wts_final must be initialized prior
to calling this function
*/

static int
mfield_calc_f(const gsl_vector *x, void *params, gsl_vector *f)
{
  int s = GSL_SUCCESS;
  mfield_workspace *w = (mfield_workspace *) params;
  size_t i, j;
  struct timeval tv0, tv1;

#if DEBUG
  fprintf(stderr, "mfield_calc_f: entering function...\n");
#endif

  gettimeofday(&tv0, NULL);

  if (f)
    gsl_vector_set_zero(f);

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
      int fit_euler = w->params.fit_euler && (mptr->global_flags & MAGDATA_GLOBFLG_EULER);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          size_t ridx = mptr->index[j]; /* residual index for this data point */
          double B_model[3];            /* internal + external */
          double B_obs[3];              /* observation vector NEC frame */
          double B_model_ns[3];         /* N/S internal + external */
          double B_obs_ns[3];           /* N/S observation vector NEC frame */

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* compute vector model for this residual */
          mfield_nonlinear_model(0, x, mptr, j, thread_id, B_model, w);

          /* compute vector model for gradient residual (N/S or E/W) */
          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
            mfield_nonlinear_model(1, x, mptr, j, thread_id, B_model_ns, w);

          if (fit_euler)
            {
              /*
               * get the Euler angles for this satellite and time period,
               * and apply rotation
               */
              size_t euler_idx = mfield_euler_idx(i, mptr->t[j], w);
              double alpha = gsl_vector_get(x, euler_idx);
              double beta = gsl_vector_get(x, euler_idx + 1);
              double gamma = gsl_vector_get(x, euler_idx + 2);
              double *q = &(mptr->q[4*j]);
              double B_vfm[3];

              B_vfm[0] = mptr->Bx_vfm[j];
              B_vfm[1] = mptr->By_vfm[j];
              B_vfm[2] = mptr->Bz_vfm[j];

              /* rotate VFM vector to NEC */
              euler_vfm2nec(EULER_FLG_ZYX, alpha, beta, gamma, q, B_vfm, B_obs);

              if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                    MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
                {
                  double *q_ns = &(mptr->q_ns[4*j]);
                  double B_vfm_ns[3];

                  B_vfm_ns[0] = mptr->Bx_vfm_ns[j];
                  B_vfm_ns[1] = mptr->By_vfm_ns[j];
                  B_vfm_ns[2] = mptr->Bz_vfm_ns[j];

                  euler_vfm2nec(EULER_FLG_ZYX, alpha, beta, gamma, q_ns, B_vfm_ns, B_obs_ns);
                }
            }
          else
            {
              /* use supplied NEC vector */
              B_obs[0] = mptr->Bx_nec[j];
              B_obs[1] = mptr->By_nec[j];
              B_obs[2] = mptr->Bz_nec[j];

              if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                    MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
                {
                  B_obs_ns[0] = mptr->Bx_nec_ns[j];
                  B_obs_ns[1] = mptr->By_nec_ns[j];
                  B_obs_ns[2] = mptr->Bz_nec_ns[j];
                }
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              /* set residual vector */
              if (f)
                gsl_vector_set(f, ridx, B_model[0] - B_obs[0]);

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              /* set residual vector */
              if (f)
                gsl_vector_set(f, ridx, B_model[1] - B_obs[1]);

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              /* set residual vector */
              if (f)
                gsl_vector_set(f, ridx, B_model[2] - B_obs[2]);
              else
                {
                  double wt = gsl_vector_get(w->wts_final, ridx);
                  wt = sqrt(wt);
                  wt = 1.0;
                  gsl_histogram_increment(w->hz, wt * (B_obs[2] - B_model[2]));
                }

              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double F = gsl_hypot3(B_model[0], B_model[1], B_model[2]);
              double F_obs = mptr->F[j];

              if (f)
                gsl_vector_set(f, ridx, F - F_obs);
              else
                {
                  double wt = gsl_vector_get(w->wts_final, ridx);
                  wt = sqrt(wt);
                  wt = 1.0;
                  gsl_histogram_increment(w->hf, wt * (F_obs - F));
                }

              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            {
              /* set residual vector */
              if (f)
                gsl_vector_set(f, ridx, (B_model_ns[0] - B_model[0]) - (B_obs_ns[0] - B_obs[0]));

              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            {
              /* set residual vector */
              if (f)
                gsl_vector_set(f, ridx, (B_model_ns[1] - B_model[1]) - (B_obs_ns[1] - B_obs[1]));

              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            {
              /* set residual vector */
              if (f)
                gsl_vector_set(f, ridx, (B_model_ns[2] - B_model[2]) - (B_obs_ns[2] - B_obs[2]));

              ++ridx;
            }
        } /* for (j = 0; j < mptr->n; ++j) */
    }

  gettimeofday(&tv1, NULL);

#if DEBUG
  fprintf(stderr, "mfield_calc_f: leaving function (%g seconds)\n", time_diff(tv0, tv1));
#endif

#if 0
  if (f)
    {
      printv_octave(x, "x");
      printv_octave(f, "f");
      exit(1);
    }
#endif

  return s;
}

/*
mfield_calc_Wf()
  Compute weighted residuals:

f~(x) = sqrt(W) f(x)

Inputs: x      - model parameters
        params - parameters
        f      - (output) f~(x)
*/

static int
mfield_calc_Wf(const gsl_vector *x, void *params, gsl_vector *f)
{
  int s;
  mfield_workspace *w = (mfield_workspace *) params;
  size_t i;
  struct timeval tv0, tv1;

#if DEBUG
  fprintf(stderr, "mfield_calc_Wf: entering function...\n");
#endif
  gettimeofday(&tv0, NULL);

  s = mfield_calc_f(x, params, f);
  if (s)
    return s;

  for (i = 0; i < w->nres; ++i)
    {
      double wi = gsl_vector_get(w->wts_final, i);
      double *fi = gsl_vector_ptr(f, i);

      *fi *= sqrt(wi);
    }

  gettimeofday(&tv1, NULL);
#if DEBUG
  fprintf(stderr, "mfield_calc_Wf: leaving function (%g seconds)\n", time_diff(tv0, tv1));
#endif

#if 0
  printv_octave(f, "f");
  printv_octave(x, "x");
  printv_octave(w->wts_final, "wts");
  exit(1);
#endif

  return GSL_SUCCESS;
}

/*
mfield_calc_df()
  Compute Jacobian matrix J(x) using OpenMP to
calculate Green's functions quickly.

Inputs: x      - parameter vector, length p
        params - mfield workspace
        J      - (output) J(x), n-by-p
*/

static int
mfield_calc_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
  int s = GSL_SUCCESS;
  mfield_workspace *w = (mfield_workspace *) params;
  size_t i, j;
  struct timeval tv0, tv1;

#if DEBUG
  fprintf(stderr, "mfield_calc_df: entering function...\n");
#endif
  gettimeofday(&tv0, NULL);

  gsl_matrix_set_zero(J);

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
      int fit_euler = w->params.fit_euler && (mptr->global_flags & MAGDATA_GLOBFLG_EULER);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          double t = mptr->ts[j];       /* use scaled time */
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          size_t ridx = mptr->index[j]; /* residual index for this data point */

          /* internal Green's functions for current point */
          gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
          gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
          gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);

          /* internal Green's functions for gradient point */
          gsl_vector_view vx_ns = gsl_matrix_row(w->omp_dX_grad, thread_id);
          gsl_vector_view vy_ns = gsl_matrix_row(w->omp_dY_grad, thread_id);
          gsl_vector_view vz_ns = gsl_matrix_row(w->omp_dZ_grad, thread_id);

          size_t euler_idx;
          double B_vfm[3], B_nec_alpha[3], B_nec_beta[3], B_nec_gamma[3];
          double B_nec_alpha_ns[3], B_nec_beta_ns[3], B_nec_gamma_ns[3];

#if MFIELD_FIT_EXTFIELD
          size_t extidx = 0;
          double extcoeff = 0.0;
          double dB_ext[3];
#endif

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* calculate internal Green's functions */
          green_calc_int(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data, w->green_array_p[thread_id]);

          /* calculate internal Green's functions for gradient point (N/S or E/W) */
          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
            {
              green_calc_int(mptr->r_ns[j], mptr->theta_ns[j], mptr->phi_ns[j],
                             vx_ns.vector.data, vy_ns.vector.data, vz_ns.vector.data,
                             w->green_array_p[thread_id]);
            }

#if MFIELD_FIT_EXTFIELD
          /* external field is fitted only to data points used for main field modeling */
          if ((MAGDATA_ExistX(mptr->flags[j]) || MAGDATA_ExistY(mptr->flags[j]) ||
               MAGDATA_ExistZ(mptr->flags[j]) || MAGDATA_ExistScalar(mptr->flags[j])) &&
              (MAGDATA_FitMF(mptr->flags[j])))
            {
              extidx = mfield_extidx(mptr->t[j], w);
              extcoeff = gsl_vector_get(x, extidx);

              /* compute external field model */
              mfield_nonlinear_model_ext(mptr->r[j], mptr->theta[j], mptr->phi[j],
                                         x, dB_ext, w);
            }
#endif

          /* compute Euler angle derivatives of B vector */
          if (fit_euler)
            {
              double *q = &(mptr->q[4*j]);
              double alpha, beta, gamma;

              euler_idx = mfield_euler_idx(i, mptr->t[j], w);
              alpha = gsl_vector_get(x, euler_idx);
              beta = gsl_vector_get(x, euler_idx + 1);
              gamma = gsl_vector_get(x, euler_idx + 2);

              /* get vector in VFM frame */
              B_vfm[0] = mptr->Bx_vfm[j];
              B_vfm[1] = mptr->By_vfm[j];
              B_vfm[2] = mptr->Bz_vfm[j];

              /* compute alpha derivative of: R_q R_3 B_vfm */
              euler_vfm2nec(EULER_FLG_ZYX|EULER_FLG_DERIV_ALPHA, alpha, beta, gamma, q, B_vfm, B_nec_alpha);

              /* compute beta derivative of: R_q R_3 B_vfm */
              euler_vfm2nec(EULER_FLG_ZYX|EULER_FLG_DERIV_BETA, alpha, beta, gamma, q, B_vfm, B_nec_beta);

              /* compute gamma derivative of: R_q R_3 B_vfm */
              euler_vfm2nec(EULER_FLG_ZYX|EULER_FLG_DERIV_GAMMA, alpha, beta, gamma, q, B_vfm, B_nec_gamma);

              if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                    MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
                {
                  q = &(mptr->q_ns[4*j]);

                  /* get vector in VFM frame */
                  B_vfm[0] = mptr->Bx_vfm_ns[j];
                  B_vfm[1] = mptr->By_vfm_ns[j];
                  B_vfm[2] = mptr->Bz_vfm_ns[j];

                  /* compute alpha derivative of: R_q R_3 B_vfm_ns */
                  euler_vfm2nec(EULER_FLG_ZYX|EULER_FLG_DERIV_ALPHA, alpha, beta, gamma, q, B_vfm, B_nec_alpha_ns);

                  /* compute beta derivative of: R_q R_3 B_vfm_ns */
                  euler_vfm2nec(EULER_FLG_ZYX|EULER_FLG_DERIV_BETA, alpha, beta, gamma, q, B_vfm, B_nec_beta_ns);

                  /* compute gamma derivative of: R_q R_3 B_vfm_ns */
                  euler_vfm2nec(EULER_FLG_ZYX|EULER_FLG_DERIV_GAMMA, alpha, beta, gamma, q, B_vfm, B_nec_gamma_ns);
                }
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);

                  jacobian_row_int(t, &vx.vector, &Jv.vector, w);

#if MFIELD_FIT_EXTFIELD
                  gsl_matrix_set(J, ridx, extidx, dB_ext[0]);
#endif
                }

              /* check if fitting Euler angles to this data point */
              if (fit_euler && MAGDATA_FitEuler(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, euler_idx, -B_nec_alpha[0]);
                  gsl_matrix_set(J, ridx, euler_idx + 1, -B_nec_beta[0]);
                  gsl_matrix_set(J, ridx, euler_idx + 2, -B_nec_gamma[0]);
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);

                  jacobian_row_int(t, &vy.vector, &Jv.vector, w);

#if MFIELD_FIT_EXTFIELD
                  gsl_matrix_set(J, ridx, extidx, dB_ext[1]);
#endif
                }

              /* check if fitting Euler angles to this data point */
              if (fit_euler && MAGDATA_FitEuler(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, euler_idx, -B_nec_alpha[1]);
                  gsl_matrix_set(J, ridx, euler_idx + 1, -B_nec_beta[1]);
                  gsl_matrix_set(J, ridx, euler_idx + 2, -B_nec_gamma[1]);
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);

                  jacobian_row_int(t, &vz.vector, &Jv.vector, w);

#if MFIELD_FIT_EXTFIELD
                  gsl_matrix_set(J, ridx, extidx, dB_ext[2]);
#endif
                }

              /* check if fitting Euler angles to this data point */
              if (fit_euler && MAGDATA_FitEuler(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, euler_idx, -B_nec_alpha[2]);
                  gsl_matrix_set(J, ridx, euler_idx + 1, -B_nec_beta[2]);
                  gsl_matrix_set(J, ridx, euler_idx + 2, -B_nec_gamma[2]);
                }

              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              gsl_vector_view Jv = gsl_matrix_row(J, ridx++);
              double X, Y, Z, F;
              size_t k;

              /* compute internal X, Y, Z */
              X = mfield_nonlinear_model_int(t, &vx.vector, x, w);
              Y = mfield_nonlinear_model_int(t, &vy.vector, x, w);
              Z = mfield_nonlinear_model_int(t, &vz.vector, x, w);

              /* add apriori (external and crustal) field */
              X += mptr->Bx_model[j];
              Y += mptr->By_model[j];
              Z += mptr->Bz_model[j];

#if MFIELD_FIT_EXTFIELD
              /* add external field correction */
              X += extcoeff * dB_ext[0];
              Y += extcoeff * dB_ext[1];
              Z += extcoeff * dB_ext[2];
#endif

              F = gsl_hypot3(X, Y, Z);

              /* compute (X dX + Y dY + Z dZ) */
              for (k = 0; k < w->nnm_mf; ++k)
                {
                  double dXk = gsl_vector_get(&vx.vector, k);
                  double dYk = gsl_vector_get(&vy.vector, k);
                  double dZk = gsl_vector_get(&vz.vector, k);
                  double val = X * dXk + Y * dYk + Z * dZk;

                  mfield_set_mf(&Jv.vector, k, val, w);
                  mfield_set_sv(&Jv.vector, k, t * val, w);
                  mfield_set_sa(&Jv.vector, k, 0.5 * t * t * val, w);
                }

#if MFIELD_FIT_EXTFIELD
              gsl_vector_set(&Jv.vector, extidx, X * dB_ext[0] + Y * dB_ext[1] + Z * dB_ext[2]);
#endif

              /* scale by 1/F */
              gsl_vector_scale(&Jv.vector, 1.0 / F);
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);
                  jacobian_row_int_grad(t, mptr->ts_ns[j], &vx.vector, &vx_ns.vector, &Jv.vector, w);
                }

              /* check if fitting Euler angles to this data point */
              if (fit_euler && MAGDATA_FitEuler(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, euler_idx, B_nec_alpha[0] - B_nec_alpha_ns[0]);
                  gsl_matrix_set(J, ridx, euler_idx + 1, B_nec_beta[0] - B_nec_beta_ns[0]);
                  gsl_matrix_set(J, ridx, euler_idx + 2, B_nec_gamma[0] - B_nec_gamma_ns[0]);
                }

              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);
                  jacobian_row_int_grad(t, mptr->ts_ns[j], &vy.vector, &vy_ns.vector, &Jv.vector, w);
                }

              /* check if fitting Euler angles to this data point */
              if (fit_euler && MAGDATA_FitEuler(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, euler_idx, B_nec_alpha[1] - B_nec_alpha_ns[1]);
                  gsl_matrix_set(J, ridx, euler_idx + 1, B_nec_beta[1] - B_nec_beta_ns[1]);
                  gsl_matrix_set(J, ridx, euler_idx + 2, B_nec_gamma[1] - B_nec_gamma_ns[1]);
                }

              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);
                  jacobian_row_int_grad(t, mptr->ts_ns[j], &vz.vector, &vz_ns.vector, &Jv.vector, w);
                }

              /* check if fitting Euler angles to this data point */
              if (fit_euler && MAGDATA_FitEuler(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, euler_idx, B_nec_alpha[2] - B_nec_alpha_ns[2]);
                  gsl_matrix_set(J, ridx, euler_idx + 1, B_nec_beta[2] - B_nec_beta_ns[2]);
                  gsl_matrix_set(J, ridx, euler_idx + 2, B_nec_gamma[2] - B_nec_gamma_ns[2]);
                }

              ++ridx;
            }
        } /* for (j = 0; j < mptr->n; ++j) */
    } /* for (i = 0; i < w->nsat; ++i) */

  gettimeofday(&tv1, NULL);
#if DEBUG
  fprintf(stderr, "mfield_calc_df: leaving function (%g seconds)\n", time_diff(tv0, tv1));
#endif

#if 0
  print_octave(J, "J");
  exit(1);
#endif

  return s;
}

/*
mfield_nonlinear_model()
  Compute total model vector for a given residual

Inputs: res_flag  - 0 = normal residual
                    1 = gradient residual
        x         - parameter vector
        mptr      - magdata structure
        idx       - index of datum in magdata
        thread_id - OpenMP thread id
        B_model   - (output) B_model (X,Y,Z) in NEC
        w         - workspace
*/

static int
mfield_nonlinear_model(const int res_flag, const gsl_vector * x, const magdata * mptr, const size_t idx,
                       const size_t thread_id, double B_model[3], mfield_workspace *w)
{
  int s = 0;
  double t, ts, r, theta, phi;
  gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
  gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
  gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);
  double B_int[3], B_prior[3], B_extcorr[3];
  size_t k;

  if (res_flag == 0)
    {
      /* residual is for this data point specified by 'idx' */
      t = mptr->t[idx];
      ts = mptr->ts[idx];
      r = mptr->r[idx];
      theta = mptr->theta[idx];
      phi = mptr->phi[idx];

      /* load apriori model of external (and possibly crustal) field */
      B_prior[0] = mptr->Bx_model[idx];
      B_prior[1] = mptr->By_model[idx];
      B_prior[2] = mptr->Bz_model[idx];
    }
  else if (res_flag == 1)
    {
      /* residual is for gradient (N/S or E/W) */
      t = mptr->t_ns[idx];
      ts = mptr->ts_ns[idx];
      r = mptr->r_ns[idx];
      theta = mptr->theta_ns[idx];
      phi = mptr->phi_ns[idx];

      /* load apriori model of external (and possibly crustal) field */
      B_prior[0] = mptr->Bx_model_ns[idx];
      B_prior[1] = mptr->By_model_ns[idx];
      B_prior[2] = mptr->Bz_model_ns[idx];
    }

  green_calc_int(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data,
                 w->green_array_p[thread_id]);

  /* compute internal field model */
  B_int[0] = mfield_nonlinear_model_int(ts, &vx.vector, x, w);
  B_int[1] = mfield_nonlinear_model_int(ts, &vy.vector, x, w);
  B_int[2] = mfield_nonlinear_model_int(ts, &vz.vector, x, w);


#if MFIELD_FIT_EXTFIELD

  /* external field is fitted only to data points used for main field modeling */
  if ((MAGDATA_ExistX(mptr->flags[idx]) || MAGDATA_ExistY(mptr->flags[idx]) ||
       MAGDATA_ExistZ(mptr->flags[idx]) || MAGDATA_ExistScalar(mptr->flags[idx])) &&
      (MAGDATA_FitMF(mptr->flags[idx])))
    {
      size_t extidx = mfield_extidx(t, w);
      double extcoeff = gsl_vector_get(x, extidx);
      double dB_ext[3];

      /* compute external field model correction */
      mfield_nonlinear_model_ext(r, theta, phi, x, dB_ext, w);
    }

  /* add correction to POMME field */
  B_extcorr[0] = extcoeff * dB_ext[0];
  B_extcorr[1] = extcoeff * dB_ext[1];
  B_extcorr[2] = extcoeff * dB_ext[2];

#else

  B_extcorr[0] = B_extcorr[1] = B_extcorr[2] = 0.0;

#endif

  /* compute total modeled field (internal + external) */
  for (k = 0; k < 3; ++k)
    B_model[k] = B_int[k] + B_prior[k] + B_extcorr[k];

  return s;
}

/*
mfield_nonlinear_model_int()
  Evaluate internal field model for a given coefficient vector

Inputs: t - scaled time
        v - vector of basis functions (dX/dg,dY/dg,dZ/dg)
        g - model coefficients
        w - workspace

Return: model = v . g_mf + t*(v . g_sv) + 1/2*t^2*(v . g_sa)
*/

static double
mfield_nonlinear_model_int(const double t, const gsl_vector *v,
                           const gsl_vector *g, const mfield_workspace *w)
{
  gsl_vector_const_view gmf = gsl_vector_const_subvector(g, 0, w->nnm_mf);
  gsl_vector_const_view vmf = gsl_vector_const_subvector(v, 0, w->nnm_mf);
  double mf, sv = 0.0, sa = 0.0, val;

  /* compute v . x_mf */
  gsl_blas_ddot(&vmf.vector, &gmf.vector, &mf);

  if (w->nnm_sv > 0)
    {
      /* compute v . x_sv */
      gsl_vector_const_view gsv = gsl_vector_const_subvector(g, w->sv_offset, w->nnm_sv);
      gsl_vector_const_view vsv = gsl_vector_const_subvector(v, 0, w->nnm_sv);
      gsl_blas_ddot(&vsv.vector, &gsv.vector, &sv);
    }

  if (w->nnm_sa > 0)
    {
      /* compute v . x_sa */
      gsl_vector_const_view gsa = gsl_vector_const_subvector(g, w->sa_offset, w->nnm_sa);
      gsl_vector_const_view vsa = gsl_vector_const_subvector(v, 0, w->nnm_sa);
      gsl_blas_ddot(&vsa.vector, &gsa.vector, &sa);
    }

  val = mf + t * sv + 0.5 * t * t * sa;

  return val;
}

/*
jacobian_row_int()
  Fill in row of Jacobian corresponding to internal Green's functions

J(row,:) = [ dB, t * dB, 1/2 t^2 dB ]

accounting for differences in spherical harmonic degrees for MF, SV, and SA

Inputs: t  - scaled time of measurement
        dB - Green's functions for desired component (X,Y,Z), size nnm_mf
        J  - (output) Jacobian row; columns 1 - w->p_int will be filled in, size p_int
        w  - workspace
*/

static inline int
jacobian_row_int(const double t, const gsl_vector * dB, gsl_vector * J, mfield_workspace * w)
{
  int s = 0;
  const double fac_sa = 0.5 * t * t;
  size_t i;

  for (i = 0; i < w->nnm_mf; ++i)
    {
      double dBi = gsl_vector_get(dB, i);

      /* main field portion */
      gsl_vector_set(J, i, dBi);

      /* secular variation portion */
      if (i < w->nnm_sv)
        gsl_vector_set(J, w->sv_offset + i, t * dBi);

      if (i < w->nnm_sa)
        gsl_vector_set(J, w->sa_offset + i, fac_sa * dBi);
    }

  return s;
}

/*
jacobian_row_int_grad()
  Fill in row of Jacobian corresponding to internal Green's functions for
a gradient residual

J(row,:) = [ dB2 - dB, t2 * dB2 - t * dB, 1/2 t2^2 dB2 - 1/2 t^2 dB ]

accounting for differences in spherical harmonic degrees for MF, SV, and SA

Inputs: t   - scaled time of measurement
        t2  - scaled time of gradient measurement
        dB  - Green's functions for desired component (X,Y,Z), size nnm_mf
        dB2 - Green's functions for desired component (X,Y,Z) of gradient point (N/S or E/W), size nnm_mf
        J   - (output) Jacobian row; columns 1 - w->p_int will be filled in, size p_int
        w   - workspace
*/

static inline int
jacobian_row_int_grad(const double t, const double t2, const gsl_vector * dB, const gsl_vector * dB2,
                      gsl_vector * J, mfield_workspace * w)
{
  int s = 0;
  const double fac1 = 0.5 * t * t;
  const double fac2 = 0.5 * t2 * t2;
  size_t i;

  for (i = 0; i < w->nnm_mf; ++i)
    {
      double dBi = gsl_vector_get(dB, i);
      double dB2i = gsl_vector_get(dB2, i);

      /* main field portion */
      gsl_vector_set(J, i, dB2i - dBi);

      /* secular variation portion */
      if (i < w->nnm_sv)
        gsl_vector_set(J, w->sv_offset + i, t2 * dB2i - t * dBi);

      if (i < w->nnm_sa)
        gsl_vector_set(J, w->sa_offset + i, fac2 * dB2i - fac1 * dBi);
    }

  return s;
}
