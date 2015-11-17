
typedef struct
{
  mfield_workspace *w;
} mfield_nonlinear_params;

static int mfield_init_nonlinear(mfield_workspace *w);
static int mfield_calc_f(const gsl_vector *x, void *params, gsl_vector *f);
static int mfield_calc_df(const gsl_vector *x, void *params, gsl_matrix *J);
static int mfield_calc_fdf(const int eval_J, const gsl_vector *x,
                           void *params, void *work);
static inline int mfield_jacobian_row(const double t, const size_t flags, gsl_vector * dX,
                                      const size_t extidx, const double dB_ext,
                                      const size_t euler_idx, const double B_nec_alpha,
                                      const double B_nec_beta, const double B_nec_gamma,
                                      gsl_vector *J, const mfield_workspace *w);
static inline int mfield_jacobian_row_F(const double t, gsl_vector * dX, gsl_vector * dY,
                                        gsl_vector * dZ, const double B_model[4],
                                        const size_t extidx, const double dB_ext[3],
                                        gsl_vector *J, const mfield_workspace *w);
static int mfield_nonlinear_matrices(gsl_matrix *dX, gsl_matrix *dY,
                                     gsl_matrix *dZ, mfield_workspace *w);
static double mfield_nonlinear_model_int(const double t, const gsl_vector *v,
                                         const gsl_vector *g, mfield_workspace *w);
static int mfield_nonlinear_model_ext(const double r, const double theta,
                                      const double phi, const gsl_vector *g,
                                      double dB[3], mfield_workspace *w);
static int mfield_nonlinear_histogram(const gsl_vector *c,
                                      mfield_workspace *w);
static int mfield_nonlinear_regularize(gsl_vector *diag,
                                       mfield_workspace *w);
static int mfield_nonlinear_driver (gsl_multifit_fdfridge * s,
                                    const size_t maxiter,
                                    const double xtol,
                                    const double gtol,
                                    const double ftol,
                                    int *info,
                                    mfield_workspace *w);
static int mfield_nonlinear_driver2 (const size_t maxiter,
                                     const double xtol, const double gtol,
                                     const double ftol, int *info,
                                     mfield_workspace *w);
void mfield_nonlinear_print_state(const size_t iter, gsl_multifit_fdfsolver *s,
                                  mfield_workspace *w);
static void mfield_nonlinear_print_state2(const size_t iter, mfield_workspace *w);

/*
mfield_calc_nonlinear()
  Solve linear least squares system, using previously stored
satellite data

Inputs: c - (input/output)
            on input, initial guess for coefficient vector
            on output, final coefficients
            units of nT, nT/year, nT/year^2
        w - workspace

Notes:
1) On input, w->data_workspace_p must be filled with satellite data
2a) mfield_init() must be called first to initialize various parameters,
    including weights
2b) this includes mfield_init_nonlinear(), called from mfield_init()
3) on output, coefficients are stored in w->c with following units:
4) static coefficients have units of nT
5) SV coefficients have units of nT/dimensionless_time
6) SA coefficients have units of nT/dimensionless_time^2
7) call mfield_coeffs() to convert coefficients to physical
   time units
*/

int
mfield_calc_nonlinear(gsl_vector *c, mfield_workspace *w)
{
  int s = 0;
  const size_t max_iter = 50;     /* maximum iterations */
  const double xtol = 1.0e-8;
  const double gtol = 1.0e-8;
  const double ftol = 1.0e-8;
  int info;
  const size_t p = w->p;          /* number of coefficients */
  const size_t n = w->nres;       /* number of residuals */
  gsl_multifit_function_fdf f;
  gsl_multilarge_function_fdf f2;
  gsl_vector *res_f = gsl_multifit_fdfridge_residual(w->fdf_s);
  struct timeval tv0, tv1;
  double res0;                    /* initial residual */
  FILE *fp_res;
  char resfile[2048];

  f.f = &mfield_calc_f;
  f.df = &mfield_calc_df;
  f.n = n;
  f.p = p;
  f.params = w;

  f2.fdf = mfield_calc_fdf;
  f2.p = p;
  f2.params = w;

  printv_octave(c, "c0");

  /* convert input vector from physical to dimensionless time units */
  mfield_coeffs(-1, c, c, w);

  /*
   * build and print residual histograms with previous coefficients
   * and previous wts_final vector
   */
  mfield_nonlinear_histogram(c, w);

  sprintf(resfile, "res.nlin.iter%zu.txt", w->niter);
  fp_res = fopen(resfile, "w");
  mfield_print_residuals(1, fp_res, NULL, NULL, NULL);

  /* compute robust weights with coefficients from previous iteration */
  {
    size_t i;

    /* compute f = Y_model - y_data with previous coefficients */
    mfield_calc_f(c, w, w->fvec);

    gsl_vector_memcpy(w->wfvec, w->fvec);

    /* compute weighted residuals r = sqrt(W) (Y - y) */
    for (i = 0; i < n; ++i)
      {
        double wi = gsl_vector_get(w->wts_spatial, i);
        double *ri = gsl_vector_ptr(w->wfvec, i);

        *ri *= sqrt(wi);
      }

    /* compute robust weights */
    gsl_multifit_robust_weights(w->wfvec, w->wts_final, w->robust_workspace_p);

    mfield_print_residuals(0, fp_res, w->fvec, w->wts_final, w->wts_spatial);

    /* compute final weights = wts_robust .* wts_spatial */
    gsl_vector_mul(w->wts_final, w->wts_spatial);
  }

  fprintf(stderr, "mfield_calc_nonlinear: wrote residuals to %s\n",
          resfile);
  fclose(fp_res);

#if MFIELD_SYNTH_DATA || MFIELD_NOWEIGHTS
  gsl_vector_set_all(w->wts_final, 1.0);
#endif

#if MFIELD_REGULARIZE && !MFIELD_SYNTH_DATA
  fprintf(stderr, "mfield_calc_nonlinear: regularizing least squares system...");
  mfield_nonlinear_regularize(w->lambda_diag, w);
  fprintf(stderr, "done\n");
#else
  gsl_vector_set_all(w->lambda_diag, 0.0);
#endif

  fprintf(stderr, "mfield_calc_nonlinear: initializing multilarge...");
  gettimeofday(&tv0, NULL);
  gsl_multilarge_nlinear_init(c, &f2, w->nlinear_workspace_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "mfield_calc_nonlinear: initializing fdfridge...");
  gettimeofday(&tv0, NULL);
  gsl_multifit_fdfridge_wset2(w->fdf_s, &f, c, w->lambda_diag, w->wts_final);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* compute initial residual */
  res0 = gsl_blas_dnrm2(res_f);

#if 1

  fprintf(stderr, "mfield_calc_nonlinear: computing nonlinear least squares solution...");
  gettimeofday(&tv0, NULL);
  s = mfield_nonlinear_driver2(max_iter, xtol, gtol, ftol, &info, w);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  if (s == GSL_SUCCESS)
    {
      fprintf(stderr, "mfield_calc_nonlinear: number of iterations: %zu\n",
              gsl_multilarge_nlinear_niter(w->nlinear_workspace_p));
      fprintf(stderr, "mfield_calc_nonlinear: function evaluations: %zu\n",
              f.nevalf);
      fprintf(stderr, "mfield_calc_nonlinear: Jacobian evaluations: %zu\n",
              f.nevaldf);
      fprintf(stderr, "mfield_calc_nonlinear: reason for stopping: %d\n", info);
      fprintf(stderr, "mfield_calc_nonlinear: initial residual: %.12e\n", res0);
      fprintf(stderr, "mfield_calc_nonlinear: final residual: %.12e\n",
              gsl_multilarge_nlinear_normf(w->nlinear_workspace_p));
    }
  else
    {
      fprintf(stderr, "mfield_calc_nonlinear: fdfridge failed: %s\n",
              gsl_strerror(s));
    }

  /* store final coefficients in physical units */
  {
    gsl_vector *x_final = gsl_multilarge_nlinear_position(w->nlinear_workspace_p);

    gsl_vector_memcpy(w->c, x_final);
    mfield_coeffs(1, w->c, c, w);

    printv_octave(c, "cfinal");
  }

#else
  fprintf(stderr, "mfield_calc_nonlinear: computing nonlinear least squares solution...");
  gettimeofday(&tv0, NULL);
  s = mfield_nonlinear_driver(w->fdf_s, max_iter, xtol, gtol, ftol, &info, w);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  if (s == GSL_SUCCESS)
    {
      fprintf(stderr, "mfield_calc_nonlinear: number of iterations: %zu\n",
              gsl_multifit_fdfridge_niter(w->fdf_s));
      fprintf(stderr, "mfield_calc_nonlinear: function evaluations: %zu\n",
              f.nevalf);
      fprintf(stderr, "mfield_calc_nonlinear: Jacobian evaluations: %zu\n",
              f.nevaldf);
      fprintf(stderr, "mfield_calc_nonlinear: reason for stopping: %d\n", info);
      fprintf(stderr, "mfield_calc_nonlinear: initial residual: %.12e\n", res0);
      fprintf(stderr, "mfield_calc_nonlinear: final residual: %.12e\n",
              gsl_blas_dnrm2(res_f));
    }
  else
    {
      fprintf(stderr, "mfield_calc_nonlinear: fdfridge failed: %s\n",
              gsl_strerror(s));
    }

  /* store final coefficients in physical units */
  {
    gsl_vector *x_final = gsl_multifit_fdfridge_position(w->fdf_s);

    gsl_vector_memcpy(w->c, x_final);
    mfield_coeffs(1, w->c, c, w);

    printv_octave(c, "cfinal");
  }
#endif

  w->niter++;

  return s;
} /* mfield_calc_nonlinear() */

/*
mfield_init_nonlinear()
  This function is called from mfield_init() to count
total number of residuals and allocate nonlinear least
squares workspaces

Notes:
1) weight_calc() must be called prior to this function to
compute spatial weights
*/

static int
mfield_init_nonlinear(mfield_workspace *w)
{
  int s = 0;
  const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmniel;
  const size_t p = w->p;        /* number of coefficients */
  const size_t nnm = w->nnm_mf; /* number of (n,m) harmonics */
  size_t ndata = 0;             /* number of distinct data points */
  size_t nres = 0;              /* number of residuals */
  struct timeval tv0, tv1;
  size_t i, j;

  /* count total number of residuals */
  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

      for (j = 0; j < mptr->n; ++j)
        {
          /* check if data point is discarded due to time interval */
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (mptr->flags[j] & MAGDATA_FLG_X)
            ++nres;
          if (mptr->flags[j] & MAGDATA_FLG_Y)
            ++nres;
          if (mptr->flags[j] & MAGDATA_FLG_Z)
            ++nres;

          /* don't increase nres if only fitting Euler angles */
          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            ++nres;

          ++ndata;
        }
    }

  fprintf(stderr, "mfield_init_nonlinear: %zu total data points\n", ndata);
  fprintf(stderr, "mfield_init_nonlinear: %zu total residuals\n", nres);
  fprintf(stderr, "mfield_init_nonlinear: %zu total parameters\n", p);

  /* precomputing these matrices make computing the residuals faster */
  w->mat_dX = gsl_matrix_alloc(ndata, nnm);
  w->mat_dY = gsl_matrix_alloc(ndata, nnm);
  w->mat_dZ = gsl_matrix_alloc(ndata, nnm);
  if (!w->mat_dX || !w->mat_dY || !w->mat_dZ)
    {
      GSL_ERROR("error allocating dX, dY, dZ", GSL_ENOMEM);
    }

  fprintf(stderr, "mfield_init_nonlinear: building matrices for nonlinear fit...");
  gettimeofday(&tv0, NULL);
  mfield_nonlinear_matrices(w->mat_dX, w->mat_dY, w->mat_dZ, w);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  w->nres = nres;
  w->ndata = ndata;

  w->fdf_s = gsl_multifit_fdfridge_alloc(T, nres, p);
  w->lambda_diag = gsl_vector_calloc(p);

  w->wts_spatial = gsl_vector_alloc(nres);
  w->wts_final = gsl_vector_alloc(nres);

  /* to save memory, allocate p-by-p workspace since a full n-by-p isn't needed */
  w->robust_workspace_p = gsl_multifit_robust_alloc(gsl_multifit_robust_huber, p, p);

  w->fvec = gsl_vector_alloc(nres);
  w->wfvec = gsl_vector_alloc(nres);

  gsl_vector_set_all(w->wts_final, 1.0);

  /* calculate spatial weights */
  {
    size_t idx = 0;
    size_t j;

    for (i = 0; i < w->nsat; ++i)
      {
        magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

        for (j = 0; j < mptr->n; ++j)
          {
            double wt; /* spatial weight */

            if (MAGDATA_Discarded(mptr->flags[j]))
              continue;

            track_weight_get(mptr->phi[j], mptr->theta[j], &wt, w->weight_workspace_p);

            if (mptr->flags[j] & MAGDATA_FLG_X)
              gsl_vector_set(w->wts_spatial, idx++, MFIELD_WEIGHT_X * wt);

            if (mptr->flags[j] & MAGDATA_FLG_Y)
              gsl_vector_set(w->wts_spatial, idx++, MFIELD_WEIGHT_Y * wt);

            if (mptr->flags[j] & MAGDATA_FLG_Z)
              gsl_vector_set(w->wts_spatial, idx++, MFIELD_WEIGHT_Z * wt);

            if (MAGDATA_ExistScalar(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
              gsl_vector_set(w->wts_spatial, idx++, MFIELD_WEIGHT_F * wt);
          }
      }
  }

  return s;
} /* mfield_init_nonlinear() */

/*
mfield_calc_f()
  Construct residual vector f using coefficients x

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
  gsl_matrix *dX = w->mat_dX;
  gsl_matrix *dY = w->mat_dY;
  gsl_matrix *dZ = w->mat_dZ;
  size_t i, j;
  size_t ridx = 0; /* index of residual */
  size_t didx = 0; /* index of data point */

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

      for (j = 0; j < mptr->n; ++j)
        {
          size_t k;
          double t = mptr->ts[j]; /* use scaled time */
          gsl_vector_view vx, vy, vz;
          double B_int[3];     /* internal field model */
          double B_extcorr[3]; /* external field correction model */
          double B_model[3];   /* a priori model (crustal/external) */
          double B_total[3];   /* internal + external */
          double B_obs[3];     /* observation vector NEC frame */

#if MFIELD_FIT_EXTFIELD
          size_t extidx = mfield_extidx(mptr->t[j], w);
          double extcoeff = gsl_vector_get(x, extidx);
          double dB_ext[3];
#endif

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          vx = gsl_matrix_row(dX, didx);
          vy = gsl_matrix_row(dY, didx);
          vz = gsl_matrix_row(dZ, didx);

          /* compute internal field model */
          B_int[0] = mfield_nonlinear_model_int(t, &vx.vector, x, w);
          B_int[1] = mfield_nonlinear_model_int(t, &vy.vector, x, w);
          B_int[2] = mfield_nonlinear_model_int(t, &vz.vector, x, w);

          /* load apriori model of external (and possibly crustal) field */
          B_model[0] = mptr->Bx_model[j];
          B_model[1] = mptr->By_model[j];
          B_model[2] = mptr->Bz_model[j];

#if MFIELD_FIT_EXTFIELD
          /* compute external field model correction */
          mfield_nonlinear_model_ext(mptr->r[j], mptr->theta[j], mptr->phi[j],
                                     x, dB_ext, w);

          /* add correction to POMME field */
          B_extcorr[0] = extcoeff * dB_ext[0];
          B_extcorr[1] = extcoeff * dB_ext[1];
          B_extcorr[2] = extcoeff * dB_ext[2];
#else
          B_extcorr[0] = B_extcorr[1] = B_extcorr[2] = 0.0;
#endif

          /* compute total modeled field (internal + external) */
          for (k = 0; k < 3; ++k)
            B_total[k] = B_int[k] + B_model[k] + B_extcorr[k];

#if MFIELD_FIT_EULER
          if (mptr->global_flags & MAGDATA_GLOBFLG_EULER)
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
            }
          else
#endif
            {
              /* use supplied NEC vector */
              B_obs[0] = mptr->Bx_nec[j];
              B_obs[1] = mptr->By_nec[j];
              B_obs[2] = mptr->Bz_nec[j];
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              /* set residual vector */
              if (f)
                gsl_vector_set(f, ridx, B_total[0] - B_obs[0]);

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              /* set residual vector */
              if (f)
                gsl_vector_set(f, ridx, B_total[1] - B_obs[1]);

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              /* set residual vector */
              if (f)
                gsl_vector_set(f, ridx, B_total[2] - B_obs[2]);
              else
                {
                  double wt = gsl_vector_get(w->wts_final, ridx);
                  wt = sqrt(wt);
                  wt = 1.0;
                  gsl_histogram_increment(w->hz, wt * (B_obs[2] - B_total[2]));
                }

              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double F = gsl_hypot3(B_total[0], B_total[1], B_total[2]);
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

          ++didx;
        } /* for (j = 0; j < mptr->n; ++j) */
    }

  if (f)
    assert(ridx == f->size);

  return s;
} /* mfield_calc_f() */

static int
mfield_calc_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
  int s = GSL_SUCCESS;
  mfield_workspace *w = (mfield_workspace *) params;
  gsl_matrix *dX = w->mat_dX;
  gsl_matrix *dY = w->mat_dY;
  gsl_matrix *dZ = w->mat_dZ;
  size_t i, j;
  size_t ridx = 0; /* index of residual */
  size_t didx = 0; /* index of data point */

  gsl_matrix_set_zero(J);

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

      for (j = 0; j < mptr->n; ++j)
        {
          double t = mptr->ts[j]; /* use scaled time */

#if MFIELD_FIT_EULER
          size_t euler_idx;
          double B_vfm[3], B_nec_alpha[3], B_nec_beta[3], B_nec_gamma[3];
#endif

#if MFIELD_FIT_EXTFIELD
          size_t extidx = mfield_extidx(mptr->t[j], w);
          double extcoeff = gsl_vector_get(x, extidx);
          double dB_ext[3];

          /* compute external field model */
          mfield_nonlinear_model_ext(mptr->r[j], mptr->theta[j], mptr->phi[j],
                                     x, dB_ext, w);
#endif

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

#if MFIELD_FIT_EULER
          /* compute Euler angle derivatives of B vector */
          if (mptr->global_flags & MAGDATA_GLOBFLG_EULER)
            {
              const double *q = &(mptr->q[4*j]);
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
            }
#endif

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv, vx;

                  /* main field portion */
                  vx = gsl_matrix_subrow(dX, didx, 0, w->nnm_mf);
                  Jv = gsl_matrix_subrow(J, ridx, 0, w->nnm_mf);
                  gsl_vector_memcpy(&Jv.vector, &vx.vector);

#if MFIELD_FIT_SECVAR
                  /* secular variation portion */
                  vx = gsl_matrix_subrow(dX, didx, 0, w->nnm_sv);
                  Jv = gsl_matrix_subrow(J, ridx, w->sv_offset, w->nnm_sv);
                  gsl_vector_memcpy(&Jv.vector, &vx.vector);
                  gsl_vector_scale(&Jv.vector, t);
#endif

#if MFIELD_FIT_SECACC
                  /* secular acceleration portion */
                  vx = gsl_matrix_subrow(dX, didx, 0, w->nnm_sa);
                  Jv = gsl_matrix_subrow(J, ridx, w->sa_offset, w->nnm_sa);
                  gsl_vector_memcpy(&Jv.vector, &vx.vector);
                  gsl_vector_scale(&Jv.vector, 0.5 * t * t);
#endif

#if MFIELD_FIT_EXTFIELD
                  gsl_matrix_set(J, ridx, extidx, dB_ext[0]);
#endif
                }

#if MFIELD_FIT_EULER
              /* check if fitting Euler angles to this data point */
              if (MAGDATA_FitEuler(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, euler_idx, -B_nec_alpha[0]);
                  gsl_matrix_set(J, ridx, euler_idx + 1, -B_nec_beta[0]);
                  gsl_matrix_set(J, ridx, euler_idx + 2, -B_nec_gamma[0]);
                }
#endif

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv, vy;

                  /* main field portion */
                  vy = gsl_matrix_subrow(dY, didx, 0, w->nnm_mf);
                  Jv = gsl_matrix_subrow(J, ridx, 0, w->nnm_mf);
                  gsl_vector_memcpy(&Jv.vector, &vy.vector);

#if MFIELD_FIT_SECVAR
                  /* secular variation portion */
                  vy = gsl_matrix_subrow(dY, didx, 0, w->nnm_sv);
                  Jv = gsl_matrix_subrow(J, ridx, w->sv_offset, w->nnm_sv);
                  gsl_vector_memcpy(&Jv.vector, &vy.vector);
                  gsl_vector_scale(&Jv.vector, t);
#endif

#if MFIELD_FIT_SECACC
                  /* secular acceleration portion */
                  vy = gsl_matrix_subrow(dY, didx, 0, w->nnm_sa);
                  Jv = gsl_matrix_subrow(J, ridx, w->sa_offset, w->nnm_sa);
                  gsl_vector_memcpy(&Jv.vector, &vy.vector);
                  gsl_vector_scale(&Jv.vector, 0.5 * t * t);
#endif

#if MFIELD_FIT_EXTFIELD
                  gsl_matrix_set(J, ridx, extidx, dB_ext[1]);
#endif
                }

#if MFIELD_FIT_EULER
              /* check if fitting Euler angles to this data point */
              if (MAGDATA_FitEuler(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, euler_idx, -B_nec_alpha[1]);
                  gsl_matrix_set(J, ridx, euler_idx + 1, -B_nec_beta[1]);
                  gsl_matrix_set(J, ridx, euler_idx + 2, -B_nec_gamma[1]);
                }
#endif

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv, vz;

                  /* main field portion */
                  vz = gsl_matrix_subrow(dZ, didx, 0, w->nnm_mf);
                  Jv = gsl_matrix_subrow(J, ridx, 0, w->nnm_mf);
                  gsl_vector_memcpy(&Jv.vector, &vz.vector);

#if MFIELD_FIT_SECVAR
                  /* secular variation portion */
                  vz = gsl_matrix_subrow(dZ, didx, 0, w->nnm_sv);
                  Jv = gsl_matrix_subrow(J, ridx, w->sv_offset, w->nnm_sv);
                  gsl_vector_memcpy(&Jv.vector, &vz.vector);
                  gsl_vector_scale(&Jv.vector, t);
#endif

#if MFIELD_FIT_SECACC
                  /* secular acceleration portion */
                  vz = gsl_matrix_subrow(dZ, didx, 0, w->nnm_sa);
                  Jv = gsl_matrix_subrow(J, ridx, w->sa_offset, w->nnm_sa);
                  gsl_vector_memcpy(&Jv.vector, &vz.vector);
                  gsl_vector_scale(&Jv.vector, 0.5 * t * t);
#endif

#if MFIELD_FIT_EXTFIELD
                  gsl_matrix_set(J, ridx, extidx, dB_ext[2]);
#endif
                }

#if MFIELD_FIT_EULER
              /* check if fitting Euler angles to this data point */
              if (MAGDATA_FitEuler(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, euler_idx, -B_nec_alpha[2]);
                  gsl_matrix_set(J, ridx, euler_idx + 1, -B_nec_beta[2]);
                  gsl_matrix_set(J, ridx, euler_idx + 2, -B_nec_gamma[2]);
                }
#endif

              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              gsl_vector_view Jv = gsl_matrix_row(J, ridx++);
              gsl_vector_view vx = gsl_matrix_row(dX, didx);
              gsl_vector_view vy = gsl_matrix_row(dY, didx);
              gsl_vector_view vz = gsl_matrix_row(dZ, didx);
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

          ++didx;
        }
    }

  assert(ridx == J->size1);

  return s;
} /* mfield_calc_df() */

static int
mfield_calc_fdf(const int eval_J, const gsl_vector *x, void *params, void *work)
{
  int s = GSL_SUCCESS;
  mfield_workspace *w = (mfield_workspace *) params;
  gsl_matrix *J = gsl_matrix_alloc(w->nres, w->p);
  gsl_vector *f = gsl_vector_alloc(w->nres);
  gsl_matrix *dX = w->mat_dX;
  gsl_matrix *dY = w->mat_dY;
  gsl_matrix *dZ = w->mat_dZ;
  double dB_ext[3] = { 0.0, 0.0, 0.0 };
  double B_extcorr[3] = { 0.0, 0.0, 0.0 }; /* external field correction model */
  size_t extidx = 0;
  double extcoeff = 0.0;
  size_t euler_idx = 0;
  double B_nec_alpha[3], B_nec_beta[3], B_nec_gamma[3];
  size_t i, j, k;
  size_t ridx = 0; /* index of residual */
  size_t didx = 0; /* index of data point */
  size_t rowidx = 0; /* row index in (J,f) */

  /* avoid unused variable warnings */
  (void)extcoeff;

  gsl_matrix_set_zero(J);

  /* loop over satellites */
  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

      /* loop over data for individual satellite */
      for (j = 0; j < mptr->n; ++j)
        {
          gsl_vector_view vx, vy, vz;
          double t = mptr->ts[j]; /* use scaled time */
          double B_int[3];        /* internal field model */
          double B_model[3];      /* a priori model (crustal/external) */
          double B_total[4];      /* internal + external */
          double B_obs[3];        /* observation vector NEC frame */
#if MFIELD_FIT_EULER
          double B_vfm[3];        /* observation vector VFM frame */
#endif

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          vx = gsl_matrix_row(dX, didx);
          vy = gsl_matrix_row(dY, didx);
          vz = gsl_matrix_row(dZ, didx);

          /* compute internal field model */
          B_int[0] = mfield_nonlinear_model_int(t, &vx.vector, x, w);
          B_int[1] = mfield_nonlinear_model_int(t, &vy.vector, x, w);
          B_int[2] = mfield_nonlinear_model_int(t, &vz.vector, x, w);

          /* load apriori model of external (and possibly crustal) field */
          B_model[0] = mptr->Bx_model[j];
          B_model[1] = mptr->By_model[j];
          B_model[2] = mptr->Bz_model[j];

#if MFIELD_FIT_EXTFIELD
          extidx = mfield_extidx(mptr->t[j], w);
          extcoeff = gsl_vector_get(x, extidx);

          /* compute external field model */
          mfield_nonlinear_model_ext(mptr->r[j], mptr->theta[j], mptr->phi[j],
                                     x, dB_ext, w);

          /* add correction to external field model */
          for (k = 0; k < 3; ++k)
            B_extcorr[k] = extcoeff * dB_ext[k];
#endif

          /* compute total modeled field (internal + external) */
          for (k = 0; k < 3; ++k)
            B_total[k] = B_int[k] + B_model[k] + B_extcorr[k];

#if MFIELD_FIT_EULER
          /* compute Euler angle derivatives of B vector */
          if (mptr->global_flags & MAGDATA_GLOBFLG_EULER)
            {
              const double *q = &(mptr->q[4*j]);
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

              /* compute observation vector in NEC frame */
              euler_vfm2nec(EULER_FLG_ZYX, alpha, beta, gamma, q, B_vfm, B_obs);
            }
          else
#endif
            {
              /* use supplied NEC vector */
              B_obs[0] = mptr->Bx_nec[j];
              B_obs[1] = mptr->By_nec[j];
              B_obs[2] = mptr->Bz_nec[j];
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              gsl_vector_view Jv = gsl_matrix_row(J, ridx);

              /* set residual vector X component */
              gsl_vector_set(f, ridx, B_total[0] - B_obs[0]);

              mfield_jacobian_row(t, mptr->flags[j], &vx.vector,
                                  extidx, dB_ext[0], euler_idx, B_nec_alpha[0],
                                  B_nec_beta[0], B_nec_gamma[0], &Jv.vector, w);

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              gsl_vector_view Jv = gsl_matrix_row(J, ridx);

              /* set residual vector Y component */
              gsl_vector_set(f, ridx, B_total[1] - B_obs[1]);

              mfield_jacobian_row(t, mptr->flags[j], &vy.vector,
                                  extidx, dB_ext[1], euler_idx, B_nec_alpha[1],
                                  B_nec_beta[1], B_nec_gamma[1], &Jv.vector, w);

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              gsl_vector_view Jv = gsl_matrix_row(J, ridx);

              /* set residual vector Z component */
              gsl_vector_set(f, ridx, B_total[2] - B_obs[2]);

              mfield_jacobian_row(t, mptr->flags[j], &vz.vector,
                                  extidx, dB_ext[2], euler_idx, B_nec_alpha[2],
                                  B_nec_beta[2], B_nec_gamma[2], &Jv.vector, w);

              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              gsl_vector_view Jv = gsl_matrix_row(J, ridx);
              double F_obs = mptr->F[j];

              B_total[3] = gsl_hypot3(B_total[0], B_total[1], B_total[2]);

              /* set scalar residual */
              gsl_vector_set(f, ridx, B_total[3] - F_obs);

              mfield_jacobian_row_F(t, &vx.vector, &vy.vector, &vz.vector,
                                    B_total, extidx, dB_ext, &Jv.vector, w);

              ++ridx;
            }

          ++didx;
        }
    }

  assert(ridx == w->nres);

  /* accumulate (J,f) into LS system */
  s = gsl_multilarge_nlinear_waccumulate(w->wts_final, J, f, work);

  gsl_matrix_free(J);
  gsl_vector_free(f);

  return s;
} /* mfield_calc_fdf() */

/*
mfield_jacobian_row()
  Construct a row of the Jacobian matrix corresponding to
a vector measurement

Inputs: t           - scaled timestamp
        flags       - MAGDATA_FLG_xxx flags for this data point
        dB_int      - Green's functions for desired vector component of
                      internal SH expansion, nnm_mf-by-1
        extidx      - index of external field coefficient
        dB_ext      - external field Green's function corresponding
                      to desired vector component
        euler_idx   - index of Euler angles
        B_nec_alpha -
        B_nec_beta  -
        B_nec_gamma -
        J           - (output) row of Jacobian, p-by-1
        w           - workspace
*/

static inline int
mfield_jacobian_row(const double t, const size_t flags, gsl_vector * dB_int,
                    const size_t extidx, const double dB_ext,
                    const size_t euler_idx, const double B_nec_alpha,
                    const double B_nec_beta, const double B_nec_gamma,
                    gsl_vector *J, const mfield_workspace *w)
{
  /* check if fitting MF to this data point */
  if (MAGDATA_FitMF(flags))
    {
      gsl_vector_view Jv, v;

      /* main field portion */
      Jv = gsl_vector_subvector(J, 0, w->nnm_mf);
      gsl_vector_memcpy(&Jv.vector, dB_int);

#if MFIELD_FIT_SECVAR
      /* secular variation portion */
      v = gsl_vector_subvector(dB_int, 0, w->nnm_sv);
      Jv = gsl_vector_subvector(J, w->sv_offset, w->nnm_sv);
      gsl_vector_memcpy(&Jv.vector, &v.vector);
      gsl_vector_scale(&Jv.vector, t);
#endif

#if MFIELD_FIT_SECACC
      /* secular acceleration portion */
      v = gsl_vector_subvector(dB_int, 0, w->nnm_sa);
      Jv = gsl_vector_subvector(J, w->sa_offset, w->nnm_sa);
      gsl_vector_memcpy(&Jv.vector, &v.vector);
      gsl_vector_scale(&Jv.vector, 0.5 * t * t);
#endif

#if MFIELD_FIT_EXTFIELD
      gsl_vector_set(J, extidx, dB_ext);
#endif
    }

#if MFIELD_FIT_EULER
  /* check if fitting Euler angles to this data point */
  if (MAGDATA_FitEuler(flags))
    {
      gsl_vector_set(J, euler_idx, -B_nec_alpha);
      gsl_vector_set(J, euler_idx + 1, -B_nec_beta);
      gsl_vector_set(J, euler_idx + 2, -B_nec_gamma);
    }
#endif

  return GSL_SUCCESS;
}

/*
mfield_jacobian_row_F()
  Construct a row of the Jacobian matrix corresponding to
a vector measurement

Inputs: t           - scaled timestamp
        flags       - MAGDATA_FLG_xxx flags for this data point
        dX          - Green's functions for X component
        dY          - Green's functions for Y component
        dZ          - Green's functions for Z component
        B_model     - total model vector
                      B_model[0] = X model
                      B_model[1] = Y model
                      B_model[2] = Z model
                      B_model[3] = F model
        extidx      - index of external field coefficient
        dB_ext      - external field vector Green's functions
        J           - (output) row of Jacobian, p-by-1
        w           - workspace
*/

static inline int
mfield_jacobian_row_F(const double t, gsl_vector * dX, gsl_vector * dY,
                      gsl_vector * dZ, const double B_model[4],
                      const size_t extidx, const double dB_ext[3],
                      gsl_vector *J, const mfield_workspace *w)
{
  size_t k;

  /* compute (X dX + Y dY + Z dZ) */
  for (k = 0; k < w->nnm_mf; ++k)
    {
      double dXk = gsl_vector_get(dX, k);
      double dYk = gsl_vector_get(dY, k);
      double dZk = gsl_vector_get(dZ, k);
      double val = B_model[0] * dXk +
                   B_model[1] * dYk +
                   B_model[2] * dZk;

      mfield_set_mf(J, k, val, w);
      mfield_set_sv(J, k, t * val, w);
      mfield_set_sa(J, k, 0.5 * t * t * val, w);
    }

#if MFIELD_FIT_EXTFIELD
  gsl_vector_set(J, extidx,
                 B_model[0] * dB_ext[0] +
                 B_model[1] * dB_ext[1] +
                 B_model[2] * dB_ext[2]);
#endif

  /* scale by 1/F */
  gsl_vector_scale(J, 1.0 / B_model[3]);

  return GSL_SUCCESS;
}

/*
mfield_nonlinear_matrices()
  Precompute matrices to evaluate residual vector and Jacobian
quickly in calc_f and calc_df

The ndata-by-nnm matrices computed are:

[dX] = dX/dg
[dY] = dY/dg
[dZ] = dZ/dg

Inputs: dX    - (output) dX/dg (main field)
        dY    - (output) dY/dg (main field)
        dZ    - (output) dZ/dg (main field)
*/

static int
mfield_nonlinear_matrices(gsl_matrix *dX, gsl_matrix *dY,
                          gsl_matrix *dZ, mfield_workspace *w)
{
  int s = GSL_SUCCESS;
  size_t i, j;
  size_t idx = 0;

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
      size_t n;
      int m;

      for (j = 0; j < mptr->n; ++j)
        {
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* compute basis functions for spherical harmonic expansions */
          mfield_green(r, theta, phi, w);

          for (n = 1; n <= w->nmax_mf; ++n)
            {
              int ni = (int) n;

              for (m = -ni; m <= ni; ++m)
                {
                  size_t cidx = mfield_coeff_nmidx(n, m);

                  gsl_matrix_set(dX, idx, cidx, w->dX[cidx]);
                  gsl_matrix_set(dY, idx, cidx, w->dY[cidx]);
                  gsl_matrix_set(dZ, idx, cidx, w->dZ[cidx]);
                }
            }

          ++idx;
        } /* for (j = 0; j < mptr->n; ++j) */
    } /* for (i = 0; i < w->nsat; ++i) */

  assert(idx == dX->size1);

  return s;
} /* mfield_nonlinear_matrices() */

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
                           const gsl_vector *g, mfield_workspace *w)
{
  gsl_vector_const_view gmf = gsl_vector_const_subvector(g, 0, w->nnm_mf);
  gsl_vector_const_view vmf = gsl_vector_const_subvector(v, 0, w->nnm_mf);
  double mf, sv = 0.0, sa = 0.0, val;

  /* compute v . x_mf */
  gsl_blas_ddot(&vmf.vector, &gmf.vector, &mf);

#if MFIELD_FIT_SECVAR
  {
    /* compute v . x_sv */
    gsl_vector_const_view gsv = gsl_vector_const_subvector(g, w->sv_offset, w->nnm_sv);
    gsl_vector_const_view vsv = gsl_vector_const_subvector(v, 0, w->nnm_sv);
    gsl_blas_ddot(&vsv.vector, &gsv.vector, &sv);
  }
#endif

#if MFIELD_FIT_SECACC
  {
    /* compute v . x_sa */
    gsl_vector_const_view gsa = gsl_vector_const_subvector(g, w->sa_offset, w->nnm_sa);
    gsl_vector_const_view vsa = gsl_vector_const_subvector(v, 0, w->nnm_sa);
    gsl_blas_ddot(&vsa.vector, &gsa.vector, &sa);
  }
#endif

  val = mf + t * sv + 0.5 * t * t * sa;

  return val;
}

/*
mfield_nonlinear_model_ext()
  Compute external field model:

r_k * (0.7 * external_dipole + 0.3 * internal_dipole)

where the dipole coefficients are aligned with the
direction of the main field dipole. r_k is the strength
of the ring current for a given day (k = doy) but is not
incorporated into the output of this function. This function
outputs the coefficient of r_k, ie, the "green's function" part.

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        g     - model coefficients
        dB    - (output) external field model green's functions
                dB[0] = X component of external field
                dB[1] = Y component of external field
                dB[2] = Z component of external field
        w     - workspace
*/

static int
mfield_nonlinear_model_ext(const double r, const double theta,
                           const double phi, const gsl_vector *g,
                           double dB[3], mfield_workspace *w)
{
#if MFIELD_FIT_EXTFIELD

  int s = 0;
  double g10 = gsl_vector_get(g, mfield_coeff_nmidx(1, 0));
  double g11 = gsl_vector_get(g, mfield_coeff_nmidx(1, 1));
  double h11 = gsl_vector_get(g, mfield_coeff_nmidx(1, -1));
  double g1 = gsl_hypot3(g10, g11, h11);
  double q[3];
  mfield_green_workspace *green_p = mfield_green_alloc(1, w->R);

  /* construct unit vector along internal dipole direction */
  q[mfield_coeff_nmidx(1, 0)] = g10 / g1;
  q[mfield_coeff_nmidx(1, 1)] = g11 / g1;
  q[mfield_coeff_nmidx(1, -1)] = h11 / g1;

  /* compute internal and external dipole field components */
  mfield_green_calc(r, theta, phi, green_p);
  mfield_green_ext(r, theta, phi, green_p);

  /* add external and induced sums */
  dB[0] = 0.7*vec_dot(q, green_p->dX_ext) + 0.3*vec_dot(q, green_p->dX);
  dB[1] = 0.7*vec_dot(q, green_p->dY_ext) + 0.3*vec_dot(q, green_p->dY);
  dB[2] = 0.7*vec_dot(q, green_p->dZ_ext) + 0.3*vec_dot(q, green_p->dZ);

  mfield_green_free(green_p);

  return s;

#else
  
  dB[0] = dB[1] = dB[2] = 0.0;
  return 0;

#endif

} /* mfield_nonlinear_model_ext() */

/*
mfield_nonlinear_histogram()
  Print residual histogram

Inputs: c - scaled/dimensionless coefficient vector
        w - workspace

Notes:
1) w->wts_final must be initialized prior to calling this function
*/

static int
mfield_nonlinear_histogram(const gsl_vector *c, mfield_workspace *w)
{
  int s = 0;
  FILE *fp;
  char filename[2048];

  /* reset histograms */
  gsl_histogram_reset(w->hf);
  gsl_histogram_reset(w->hz);

  /* loop through data and construct residual histograms */
  mfield_calc_f(c, w, NULL);

  /* scale histograms */
  gsl_histogram_scale(w->hf, 1.0 / gsl_histogram_sum(w->hf));
  gsl_histogram_scale(w->hz, 1.0 / gsl_histogram_sum(w->hz));

  /* print histograms to file */

  sprintf(filename, "reshistF.nlin.iter%zu.dat", w->niter);
  fprintf(stderr, "mfield_nonlinear_histogram: writing %s...", filename);
  fp = fopen(filename, "w");
  mfield_print_histogram(fp, w->hf);
  fclose(fp);
  fprintf(stderr, "done\n");

  sprintf(filename, "reshistZ.nlin.iter%zu.dat", w->niter);
  fprintf(stderr, "mfield_nonlinear_histogram: writing %s...", filename);
  fp = fopen(filename, "w");
  mfield_print_histogram(fp, w->hz);
  fclose(fp);
  fprintf(stderr, "done\n");

  return s;
} /* mfield_nonlinear_histogram() */

/*
mfield_nonlinear_regularize()
  Construct diag = diag(L) for regularized fit using frozen
flux assumption of Gubbins, 1983
*/

static int
mfield_nonlinear_regularize(gsl_vector *diag, mfield_workspace *w)
{
  int s = 0;
  const size_t nmin = 9;
  const size_t nmax = w->nmax_mf;
  const double c = 3485.0;       /* Earth core radius */
  const double a = MFIELD_RE_KM; /* Earth surface radius */
  const double ratio = a / c;
  size_t n;
  int m;
  double lambda_mf = w->lambda_mf;
  double lambda_sv = w->lambda_sv;
  double lambda_sa = w->lambda_sa;

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;
      double term = (n + 1.0) / sqrt(2.0*n + 1.0) * pow(ratio, n + 2.0);

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);

          mfield_set_mf(diag, cidx, lambda_mf, w);

          if (n >= nmin)
            {
              mfield_set_sv(diag, cidx, lambda_sv, w);
              mfield_set_sa(diag, cidx, lambda_sa * term, w);
            }
          else
            {
              mfield_set_sv(diag, cidx, lambda_mf, w);
              mfield_set_sa(diag, cidx, lambda_mf, w);
            }
        }
    }

  return s;
} /* mfield_nonlinear_regularize() */

/*
mfield_nonlinear_driver()
  Iterate the nonlinear least squares solver until completion

Inputs: s       - fdfridge workspace
        maxiter - maximum iterations to allow
        xtol    - tolerance in step x
        gtol    - tolerance in gradient
        ftol    - tolerance in ||f||
        info    - (output) info flag on why iteration terminated
                  1 = stopped due to small step size ||dx|
                  2 = stopped due to small gradient
                  3 = stopped due to small change in f
                  GSL_ETOLX = ||dx|| has converged to within machine
                              precision (and xtol is too small)
                  GSL_ETOLG = ||g||_inf is smaller than machine
                              precision (gtol is too small)
                  GSL_ETOLF = change in ||f|| is smaller than machine
                              precision (ftol is too small)

Return: GSL_SUCCESS if converged, GSL_MAXITER if maxiter exceeded without
converging
*/

static int
mfield_nonlinear_driver (gsl_multifit_fdfridge * s,
                         const size_t maxiter,
                         const double xtol, const double gtol,
                         const double ftol, int *info,
                         mfield_workspace *w)
{
  int status;
  gsl_multifit_fdfsolver *fdf_s = s->s;
  size_t iter = 0;

  do
    {
      if (iter % 5 == 0 || iter == 1)
        mfield_nonlinear_print_state(iter, fdf_s, w);

      status = gsl_multifit_fdfsolver_iterate (fdf_s);

      /*
       * if status is GSL_ENOPROG or GSL_SUCCESS, continue iterating,
       * otherwise the method has converged with a GSL_ETOLx flag
       */
      if (status != GSL_SUCCESS && status != GSL_ENOPROG)
        break;

      /* test for convergence */
      status = gsl_multifit_fdfsolver_test(fdf_s, xtol, gtol, ftol, info);
    }
  while (status == GSL_CONTINUE && ++iter < maxiter);

  /*
   * the following error codes mean that the solution has converged
   * to within machine precision, so record the error code in info
   * and return success
   */
  if (status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG)
    {
      *info = status;
      status = GSL_SUCCESS;
    }

  /* check if max iterations reached */
  if (iter >= maxiter && status != GSL_SUCCESS)
    status = GSL_EMAXITER;

  return status;
} /* mfield_nonlinear_driver() */

/*
mfield_nonlinear_driver2()
  Iterate the nonlinear least squares solver until completion

Inputs: maxiter - maximum iterations to allow
        xtol    - tolerance in step x
        gtol    - tolerance in gradient
        ftol    - tolerance in ||f||
        info    - (output) info flag on why iteration terminated
                  1 = stopped due to small step size ||dx|
                  2 = stopped due to small gradient
                  3 = stopped due to small change in f
                  GSL_ETOLX = ||dx|| has converged to within machine
                              precision (and xtol is too small)
                  GSL_ETOLG = ||g||_inf is smaller than machine
                              precision (gtol is too small)
                  GSL_ETOLF = change in ||f|| is smaller than machine
                              precision (ftol is too small)

Return: GSL_SUCCESS if converged, GSL_MAXITER if maxiter exceeded without
converging
*/

static int
mfield_nonlinear_driver2 (const size_t maxiter,
                          const double xtol, const double gtol,
                          const double ftol, int *info,
                          mfield_workspace *w)
{
  int status;
  size_t iter = 0;

  do
    {
      if (iter % 5 == 0 || iter == 1)
        mfield_nonlinear_print_state2(iter, w);

      status = gsl_multilarge_nlinear_iterate (w->nlinear_workspace_p);

      /*
       * if status is GSL_ENOPROG or GSL_SUCCESS, continue iterating,
       * otherwise the method has converged with a GSL_ETOLx flag
       */
      if (status != GSL_SUCCESS && status != GSL_ENOPROG)
        break;

      /* test for convergence */
      status = gsl_multilarge_nlinear_test(xtol, gtol, ftol, info,
                                           w->nlinear_workspace_p);
    }
  while (status == GSL_CONTINUE && ++iter < maxiter);

  /*
   * the following error codes mean that the solution has converged
   * to within machine precision, so record the error code in info
   * and return success
   */
  if (status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG)
    {
      *info = status;
      status = GSL_SUCCESS;
    }

  /* check if max iterations reached */
  if (iter >= maxiter && status != GSL_SUCCESS)
    status = GSL_EMAXITER;

  return status;
} /* mfield_nonlinear_driver() */

void
mfield_nonlinear_print_state(const size_t iter, gsl_multifit_fdfsolver *s,
                             mfield_workspace *w)
{
  fprintf(stderr, "iteration %zu:\n", iter);

  fprintf(stderr, "\t dipole: %12.4f %12.4f %12.4f [nT]\n",
          gsl_vector_get(s->x, mfield_coeff_nmidx(1, 0)),
          gsl_vector_get(s->x, mfield_coeff_nmidx(1, 1)),
          gsl_vector_get(s->x, mfield_coeff_nmidx(1, -1)));

#if MFIELD_FIT_EULER
  {
    size_t i;

    for (i = 0; i < w->nsat; ++i)
      {
        magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

        if (mptr->n == 0)
          continue;

        if (mptr->global_flags & MAGDATA_GLOBFLG_EULER)
          {
            double t0 = w->data_workspace_p->t0[i];
            size_t euler_idx = mfield_euler_idx(i, t0, w);

            fprintf(stderr, "\t euler : %12.4f %12.4f %12.4f [deg]\n",
                    gsl_vector_get(s->x, euler_idx) * 180.0 / M_PI,
                    gsl_vector_get(s->x, euler_idx + 1) * 180.0 / M_PI,
                    gsl_vector_get(s->x, euler_idx + 2) * 180.0 / M_PI);
          }
      }
  }
#endif

  fprintf(stderr, "\t |f(x)|: %12g\n", gsl_blas_dnrm2(s->f));
}

static void
mfield_nonlinear_print_state2(const size_t iter, mfield_workspace *w)
{
  gsl_vector *x = gsl_multilarge_nlinear_position(w->nlinear_workspace_p);
  double normf = gsl_multilarge_nlinear_normf(w->nlinear_workspace_p);

  fprintf(stderr, "iteration %zu:\n", iter);

  fprintf(stderr, "\t dipole: %12.4f %12.4f %12.4f [nT]\n",
          gsl_vector_get(x, mfield_coeff_nmidx(1, 0)),
          gsl_vector_get(x, mfield_coeff_nmidx(1, 1)),
          gsl_vector_get(x, mfield_coeff_nmidx(1, -1)));

#if MFIELD_FIT_EULER
  {
    size_t i;

    for (i = 0; i < w->nsat; ++i)
      {
        magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

        if (mptr->n == 0)
          continue;

        if (mptr->global_flags & MAGDATA_GLOBFLG_EULER)
          {
            double t0 = w->data_workspace_p->t0[i];
            size_t euler_idx = mfield_euler_idx(i, t0, w);

            fprintf(stderr, "\t euler : %12.4f %12.4f %12.4f [deg]\n",
                    gsl_vector_get(x, euler_idx) * 180.0 / M_PI,
                    gsl_vector_get(x, euler_idx + 1) * 180.0 / M_PI,
                    gsl_vector_get(x, euler_idx + 2) * 180.0 / M_PI);
          }
      }
  }
#endif

  fprintf(stderr, "\t ||f(x)||: %12g\n", normf);
}
