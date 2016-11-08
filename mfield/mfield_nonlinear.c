#define OLD_FDF     0
#define DEBUG       0

typedef struct
{
  mfield_workspace *w;
} mfield_nonlinear_params;

static int mfield_init_nonlinear(mfield_workspace *w);
static int mfield_calc_f(const gsl_vector *x, void *params, gsl_vector *f);
static int mfield_calc_Wf(const gsl_vector *x, void *params, gsl_vector *f);
static int mfield_calc_df(const gsl_vector *x, void *params, gsl_matrix *J);
static int mfield_calc_df2(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u,
                           void *params, gsl_vector * v, gsl_matrix * JTJ);
static inline int mfield_jacobian_JTu(const double t, const size_t flags, const double weight,
                                      const gsl_vector * u, const size_t ridx, gsl_vector * dB_int, const size_t extidx,
                                      const double dB_ext, const size_t euler_idx, const double B_nec_alpha,
                                      const double B_nec_beta, const double B_nec_gamma,
                                      gsl_vector *JTy, const mfield_workspace *w);
static inline int mfield_jacobian_Ju(const double t, const size_t flags, const double weight,
                                     const gsl_vector * u, const size_t ridx, gsl_vector * dB_int, const size_t extidx,
                                     const double dB_ext, const size_t euler_idx, const double B_nec_alpha,
                                     const double B_nec_beta, const double B_nec_gamma,
                                     gsl_vector *Ju, const mfield_workspace *w);
static inline int mfield_jacobian_JTJ(const double t, const size_t flags, const double weight,
                                      gsl_vector * dB_int, const size_t extidx,
                                      const double dB_ext, const size_t euler_idx, const double B_nec_alpha,
                                      const double B_nec_beta, const double B_nec_gamma,
                                      gsl_matrix *JTJ, const mfield_workspace *w);
static inline int mfield_jacobian_row_F(CBLAS_TRANSPOSE_t TransJ, const double t, const double weight,
                                        const gsl_vector * u, const size_t ridx,
                                        gsl_vector * dX, gsl_vector * dY, gsl_vector * dZ,
                                        const double B_model[4], const size_t extidx, const double dB_ext[3],
                                        gsl_vector *J_int, gsl_matrix *JTJ, gsl_vector *v,
                                        const mfield_workspace *w);
static int mfield_nonlinear_vector_precompute(const gsl_vector *weights, mfield_workspace *w);
static int mfield_vector_green(const double t, const double weight, const gsl_vector *g,
                               gsl_vector *G, mfield_workspace *w);
static double mfield_nonlinear_model_int(const double t, const gsl_vector *v,
                                         const gsl_vector *g, const mfield_workspace *w);
static int mfield_nonlinear_model_ext(const double r, const double theta,
                                      const double phi, const gsl_vector *g,
                                      double dB[3], const mfield_workspace *w);
static int mfield_nonlinear_histogram(const gsl_vector *c,
                                      mfield_workspace *w);
static int mfield_nonlinear_regularize(gsl_vector *diag,
                                       mfield_workspace *w);
static void mfield_nonlinear_callback(const size_t iter, void *params,
                                      const gsl_multifit_nlinear_workspace *multifit_p);
static void mfield_nonlinear_callback2(const size_t iter, void *params,
                                       const gsl_multilarge_nlinear_workspace *multifit_p);

/*
mfield_calc_nonlinear()
  Solve linear least squares system, using previously stored
satellite data

Inputs: c    - (input/output)
               on input, initial guess for coefficient vector
               on output, final coefficients
               units of nT, nT/year, nT/year^2
        w    - workspace

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
  const double xtol = 1.0e-4;
  const double gtol = 1.5e-3;
  const double ftol = 1.0e-6;
  int info;
  const size_t p = w->p;          /* number of coefficients */
  const size_t n = w->nres;       /* number of residuals */
  gsl_multifit_nlinear_fdf fdf;
  gsl_multilarge_nlinear_fdf fdf2;
  gsl_vector *f;
  struct timeval tv0, tv1;
  double res0;                    /* initial residual */
  FILE *fp_res;
  char resfile[2048];

  fdf.f = mfield_calc_f;
  fdf.df = mfield_calc_df;
  fdf.fvv = NULL;
  fdf.n = n;
  fdf.p = p;
  fdf.params = w;

  fdf2.f = mfield_calc_Wf;
  fdf2.df = mfield_calc_df2;
  fdf2.fvv = NULL;
  fdf2.n = n;
  fdf2.p = p;
  fdf2.params = w;

  printv_octave(c, "c0");

  /* convert input vector from physical to dimensionless time units */
  mfield_coeffs(-1, c, c, w);

#if OLD_FDF
  /*
   * build and print residual histograms with previous coefficients
   * and previous wts_final vector
   */
  mfield_nonlinear_histogram(c, w);
#endif

  sprintf(resfile, "res.nlin.iter%zu.txt", w->niter);
  fp_res = fopen(resfile, "w");
  mfield_print_residuals(1, fp_res, NULL, NULL, NULL);

  /* compute robust weights with coefficients from previous iteration */
  if (w->niter > 0)
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
#if 0
      gsl_multifit_robust_weights(w->wfvec, w->wts_final, w->robust_workspace_p);
#else
      gsl_multifit_robust_weights(w->fvec, w->wts_final, w->robust_workspace_p);
#endif

      mfield_print_residuals(0, fp_res, w->fvec, w->wts_final, w->wts_spatial);

      /* compute final weights = wts_robust .* wts_spatial */
      gsl_vector_mul(w->wts_final, w->wts_spatial);
    }
  else
    {
      gsl_vector_memcpy(w->wts_final, w->wts_spatial);
    }

  fprintf(stderr, "mfield_calc_nonlinear: wrote residuals to %s\n",
          resfile);
  fclose(fp_res);

#if MFIELD_SYNTH_DATA || MFIELD_EMAG2 || MFIELD_NOWEIGHTS
  gsl_vector_set_all(w->wts_final, 1.0);
#endif

#if MFIELD_REGULARIZE && !MFIELD_SYNTH_DATA
  fprintf(stderr, "mfield_calc_nonlinear: regularizing least squares system...");

  /* compute diag(L) */
  mfield_nonlinear_regularize(w->lambda_diag, w);

  /* compute L^T L */
  gsl_vector_memcpy(w->LTL, w->lambda_diag);
  gsl_vector_mul(w->LTL, w->lambda_diag);

  fprintf(stderr, "done\n");
#else
  gsl_vector_set_all(w->lambda_diag, 0.0);
  gsl_vector_set_all(w->LTL, 0.0);
#endif

#if !OLD_FDF

  fprintf(stderr, "mfield_calc_nonlinear: precomputing vector J_int^T W J_int...");
  gettimeofday(&tv0, NULL);
  mfield_nonlinear_vector_precompute(w->wts_final, w);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "mfield_calc_nonlinear: initializing multilarge...");
  gettimeofday(&tv0, NULL);
  gsl_multilarge_nlinear_init(c, &fdf2, w->nlinear_workspace_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* compute initial residual */
  f = gsl_multilarge_nlinear_residual(w->nlinear_workspace_p);
  res0 = gsl_blas_dnrm2(f);

  fprintf(stderr, "mfield_calc_nonlinear: computing nonlinear least squares solution...");
  gettimeofday(&tv0, NULL);
  s = gsl_multilarge_nlinear_driver(max_iter, xtol, gtol, ftol,
                                    mfield_nonlinear_callback2, (void *) w,
                                    &info, w->nlinear_workspace_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  if (s == GSL_SUCCESS)
    {
      fprintf(stderr, "mfield_calc_nonlinear: NITER  = %zu\n",
              gsl_multilarge_nlinear_niter(w->nlinear_workspace_p));
      fprintf(stderr, "mfield_calc_nonlinear: NFEV   = %zu\n", fdf2.nevalf);
      fprintf(stderr, "mfield_calc_nonlinear: NJUEV  = %zu\n", fdf2.nevaldfu);
      fprintf(stderr, "mfield_calc_nonlinear: NJTJEV = %zu\n", fdf2.nevaldf2);
      fprintf(stderr, "mfield_calc_nonlinear: NAEV   = %zu\n", fdf2.nevalfvv);
      fprintf(stderr, "mfield_calc_nonlinear: reason for stopping: %d\n", info);
      fprintf(stderr, "mfield_calc_nonlinear: initial |f(x)|: %.12e\n", res0);
      fprintf(stderr, "mfield_calc_nonlinear: final   |f(x)|: %.12e\n",
              gsl_blas_dnrm2(f));
    }
  else
    {
      fprintf(stderr, "mfield_calc_nonlinear: failed: %s\n",
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

  fprintf(stderr, "mfield_calc_nonlinear: initializing multifit...");
  gettimeofday(&tv0, NULL);
  gsl_multifit_nlinear_winit(c, w->wts_final, &fdf, w->multifit_nlinear_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* compute initial residual */
  f = gsl_multifit_nlinear_residual(w->multifit_nlinear_p);
  res0 = gsl_blas_dnrm2(f);

  fprintf(stderr, "mfield_calc_nonlinear: computing nonlinear least squares solution...");
  gettimeofday(&tv0, NULL);
  s = gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,
                                  mfield_nonlinear_callback, (void *) w,
                                  &info, w->multifit_nlinear_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  if (s == GSL_SUCCESS)
    {
      fprintf(stderr, "mfield_calc_nonlinear: NITER = %zu\n",
              gsl_multifit_nlinear_niter(w->multifit_nlinear_p));
      fprintf(stderr, "mfield_calc_nonlinear: NFEV  = %zu\n", fdf.nevalf);
      fprintf(stderr, "mfield_calc_nonlinear: NJEV  = %zu\n", fdf.nevaldf);
      fprintf(stderr, "mfield_calc_nonlinear: NAEV  = %zu\n", fdf.nevalfvv);
      fprintf(stderr, "mfield_calc_nonlinear: reason for stopping: %d\n", info);
      fprintf(stderr, "mfield_calc_nonlinear: initial |f(x)|: %.12e\n", res0);
      fprintf(stderr, "mfield_calc_nonlinear: final   |f(x)|: %.12e\n",
              gsl_blas_dnrm2(f));
    }
  else
    {
      fprintf(stderr, "mfield_calc_nonlinear: multifit failed: %s\n",
              gsl_strerror(s));
    }

  /* store final coefficients in physical units */
  {
    gsl_vector *x_final = gsl_multifit_nlinear_position(w->multifit_nlinear_p);

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
  const size_t p = w->p;        /* number of coefficients */
  size_t ndata = 0;             /* number of distinct data points */
  size_t nres = 0;              /* total number of residuals */
  size_t nres_scal = 0;         /* number of scalar residuals */
  size_t nres_X = 0;            /* number of X residuals */
  size_t nres_Y = 0;            /* number of Y residuals */
  size_t nres_Z = 0;            /* number of Z residuals */
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

          /* store starting residual index for this data point */
          mptr->index[j] = nres_X + nres_Y + nres_Z + nres_scal;

          if (MAGDATA_ExistX(mptr->flags[j]))
            ++nres_X;
          if (MAGDATA_ExistY(mptr->flags[j]))
            ++nres_Y;
          if (MAGDATA_ExistZ(mptr->flags[j]))
            ++nres_Z;

          /* don't increase nres if only fitting Euler angles */
          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            ++nres_scal;

          ++ndata;
        }
    }

  nres = nres_scal + nres_X + nres_Y + nres_Z;

  w->nres_vec = nres_X + nres_Y + nres_Z;
  w->nres_scal = nres_scal;
  w->nres = nres;
  w->ndata = ndata;

  fprintf(stderr, "mfield_init_nonlinear: %zu distinct data points\n", ndata);
  fprintf(stderr, "mfield_init_nonlinear: %zu scalar residuals\n", w->nres_scal);
  fprintf(stderr, "mfield_init_nonlinear: %zu X residuals\n", nres_X);
  fprintf(stderr, "mfield_init_nonlinear: %zu Y residuals\n", nres_Y);
  fprintf(stderr, "mfield_init_nonlinear: %zu Z residuals\n", nres_Z);
  fprintf(stderr, "mfield_init_nonlinear: %zu total residuals\n", w->nres);
  fprintf(stderr, "mfield_init_nonlinear: %zu total parameters\n", p);

#if OLD_FDF

  /* allocate fit workspace */
  {
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_parameters fdf_params =
      gsl_multifit_nlinear_default_parameters();

    fdf_params.solver = gsl_multifit_nlinear_solver_cholesky;
    fdf_params.scale = gsl_multifit_nlinear_scale_levenberg;
    fdf_params.trs = gsl_multifit_nlinear_trs_dogleg;
    fdf_params.h_fvv = 0.5;
    w->multifit_nlinear_p = gsl_multifit_nlinear_alloc(T, &fdf_params, nres, p);
  }

#else

  /* allocate fit workspace */
  {
    const gsl_multilarge_nlinear_type *T = gsl_multilarge_nlinear_trust;
    gsl_multilarge_nlinear_parameters fdf_params =
      gsl_multilarge_nlinear_default_parameters();

    fdf_params.trs = gsl_multilarge_nlinear_trs_ddogleg;
    fdf_params.h_fvv = 0.5;
    w->nlinear_workspace_p = gsl_multilarge_nlinear_alloc(T, &fdf_params, nres, p);
  }

#endif

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

            if (MAGDATA_ExistX(mptr->flags[j]))
              gsl_vector_set(w->wts_spatial, idx++, MFIELD_WEIGHT_X * wt);

            if (MAGDATA_ExistY(mptr->flags[j]))
              gsl_vector_set(w->wts_spatial, idx++, MFIELD_WEIGHT_Y * wt);

            if (MAGDATA_ExistZ(mptr->flags[j]))
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

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          size_t k;
          double t = mptr->ts[j];       /* use scaled time */
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          size_t ridx = mptr->index[j]; /* residual index for this data point */
          gsl_vector_view vx, vy, vz;
          double B_int[3];              /* internal field model */
          double B_extcorr[3];          /* external field correction model */
          double B_model[3];            /* a priori model (crustal/external) */
          double B_total[3];            /* internal + external */
          double B_obs[3];              /* observation vector NEC frame */

#if MFIELD_FIT_EXTFIELD
          size_t extidx = 0;
          double extcoeff = 0.0;
          double dB_ext[3];
#endif

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          vx = gsl_matrix_row(w->omp_dX, thread_id);
          vy = gsl_matrix_row(w->omp_dY, thread_id);
          vz = gsl_matrix_row(w->omp_dZ, thread_id);

          green_calc_int(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data,
                         w->green_array_p[thread_id]);

          /* compute internal field model */
          B_int[0] = mfield_nonlinear_model_int(t, &vx.vector, x, w);
          B_int[1] = mfield_nonlinear_model_int(t, &vy.vector, x, w);
          B_int[2] = mfield_nonlinear_model_int(t, &vz.vector, x, w);

          /* load apriori model of external (and possibly crustal) field */
          B_model[0] = mptr->Bx_model[j];
          B_model[1] = mptr->By_model[j];
          B_model[2] = mptr->Bz_model[j];

#if MFIELD_FIT_EXTFIELD
          /* external field is fitted only to data points used for main field modeling */
          if ((MAGDATA_ExistX(mptr->flags[j]) || MAGDATA_ExistY(mptr->flags[j]) ||
               MAGDATA_ExistZ(mptr->flags[j]) || MAGDATA_ExistScalar(mptr->flags[j])) &&
              (MAGDATA_FitMF(mptr->flags[j])))
            {
              extidx = mfield_extidx(mptr->t[j], w);
              extcoeff = gsl_vector_get(x, extidx);

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
        } /* for (j = 0; j < mptr->n; ++j) */
    }

  gettimeofday(&tv1, NULL);
#if DEBUG
  fprintf(stderr, "mfield_calc_f: leaving function (%g seconds)\n", time_diff(tv0, tv1));
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

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          double t = mptr->ts[j];       /* use scaled time */
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          size_t ridx = mptr->index[j]; /* residual index for this data point */
          gsl_vector_view vx, vy, vz;

#if MFIELD_FIT_EULER
          size_t euler_idx;
          double B_vfm[3], B_nec_alpha[3], B_nec_beta[3], B_nec_gamma[3];
#endif

#if MFIELD_FIT_EXTFIELD
          size_t extidx = 0;
          double extcoeff = 0.0;
          double dB_ext[3];
#endif

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          vx = gsl_matrix_row(w->omp_dX, thread_id);
          vy = gsl_matrix_row(w->omp_dY, thread_id);
          vz = gsl_matrix_row(w->omp_dZ, thread_id);

          green_calc_int(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data,
                         w->green_array_p[thread_id]);

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
                  gsl_vector_view Jv, v;

                  /* main field portion */
                  Jv = gsl_matrix_subrow(J, ridx, 0, w->nnm_mf);
                  gsl_vector_memcpy(&Jv.vector, &vx.vector);

#if MFIELD_FIT_SECVAR
                  /* secular variation portion */
                  v = gsl_vector_subvector(&vx.vector, 0, w->nnm_sv);
                  Jv = gsl_matrix_subrow(J, ridx, w->sv_offset, w->nnm_sv);
                  gsl_blas_daxpy(t, &v.vector, &Jv.vector);
#endif

#if MFIELD_FIT_SECACC
                  /* secular acceleration portion */
                  v = gsl_vector_subvector(&vx.vector, 0, w->nnm_sa);
                  Jv = gsl_matrix_subrow(J, ridx, w->sa_offset, w->nnm_sa);
                  gsl_blas_daxpy(0.5 * t * t, &v.vector, &Jv.vector);
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
                  gsl_vector_view Jv, v;

                  /* main field portion */
                  Jv = gsl_matrix_subrow(J, ridx, 0, w->nnm_mf);
                  gsl_vector_memcpy(&Jv.vector, &vy.vector);

#if MFIELD_FIT_SECVAR
                  /* secular variation portion */
                  v = gsl_vector_subvector(&vy.vector, 0, w->nnm_sv);
                  Jv = gsl_matrix_subrow(J, ridx, w->sv_offset, w->nnm_sv);
                  gsl_blas_daxpy(t, &v.vector, &Jv.vector);
#endif

#if MFIELD_FIT_SECACC
                  /* secular acceleration portion */
                  v = gsl_vector_subvector(&vy.vector, 0, w->nnm_sa);
                  Jv = gsl_matrix_subrow(J, ridx, w->sa_offset, w->nnm_sa);
                  gsl_blas_daxpy(0.5 * t * t, &v.vector, &Jv.vector);
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
                  gsl_vector_view Jv, v;

                  /* main field portion */
                  Jv = gsl_matrix_subrow(J, ridx, 0, w->nnm_mf);
                  gsl_vector_memcpy(&Jv.vector, &vz.vector);

#if MFIELD_FIT_SECVAR
                  /* secular variation portion */
                  v = gsl_vector_subvector(&vz.vector, 0, w->nnm_sv);
                  Jv = gsl_matrix_subrow(J, ridx, w->sv_offset, w->nnm_sv);
                  gsl_blas_daxpy(t, &v.vector, &Jv.vector);
#endif

#if MFIELD_FIT_SECACC
                  /* secular acceleration portion */
                  v = gsl_vector_subvector(&vz.vector, 0, w->nnm_sa);
                  Jv = gsl_matrix_subrow(J, ridx, w->sa_offset, w->nnm_sa);
                  gsl_blas_daxpy(0.5 * t * t, &v.vector, &Jv.vector);
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
        }
    }

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
mfield_calc_df()
  Compute J^T J matrix and J^T u or J u vector using OpenMP
for speed improvement
*/

static int
mfield_calc_df2(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u,
                void *params, gsl_vector * v, gsl_matrix * JTJ)
{
  mfield_workspace *w = (mfield_workspace *) params;
  size_t i, j;
  gsl_matrix_view JTJ_int; /* internal field portion of J^T J */
  struct timeval tv0, tv1;

  gettimeofday(&tv0, NULL);
#if DEBUG
  fprintf(stderr, "mfield_calc_df2: entering function...\n");

  fprintf(stderr, "mfield_calc_df2: TransJ = %s...\n",
          TransJ == CblasTrans ? "trans" : "notrans");
#endif

  /* initialize outputs to 0 */
  if (v)
    gsl_vector_set_zero(v);

  if (JTJ)
    {
      gsl_matrix_set_zero(JTJ);

      /* copy previously computed vector internal field portion of J^T J
       * (doesn't depend on x) */
      JTJ_int = gsl_matrix_submatrix(JTJ, 0, 0, w->p_int, w->p_int);
      gsl_matrix_tricpy('L', 1, &JTJ_int.matrix, w->JTJ_vec);
    }

  /*
   * omp_rowidx[thread_id] contains the number of currently filled rows
   * of omp_J[thread_id]. When omp_J[thread_id] is filled, it is folded
   * into the JTJ_vec matrix and then omp_rowidx[thread_id] is reset to 0
   */
  for (i = 0; i < w->max_threads; ++i)
    w->omp_rowidx[i] = 0;

  /* loop over satellites */
  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

      /* loop over data for individual satellite */
#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          gsl_vector_view vx, vy, vz;
          size_t k;
          double t = mptr->ts[j];       /* use scaled time */
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          size_t ridx = mptr->index[j]; /* residual index for this data point */
          double B_int[3];              /* internal field model */
          double B_model[3];            /* a priori model (crustal/external) */
          double B_total[4];            /* internal + external */
          double dB_ext[3] = { 0.0, 0.0, 0.0 };
          double B_extcorr[3] = { 0.0, 0.0, 0.0 }; /* external field correction model */
          double B_nec_alpha[3], B_nec_beta[3], B_nec_gamma[3];
          size_t extidx = 0;
          size_t euler_idx = 0;
#if MFIELD_FIT_EXTFIELD
          double extcoeff = 0.0;
#endif
#if MFIELD_FIT_EULER
          double B_vfm[3];        /* observation vector VFM frame */
#endif

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          vx = gsl_matrix_row(w->omp_dX, thread_id);
          vy = gsl_matrix_row(w->omp_dY, thread_id);
          vz = gsl_matrix_row(w->omp_dZ, thread_id);

          green_calc_int(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data,
                         w->green_array_p[thread_id]);

          /* compute internal field model */
          B_int[0] = mfield_nonlinear_model_int(t, &vx.vector, x, w);
          B_int[1] = mfield_nonlinear_model_int(t, &vy.vector, x, w);
          B_int[2] = mfield_nonlinear_model_int(t, &vz.vector, x, w);

          /* load apriori model of external (and possibly crustal) field */
          B_model[0] = mptr->Bx_model[j];
          B_model[1] = mptr->By_model[j];
          B_model[2] = mptr->Bz_model[j];

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

              /* add correction to external field model */
              for (k = 0; k < 3; ++k)
                B_extcorr[k] = extcoeff * dB_ext[k];
            }
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
            }
#endif

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              double wj = gsl_vector_get(w->wts_final, ridx);

              if (TransJ == CblasTrans)
                {
#pragma omp critical
                  {
                    mfield_jacobian_JTu(t, mptr->flags[j], wj, u, ridx, &vx.vector,
                                        extidx, dB_ext[0], euler_idx, B_nec_alpha[0],
                                        B_nec_beta[0], B_nec_gamma[0], v, w);
                  }
                }
              else
                {
#pragma omp critical
                  {
                    mfield_jacobian_Ju(t, mptr->flags[j], wj, u, ridx, &vx.vector,
                                       extidx, dB_ext[0], euler_idx, B_nec_alpha[0],
                                       B_nec_beta[0], B_nec_gamma[0], v, w);
                  }
                }

#pragma omp critical
              {
                mfield_jacobian_JTJ(t, mptr->flags[j], wj, &vx.vector,
                                    extidx, dB_ext[0], euler_idx, B_nec_alpha[0],
                                    B_nec_beta[0], B_nec_gamma[0], JTJ, w);
              }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              double wj = gsl_vector_get(w->wts_final, ridx);

              if (TransJ == CblasTrans)
                {
#pragma omp critical
                  {
                    mfield_jacobian_JTu(t, mptr->flags[j], wj, u, ridx, &vy.vector,
                                        extidx, dB_ext[1], euler_idx, B_nec_alpha[1],
                                        B_nec_beta[1], B_nec_gamma[1], v, w);
                  }
                }
              else
                {
#pragma omp critical
                  {
                    mfield_jacobian_Ju(t, mptr->flags[j], wj, u, ridx, &vy.vector,
                                       extidx, dB_ext[1], euler_idx, B_nec_alpha[1],
                                       B_nec_beta[1], B_nec_gamma[1], v, w);
                  }
                }

#pragma omp critical
              {
                mfield_jacobian_JTJ(t, mptr->flags[j], wj, &vy.vector,
                                    extidx, dB_ext[1], euler_idx, B_nec_alpha[1],
                                    B_nec_beta[1], B_nec_gamma[1], JTJ, w);
              }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              double wj = gsl_vector_get(w->wts_final, ridx);

              if (TransJ == CblasTrans)
                {
#pragma omp critical
                  {
                    mfield_jacobian_JTu(t, mptr->flags[j], wj, u, ridx, &vz.vector,
                                        extidx, dB_ext[2], euler_idx, B_nec_alpha[2],
                                        B_nec_beta[2], B_nec_gamma[2], v, w);
                  }
                }
              else
                {
#pragma omp critical
                  {
                    mfield_jacobian_Ju(t, mptr->flags[j], wj, u, ridx, &vz.vector,
                                       extidx, dB_ext[2], euler_idx, B_nec_alpha[2],
                                       B_nec_beta[2], B_nec_gamma[2], v, w);
                  }
                }

#pragma omp critical
              {
                mfield_jacobian_JTJ(t, mptr->flags[j], wj, &vz.vector,
                                    extidx, dB_ext[2], euler_idx, B_nec_alpha[2],
                                    B_nec_beta[2], B_nec_gamma[2], JTJ, w);
              }

              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double wj = gsl_vector_get(w->wts_final, ridx);
              gsl_vector_view Jv = gsl_matrix_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id]++, 0, w->p_int);

              B_total[3] = gsl_hypot3(B_total[0], B_total[1], B_total[2]);

#pragma omp critical
              {
                mfield_jacobian_row_F(TransJ, t, wj, u, ridx, &vx.vector, &vy.vector, &vz.vector,
                                      B_total, extidx, dB_ext, &Jv.vector, JTJ, v, w);
              }

              ++ridx;
            }

          /* check if omp_J[thread_id] is full and should be folded into JTJ */
          if (w->omp_rowidx[thread_id] == 4 * w->data_block)
            {
              if (JTJ)
                {
                  /* accumulate scalar J_int^T J_int into J^T J; it is much faster to do this
                   * with blocks and dsyrk() rather than individual rows with dsyr() */
                  gsl_matrix_view Jm = gsl_matrix_submatrix(w->omp_J[thread_id], 0, 0, w->omp_rowidx[thread_id], w->p_int);

#pragma omp critical
                  {
                    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &Jm.matrix, 1.0, &JTJ_int.matrix);
                  }
                }

              /* reset for new block of rows */
              w->omp_rowidx[thread_id] = 0;
            }
        }
    }

  /* accumulate any last rows of internal field Green's functions */
  for (i = 0; i < w->max_threads; ++i)
    {
      if (JTJ && w->omp_rowidx[i] > 0)
        {
          gsl_matrix_view Jm = gsl_matrix_submatrix(w->omp_J[i], 0, 0, w->omp_rowidx[i], w->p_int);
          gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &Jm.matrix, 1.0, &JTJ_int.matrix);
        }
    }

#if 0
  /* regularize by adding L^T L to diag(J^T J) */
  if (JTJ)
    {
      gsl_vector_view v = gsl_matrix_diagonal(JTJ);
      gsl_vector_add(&v.vector, w->LTL);
    }
#endif

  gettimeofday(&tv1, NULL);
#if DEBUG
  fprintf(stderr, "mfield_calc_df2: leaving function... (%g seconds)\n", time_diff(tv0, tv1));
#endif

  return GSL_SUCCESS;
}

/*
mfield_jacobian_JTu()
  Update the J^T u vector with a new row of the Jacobian matrix,
corresponding to a vector residual.

Inputs: t           - scaled timestamp
        flags       - MAGDATA_FLG_xxx flags for this data point
        weight      - weight for this data point
        u           - input vector u, size n
        ridx        - residual index of this row in [0,nres-1]
        dB_int      - Green's functions for desired vector component of
                      internal SH expansion, nnm_mf-by-1
        extidx      - index of external field coefficient in [0,next-1]
        dB_ext      - external field Green's function corresponding
                      to desired vector component
        euler_idx   - index of Euler angles
        B_nec_alpha - Green's function for alpha Euler angle
        B_nec_beta  - Green's function for beta Euler angle
        B_nec_gamma - Green's function for gamma Euler angle
        JTu         - (output) J^T y vector
        w           - workspace
*/

static inline int
mfield_jacobian_JTu(const double t, const size_t flags, const double weight,
                    const gsl_vector * u, const size_t ridx, gsl_vector * dB_int, const size_t extidx,
                    const double dB_ext, const size_t euler_idx, const double B_nec_alpha,
                    const double B_nec_beta, const double B_nec_gamma,
                    gsl_vector *JTu, const mfield_workspace *w)
{
  const double y = gsl_vector_get(u, ridx);
  const double sWy = sqrt(weight) * y;
  gsl_vector_view g_mf = gsl_vector_subvector(dB_int, 0, w->nnm_mf);
#if MFIELD_FIT_SECVAR
  gsl_vector_view g_sv = gsl_vector_subvector(dB_int, 0, w->nnm_sv);
#endif
#if MFIELD_FIT_SECACC
  gsl_vector_view g_sa = gsl_vector_subvector(dB_int, 0, w->nnm_sa);
#endif

  (void) extidx;
  (void) dB_ext;
  (void) euler_idx;
  (void) B_nec_alpha;
  (void) B_nec_beta;
  (void) B_nec_gamma;

  /* check if fitting MF to this data point */
  if (MAGDATA_FitMF(flags))
    {
      gsl_vector_view v;

      /* update J^T y */
      v = gsl_vector_subvector(JTu, 0, w->nnm_mf);
      gsl_blas_daxpy(sWy, &g_mf.vector, &v.vector);

#if MFIELD_FIT_SECVAR
      v = gsl_vector_subvector(JTu, w->sv_offset, w->nnm_sv);
      gsl_blas_daxpy(t * sWy, &g_sv.vector, &v.vector);
#endif

#if MFIELD_FIT_SECACC
      v = gsl_vector_subvector(JTu, w->sa_offset, w->nnm_sa);
      gsl_blas_daxpy(0.5 * t * t * sWy, &g_sa.vector, &v.vector);
#endif

#if MFIELD_FIT_EXTFIELD
      {
        double *ptr = gsl_vector_ptr(JTu, extidx);

        /* update J^T y */
        *ptr += dB_ext * sWy;
      }
#endif /* MFIELD_FIT_EXTFIELD */
    }

#if MFIELD_FIT_EULER
  /* check if fitting Euler angles to this data point */
  if (MAGDATA_FitEuler(flags))
    {
      double x_data[3];
      gsl_vector_view vJTu = gsl_vector_subvector(JTu, euler_idx, 3);
      gsl_vector_view v = gsl_vector_view_array(x_data, 3);

      x_data[0] = -B_nec_alpha;
      x_data[1] = -B_nec_beta;
      x_data[2] = -B_nec_gamma;

      /* update J^T y */
      gsl_blas_daxpy(sWy, &v.vector, &vJTu.vector);
    }
#endif /* MFIELD_FIT_EULER */

  return GSL_SUCCESS;
}

/*
mfield_jacobian_Ju()
  Update the J u vector with a new row of the Jacobian matrix,
corresponding to a vector residual.

For this row of J u (specified by ridx), the value is given
by:

(J u)_i = sqrt(w_i) J_i^T u

where J_i^T is row i of the Jacobian (1-by-p)

Inputs: t           - scaled timestamp
        flags       - MAGDATA_FLG_xxx flags for this data point
        weight      - weight for this data point
        u           - input vector u, size p
        ridx        - residual index of this row in [0,nres-1]
        dB_int      - Green's functions for desired vector component of
                      internal SH expansion, nnm_mf-by-1
        extidx      - index of external field coefficient in [0,next-1]
        dB_ext      - external field Green's function corresponding
                      to desired vector component
        euler_idx   - index of Euler angles
        B_nec_alpha - Green's function for alpha Euler angle
        B_nec_beta  - Green's function for beta Euler angle
        B_nec_gamma - Green's function for gamma Euler angle
        Ju          - (output) J u vector
        w           - workspace
*/

static inline int
mfield_jacobian_Ju(const double t, const size_t flags, const double weight,
                   const gsl_vector * u, const size_t ridx, gsl_vector * dB_int, const size_t extidx,
                   const double dB_ext, const size_t euler_idx, const double B_nec_alpha,
                   const double B_nec_beta, const double B_nec_gamma,
                   gsl_vector *Ju, const mfield_workspace *w)
{
  double *Ju_ptr = gsl_vector_ptr(Ju, ridx);
  gsl_vector_view g_mf = gsl_vector_subvector(dB_int, 0, w->nnm_mf);
#if MFIELD_FIT_SECVAR
  gsl_vector_view g_sv = gsl_vector_subvector(dB_int, 0, w->nnm_sv);
#endif
#if MFIELD_FIT_SECACC
  gsl_vector_view g_sa = gsl_vector_subvector(dB_int, 0, w->nnm_sa);
#endif

  (void) extidx;
  (void) dB_ext;
  (void) euler_idx;
  (void) B_nec_alpha;
  (void) B_nec_beta;
  (void) B_nec_gamma;

  /* check if fitting MF to this data point */
  if (MAGDATA_FitMF(flags))
    {
      gsl_vector_const_view u_mf = gsl_vector_const_subvector(u, 0, w->nnm_mf);
      double tmp;

      /* update J u */
      gsl_blas_ddot(&g_mf.vector, &u_mf.vector, &tmp);
      *Ju_ptr = tmp;

#if MFIELD_FIT_SECVAR
      {
        gsl_vector_const_view u_sv = gsl_vector_const_subvector(u, w->sv_offset, w->nnm_sv);
        gsl_blas_ddot(&g_sv.vector, &u_sv.vector, &tmp);
        *Ju_ptr += t * tmp;
      }
#endif

#if MFIELD_FIT_SECACC
      {
        gsl_vector_const_view u_sa = gsl_vector_const_subvector(u, w->sa_offset, w->nnm_sa);
        gsl_blas_ddot(&g_sa.vector, &u_sa.vector, &tmp);
        *Ju_ptr += 0.5 * t * t * tmp;
      }
#endif

#if MFIELD_FIT_EXTFIELD

      /* update J u */
      *Ju_ptr += dB_ext * gsl_vector_get(u, extidx);

#endif /* MFIELD_FIT_EXTFIELD */
    }

#if MFIELD_FIT_EULER
  /* check if fitting Euler angles to this data point */
  if (MAGDATA_FitEuler(flags))
    {
      double x_data[3];
      gsl_vector_const_view vu = gsl_vector_const_subvector(u, euler_idx, 3);
      gsl_vector_view v = gsl_vector_view_array(x_data, 3);
      double tmp;

      x_data[0] = -B_nec_alpha;
      x_data[1] = -B_nec_beta;
      x_data[2] = -B_nec_gamma;

      gsl_blas_ddot(&vu.vector, &v.vector, &tmp);
      *Ju_ptr += tmp;
    }
#endif /* MFIELD_FIT_EULER */

  *Ju_ptr *= sqrt(weight);

  return GSL_SUCCESS;
}

/*
mfield_jacobian_JTJ()
  Update the J^T J matrix with a new row of the Jacobian matrix,
corresponding to a vector residual. The internal field portion of J^T J
does not need to be computed, since it is independent of the model parameters
and is pre-computed. Only the Euler and external field
portion of J^T J must be updated.

Inputs: t           - scaled timestamp
        flags       - MAGDATA_FLG_xxx flags for this data point
        weight      - weight for this data point
        u           - input vector u
        ridx        - residual index of this row in [0,nres-1]
        dB_int      - Green's functions for desired vector component of
                      internal SH expansion, nnm_mf-by-1
        extidx      - index of external field coefficient in [0,next-1]
        dB_ext      - external field Green's function corresponding
                      to desired vector component
        euler_idx   - index of Euler angles
        B_nec_alpha - Green's function for alpha Euler angle
        B_nec_beta  - Green's function for beta Euler angle
        B_nec_gamma - Green's function for gamma Euler angle
        JTJ         - (output) J^T J matrix, possibly NULL
        w           - workspace
*/

static inline int
mfield_jacobian_JTJ(const double t, const size_t flags, const double weight,
                    gsl_vector * dB_int, const size_t extidx,
                    const double dB_ext, const size_t euler_idx, const double B_nec_alpha,
                    const double B_nec_beta, const double B_nec_gamma,
                    gsl_matrix *JTJ, const mfield_workspace *w)
{
  gsl_vector_view g_mf = gsl_vector_subvector(dB_int, 0, w->nnm_mf);
#if MFIELD_FIT_SECVAR
  gsl_vector_view g_sv = gsl_vector_subvector(dB_int, 0, w->nnm_sv);
#endif
#if MFIELD_FIT_SECACC
  gsl_vector_view g_sa = gsl_vector_subvector(dB_int, 0, w->nnm_sa);
#endif

  (void) extidx;
  (void) dB_ext;
  (void) euler_idx;
  (void) B_nec_alpha;
  (void) B_nec_beta;
  (void) B_nec_gamma;

  /* check if fitting MF to this data point */
  if (MAGDATA_FitMF(flags))
    {
#if MFIELD_FIT_EXTFIELD
      /* update J^T J */
      gsl_vector_view v;
      double *ptr33 = gsl_matrix_ptr(JTJ, extidx, extidx);

      /* update (J^T J)_33 */
      *ptr33 += dB_ext * dB_ext * weight;

      /* update (J^T J)_31 = J_ext^T J_int */
      v = gsl_matrix_subrow(JTJ, extidx, 0, w->nnm_mf);
      gsl_blas_daxpy(dB_ext * weight, &g_mf.vector, &v.vector);

#if MFIELD_FIT_SECVAR
      v = gsl_matrix_subrow(JTJ, extidx, w->sv_offset, w->nnm_sv);
      gsl_blas_daxpy(t * dB_ext * weight, &g_sv.vector, &v.vector);
#endif

#if MFIELD_FIT_SECACC
      v = gsl_matrix_subrow(JTJ, extidx, w->sa_offset, w->nnm_sa);
      gsl_blas_daxpy(0.5 * t * t * dB_ext * weight, &g_sa.vector, &v.vector);
#endif

#endif /* MFIELD_FIT_EXTFIELD */
    }

#if MFIELD_FIT_EULER
  /* check if fitting Euler angles to this data point */
  if (MAGDATA_FitEuler(flags))
    {
      double x_data[3];
      gsl_vector_view v = gsl_vector_view_array(x_data, 3);
      gsl_matrix_view m;

      x_data[0] = -B_nec_alpha;
      x_data[1] = -B_nec_beta;
      x_data[2] = -B_nec_gamma;

      /* update (J^T J)_22 */
      m = gsl_matrix_submatrix(JTJ, euler_idx, euler_idx, 3, 3);
      gsl_blas_dsyr(CblasLower, weight, &v.vector, &m.matrix);

      if (MAGDATA_FitMF(flags))
        {
          /* update (J^T J)_21 */

          m = gsl_matrix_submatrix(JTJ, euler_idx, 0, 3, w->nnm_mf);
          gsl_blas_dger(weight, &v.vector, &g_mf.vector, &m.matrix);

#if MFIELD_FIT_SECVAR
          m = gsl_matrix_submatrix(JTJ, euler_idx, w->sv_offset, 3, w->nnm_sv);
          gsl_blas_dger(t * weight, &v.vector, &g_sv.vector, &m.matrix);
#endif

#if MFIELD_FIT_SECACC
          m = gsl_matrix_submatrix(JTJ, euler_idx, w->sa_offset, 3, w->nnm_sa);
          gsl_blas_dger(0.5 * t * t * weight, &v.vector, &g_sa.vector, &m.matrix);
#endif

#if MFIELD_FIT_EXTFIELD
          /* update (J^T J)_32 */
          {
            gsl_vector_view v32 = gsl_matrix_subrow(JTJ, extidx, euler_idx, 3);
            gsl_blas_daxpy(dB_ext * weight, &v.vector, &v32.vector);
          }
#endif
        }
    }
#endif /* MFIELD_FIT_EULER */

  return GSL_SUCCESS;
}

/*
mfield_jacobian_row_F()
  Construct a row of the Jacobian matrix corresponding to
a scalar measurement and update J^T J matrix and op(J) u
vector

Inputs: TransJ      - op(J)
        t           - scaled timestamp
        flags       - MAGDATA_FLG_xxx flags for this data point
        weight      - weight for this data point
        u           - input vector
        ridx        - index of this row in [0,nres-1]
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
        J_int       - (output) row of Jacobian (weighted) for internal
                               Green's functions, p_int-by-1
        JTJ         - (output) updated J^T W J matrix, possibly NULL
        v           - (output) updated op(J) u vector
        w           - workspace
*/

static inline int
mfield_jacobian_row_F(CBLAS_TRANSPOSE_t TransJ, const double t, const double weight,
                      const gsl_vector * u, const size_t ridx,
                      gsl_vector * dX, gsl_vector * dY, gsl_vector * dZ,
                      const double B_model[4], const size_t extidx, const double dB_ext[3],
                      gsl_vector *J_int, gsl_matrix *JTJ, gsl_vector *v,
                      const mfield_workspace *w)
{
  const double sqrt_weight = sqrt(weight);
  size_t k;
  double b[3], ui, *Ju_ptr;

  (void) extidx;
  (void) dB_ext;
  (void) JTJ;

  if (TransJ == CblasTrans)
    {
      ui = gsl_vector_get(u, ridx);
    }
  else
    {
      Ju_ptr = gsl_vector_ptr(v, ridx);
    }

  /* compute unit vector in model direction */
  for (k = 0; k < 3; ++k)
    b[k] = B_model[k] / B_model[3];

  /* compute (X dX + Y dY + Z dZ) */
  for (k = 0; k < w->nnm_mf; ++k)
    {
      double dXk = gsl_vector_get(dX, k);
      double dYk = gsl_vector_get(dY, k);
      double dZk = gsl_vector_get(dZ, k);
      double val = sqrt_weight * (b[0] * dXk +
                                  b[1] * dYk +
                                  b[2] * dZk);

      mfield_set_mf(J_int, k, val, w);
      mfield_set_sv(J_int, k, t * val, w);
      mfield_set_sa(J_int, k, 0.5 * t * t * val, w);
    }

  if (TransJ == CblasTrans)
    {
      /* update J^T u */
      gsl_vector_view vJTu = gsl_vector_subvector(v, 0, w->p_int);
      gsl_blas_daxpy(ui, J_int, &vJTu.vector);
    }
  else
    {
      /* update (J u)_i = J_int . u(1:pint) */
      gsl_vector_const_view z = gsl_vector_const_subvector(u, 0, w->p_int);
      gsl_blas_ddot(J_int, &z.vector, Ju_ptr);
    }

#if MFIELD_FIT_EXTFIELD
  {
    double val = sqrt_weight * (b[0] * dB_ext[0] +
                                b[1] * dB_ext[1] +
                                b[2] * dB_ext[2]);

    if (TransJ == CblasTrans)
      {
        /* update J^T u */
        double *ptr = gsl_vector_ptr(v, extidx);
        *ptr += val * ui;
      }
    else
      {
        double tmp = gsl_vector_get(u, extidx);
        *Ju_ptr += val * tmp;
      }

    /* update J^T J */
    if (JTJ)
      {
        gsl_vector_view v31 = gsl_matrix_subrow(JTJ, extidx, 0, w->p_int);
        double *ptr33 = gsl_matrix_ptr(JTJ, extidx, extidx);

        /* update (J^T J)_33 */
        *ptr33 += val * val;

        /* update (J^T J)_31 */
        gsl_blas_daxpy(val, J_int, &v31.vector);
      }
  }
#endif

  return GSL_SUCCESS;
}

/*
mfield_nonlinear_vector_precompute()
  Precompute J_int^T W J_int for vector measurements, since
this submatrix is independent of the model parameters and
only needs to be computed once per iteration. This function
uses OpenMP to speed up the calculation

Inputs: weights - weight vector
        w       - workspace

Notes:
1) w->JTJ_vec is updated with J_int^T W J_int for vector
residuals
*/

static int
mfield_nonlinear_vector_precompute(const gsl_vector *weights, mfield_workspace *w)
{
  int s = GSL_SUCCESS;
  size_t i, j;

  gsl_matrix_set_zero(w->JTJ_vec);

  /* check for quick return */
  if (w->nres_vec == 0)
    return GSL_SUCCESS;

  /*
   * omp_rowidx[thread_id] contains the number of currently filled rows
   * of omp_J[thread_id]. When omp_J[thread_id] is filled, it is folded
   * into the JTJ_vec matrix and then omp_rowidx[thread_id] is reset to 0
   */
  for (i = 0; i < w->max_threads; ++i)
    w->omp_rowidx[i] = 0;

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          size_t ridx = mptr->index[j]; /* residual index for this data point in [0:nres-1] */
          double t = mptr->ts[j];
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          gsl_vector_view vx, vy, vz;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          vx = gsl_matrix_row(w->omp_dX, thread_id);
          vy = gsl_matrix_row(w->omp_dY, thread_id);
          vz = gsl_matrix_row(w->omp_dZ, thread_id);

          green_calc_int(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data,
                         w->green_array_p[thread_id]);

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              double wj = gsl_vector_get(weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green(t, wj, &vx.vector, &v.vector, w);
                }
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              double wj = gsl_vector_get(weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green(t, wj, &vy.vector, &v.vector, w);
                }
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              double wj = gsl_vector_get(weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green(t, wj, &vz.vector, &v.vector, w);
                }
            }

          if (w->omp_rowidx[thread_id] >= 4*w->data_block - 3)
            {
              /* fold current matrix block into JTJ_vec, one thread at a time */
              gsl_matrix_view m = gsl_matrix_submatrix(w->omp_J[thread_id], 0, 0, w->omp_rowidx[thread_id], w->p_int);

#pragma omp critical
              {
                gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &m.matrix, 1.0, w->JTJ_vec);
              }

              w->omp_rowidx[thread_id] = 0;
            }
        } /* for (j = 0; j < mptr->n; ++j) */
    } /* for (i = 0; i < w->nsat; ++i) */

  /* now loop through to see if any rows were not accumulated into JTJ_vec */
  for (i = 0; i < w->max_threads; ++i)
    {
      if (w->omp_rowidx[i] > 0)
        {
          /* accumulate final Green's functions into JTJ_vec */
          gsl_matrix_view m = gsl_matrix_submatrix(w->omp_J[i], 0, 0, w->omp_rowidx[i], w->p_int);
          gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &m.matrix, 1.0, w->JTJ_vec);
        }
    }

  return s;
}

/*
mfield_vector_green()
  Function to compute sqrt(w) [ J_mf J_sv J_sa ] for a given set of
vector Green's functions

Inputs: t      - scaled timestamp
        weight - weight for this measurement
        g      - vector Green's functions J_mf, size nnm_mf
        G      - (output) combined vector G = sqrt(w) [ g ; t*g ; 0.5*t*t*g ],
                 size w->p_int
        w      - workspace
*/

static int
mfield_vector_green(const double t, const double weight, const gsl_vector *g,
                    gsl_vector *G, mfield_workspace *w)
{
  const double sqrt_weight = sqrt(weight);
  size_t i;

  /* form G */
  for (i = 0; i < w->nnm_mf; ++i)
    {
      double gi = sqrt_weight * gsl_vector_get(g, i);

      mfield_set_mf(G, i, gi, w);
      mfield_set_sv(G, i, t * gi, w);
      mfield_set_sa(G, i, 0.5 * t * t * gi, w);
    }

  return GSL_SUCCESS;
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
                           double dB[3], const mfield_workspace *w)
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

  (void) r;
  (void) theta;
  (void) phi;
  (void) g;
  (void) w;
  
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
  const size_t nmin_sv = 13;
  const size_t nmin_sa = 8;
  const size_t nmax = w->nmax_mf;
  const double c = 3485.0;       /* Earth core radius */
  const double a = MFIELD_RE_KM; /* Earth surface radius */
  const double ratio = a / c;
  size_t n;
  int m;
  double lambda_mf = w->lambda_mf;
  double lambda_sv = w->lambda_sv;
  double lambda_sa = w->lambda_sa;

  gsl_vector_set_zero(diag);

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;
      double term = (n + 1.0) / sqrt(2.0*n + 1.0) * pow(ratio, n + 2.0);

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);

          mfield_set_mf(diag, cidx, lambda_mf, w);

#if 0
          if (n >= nmin_sv)
            mfield_set_sv(diag, cidx, lambda_sv, w);
#else
          if (n >= nmin_sv)
            mfield_set_sv(diag, cidx, lambda_sv * n * n * n, w);
#endif

          if (n >= nmin_sa)
            mfield_set_sa(diag, cidx, lambda_sa * term, w);
          else
            mfield_set_sa(diag, cidx, lambda_mf, w);
        }
    }

  return s;
} /* mfield_nonlinear_regularize() */

static void
mfield_nonlinear_callback(const size_t iter, void *params,
                          const gsl_multifit_nlinear_workspace *multifit_p)
{
  mfield_workspace *w = (mfield_workspace *) params;
  gsl_vector *x = gsl_multifit_nlinear_position(multifit_p);
  gsl_vector *f = gsl_multifit_nlinear_residual(multifit_p);
  double avratio = gsl_multifit_nlinear_avratio(multifit_p);
  double rcond;

  /* print out state every 5 iterations */
  if (iter % 5 != 0 && iter != 1)
    return;

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

  fprintf(stderr, "\t |a|/|v|:    %12g\n", avratio);
  fprintf(stderr, "\t ||f(x)||:   %12g\n", gsl_blas_dnrm2(f));

  gsl_multifit_nlinear_rcond(&rcond, multifit_p);
  fprintf(stderr, "\t cond(J(x)): %12g\n", 1.0 / rcond);
}

static void
mfield_nonlinear_callback2(const size_t iter, void *params,
                           const gsl_multilarge_nlinear_workspace *multilarge_p)
{
  mfield_workspace *w = (mfield_workspace *) params;
  gsl_vector *x = gsl_multilarge_nlinear_position(multilarge_p);
  gsl_vector *dx = gsl_multilarge_nlinear_step(multilarge_p);
  gsl_vector *f = gsl_multilarge_nlinear_residual(multilarge_p);
  double avratio = gsl_multilarge_nlinear_avratio(multilarge_p);
  double rcond, max_dx = 0.0;
  size_t i;

  /* print out state every 5 iterations */
  if (iter % 5 != 0 && iter != 1)
    return;

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

  /* compute max | dx_i / x_i | */
  for (i = 0; i < w->p; ++i)
    {
      double dxi = gsl_vector_get(dx, i);
      double xi = gsl_vector_get(x, i);

      if (fabs(xi) < GSL_DBL_EPSILON)
        continue;

      max_dx = GSL_MAX(max_dx, fabs(dxi / xi));
    }

  fprintf(stderr, "\t max |dx/x|:  %12g\n", max_dx);
  fprintf(stderr, "\t |a|/|v|:     %12g\n", avratio);
  fprintf(stderr, "\t |f(x)|:      %12g\n", gsl_blas_dnrm2(f));

  gsl_multilarge_nlinear_rcond(&rcond, multilarge_p);
  fprintf(stderr, "\t cond(J(x)):  %12g\n", 1.0 / rcond);
}
