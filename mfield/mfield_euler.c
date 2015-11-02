
static int mfield_euler_R3(const double alpha, const double beta,
                           const double gamma, gsl_matrix *R3);
static int mfield_euler_Rq(const double *q, gsl_matrix *Rq);
static int mfield_euler_apply_Rx(int deriv, const double alpha,
                                 const double A_in[3], double A_out[3]);
static int mfield_euler_apply_Ry(int deriv, const double alpha,
                                 const double A_in[3], double A_out[3]);
static int mfield_euler_apply_Rz(int deriv, const double alpha,
                                 const double A_in[3], double A_out[3]);
static size_t mfield_euler_bin(const size_t sat_idx, const double t,
                               const mfield_workspace *w);

/*
mfield_euler_idx()
  Return index of Euler angles for a given satellite and time period

Inputs: sat_idx - satellite index in [0,nsat-1]
        t       - time (CDF_EPOCH)
        w       - workspace

Return: index so that:
alpha = w->c[index]
beta = w->c[index + 1]
gamma = w->c[index + 2]

Notes:
1) The Euler angles are organized inside the coefficient vector as:

c = [ MF | SV | SA | Euler | Ext ]

and

Euler = [ Euler_0 | Euler_1 | ... | Euler_n ]

where Euler_k are the Euler angles for satellite k. Now

Euler_k = [ Euler_{k0} | Euler_{k1} | ... | Euler_{km} ]

and m = nbins - 1. The appropriate bin is selected according to
the time, as Euler angles are computed for regular time intervals,
specified by the euler_period parameter. Finally, for satellite k and
bin l, the Euler angles are:

Euler_{kl} = [ alpha beta gamma ]
*/

size_t
mfield_euler_idx(const size_t sat_idx, const double t, const mfield_workspace *w)
{
  /*size_t idx = w->euler_offset + 3 * sat_idx;*/
  size_t bin = mfield_euler_bin(sat_idx, t, w);
  size_t idx = w->euler_offset + w->offset_euler[sat_idx] + 3 * bin;

  assert(idx >= w->euler_offset && idx < w->euler_offset + w->neuler);

  return idx;
} /* mfield_euler_idx() */

/*
mfield_nec2vfm()
  Convert a vector from NEC to VFM frame

Inputs: alpha - Euler angle alpha (radians)
        beta  - Euler angle beta (radians)
        gamma - Euler angle gamma (radians)
        q     - quaternions for CRF to NEC
        B_in  - input vector (NEC frame)
        B_out - output vector (VFM frame)

Return: success or error

Notes:
1) In place transform is allowed (ie: B_in == B_out)
*/

int
mfield_nec2vfm(const double alpha, const double beta,
               const double gamma, const double q[],
               const double B_in[3], double B_out[3])
{
  double R3_data[9], Rq_data[9], tmp_data[3];
  gsl_matrix_view R3 = gsl_matrix_view_array(R3_data, 3, 3);
  gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
  gsl_vector_const_view in = gsl_vector_const_view_array(B_in, 3);
  gsl_vector_view out = gsl_vector_view_array(B_out, 3);
  gsl_vector_view tmp = gsl_vector_view_array(tmp_data, 3);

  /* construct Euler rotation matrix R3 (Olsen et l, 2013, eq 2) */
  mfield_euler_R3(alpha, beta, gamma, &R3.matrix);

  /* construct quaternion rotation matrix R_q (Olsen et al, 2013, eq 3) */
  mfield_euler_Rq(q, &Rq.matrix);

  /* apply (inverse) quaternion rotation for NEC to CRF: B_tmp = R_q^T B_in */
  gsl_blas_dgemv(CblasTrans, 1.0, &Rq.matrix, &in.vector, 0.0, &tmp.vector);

  /* apply (inverse) Euler rotations for CRF to VFM: B_out = R_3^T B_tmp */
  gsl_blas_dgemv(CblasTrans, 1.0, &R3.matrix, &tmp.vector, 0.0, &out.vector);

  return 0;
} /* mfield_nec2vfm() */

/*
mfield_vfm2nec()
  Convert a vector from VFM to NEC frame

Inputs: alpha - Euler angle alpha (radians)
        beta  - Euler angle beta (radians)
        gamma - Euler angle gamma (radians)
        q     - quaternions for CRF to NEC
        B_in  - input vector (VFM frame)
        B_out - output vector (NEC frame)

Return: success or error

Notes:
1) In place transform is allowed (ie: B_in == B_out)
*/

int
mfield_vfm2nec(const double alpha, const double beta,
               const double gamma, const double q[],
               const double B_in[3], double B_out[3])
{
  double R3_data[9], Rq_data[9], tmp_data[3];
  gsl_matrix_view R3 = gsl_matrix_view_array(R3_data, 3, 3);
  gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
  gsl_vector_const_view in = gsl_vector_const_view_array(B_in, 3);
  gsl_vector_view out = gsl_vector_view_array(B_out, 3);
  gsl_vector_view tmp = gsl_vector_view_array(tmp_data, 3);

  /* construct Euler rotation matrix R3 (Olsen et l, 2013, eq 2) */
  mfield_euler_R3(alpha, beta, gamma, &R3.matrix);

  /* construct quaternion rotation matrix R_q (Olsen et al, 2013, eq 3) */
  mfield_euler_Rq(q, &Rq.matrix);

  /* apply Euler rotations for VFM to CRF: B_tmp = R_3 B_in */
  gsl_blas_dgemv(CblasNoTrans, 1.0, &R3.matrix, &in.vector, 0.0, &tmp.vector);

  /* apply quaternion rotation for CRF to NEC: B_out = R_q B_tmp */
  gsl_blas_dgemv(CblasNoTrans, 1.0, &Rq.matrix, &tmp.vector, 0.0, &out.vector);

  return 0;
} /* mfield_vfm2nec() */

/*
mfield_vfm2nec_deriv()
  Convert a vector from VFM to NEC frame

Inputs: flags - derivative flags
        alpha - Euler angle alpha (radians)
        beta  - Euler angle beta (radians)
        gamma - Euler angle gamma (radians)
        q     - quaternions for CRF to NEC
        B_in  - input vector (VFM frame)
        B_out - output vector (NEC frame)

Return: success or error

Notes:
1) In place transform is allowed (ie: B_in == B_out)
*/

int
mfield_vfm2nec_deriv(const size_t flags, const double alpha, const double beta,
                     const double gamma, const double q[],
                     const double B_in[3], double B_out[3])
{
  double Rq_data[9], tmp_data[3];
  gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
  gsl_vector_view out = gsl_vector_view_array(B_out, 3);
  gsl_vector_view tmp = gsl_vector_view_array(tmp_data, 3);
  double T1[3], T2[3];
  int deriv_alpha = flags & MFIELD_EULER_DERIV_ALPHA;
  int deriv_beta = flags & MFIELD_EULER_DERIV_BETA;
  int deriv_gamma = flags & MFIELD_EULER_DERIV_GAMMA;

  /* construct quaternion rotation matrix R_q (Olsen et al, 2013, eq 3) */
  mfield_euler_Rq(q, &Rq.matrix);

  /* apply Euler rotations for VFM to CRF: B_tmp = R_3 B_in */
  mfield_euler_apply_Rx(deriv_alpha, alpha, B_in, T1);
  mfield_euler_apply_Ry(deriv_beta, beta, T1, T2);
  mfield_euler_apply_Rz(deriv_gamma, gamma, T2, tmp_data);

  /* apply quaternion rotation for CRF to NEC: B_out = R_q B_tmp */
  gsl_blas_dgemv(CblasNoTrans, 1.0, &Rq.matrix, &tmp.vector, 0.0, &out.vector);

  return 0;
} /* mfield_vfm2nec_deriv() */

int
mfield_euler_print(const char *filename, const size_t sat_idx,
                   const mfield_workspace *w)
{
  FILE *fp;
  size_t i;
  double t0 = w->data_workspace_p->t0[sat_idx];
  double t1 = w->data_workspace_p->t1[sat_idx];
  double t, dt;
  double alpha_mean = 0.0,
         beta_mean = 0.0,
         gamma_mean = 0.0;
  size_t nmean = 0;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "mfield_euler_print: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Euler angles for satellite %zu\n", sat_idx);
  fprintf(fp, "# Field %zu: timestamp (CDF_EPOCH)\n", i++);
  fprintf(fp, "# Field %zu: time (decimal years)\n", i++);
  fprintf(fp, "# Field %zu: alpha (degrees)\n", i++);
  fprintf(fp, "# Field %zu: beta (degrees)\n", i++);
  fprintf(fp, "# Field %zu: gamma (degrees)\n", i++);
  fprintf(fp, "# Field %zu: alpha - mean (arcseconds)\n", i++);
  fprintf(fp, "# Field %zu: beta - mean (arcseconds)\n", i++);
  fprintf(fp, "# Field %zu: gamma - mean (arcseconds)\n", i++);

  if (w->params.euler_period > 0.0)
    dt = w->params.euler_period * 86400000; /* convert to ms */
  else
    dt = t1 - t0;

  /* compute means of Euler angles */
  for (t = t0; t < t1; t += dt)
    {
      size_t euler_idx = mfield_euler_idx(sat_idx, t, w);
      double alpha = gsl_vector_get(w->c, euler_idx);
      double beta = gsl_vector_get(w->c, euler_idx + 1);
      double gamma = gsl_vector_get(w->c, euler_idx + 2);

      if (alpha == 0.0 || beta == 0.0 || gamma == 0.0)
        continue;

      alpha_mean += alpha;
      beta_mean += beta;
      gamma_mean += gamma;
      ++nmean;
    }

  if (nmean > 0)
    {
      alpha_mean /= (double) nmean;
      beta_mean /= (double) nmean;
      gamma_mean /= (double) nmean;
    }

  for (t = t0; t < t1; t += dt)
    {
      size_t euler_idx = mfield_euler_idx(sat_idx, t, w);
      double alpha = gsl_vector_get(w->c, euler_idx);
      double beta = gsl_vector_get(w->c, euler_idx + 1);
      double gamma = gsl_vector_get(w->c, euler_idx + 2);

      fprintf(fp, "%f %f %.10f %.10f %.10f %f %f %f\n",
              t,
              satdata_epoch2year(t),
              wrap180(alpha * 180.0 / M_PI),
              wrap180(beta * 180.0 / M_PI),
              wrap180(gamma * 180.0 / M_PI),
              (alpha - alpha_mean) * 180.0 / M_PI * 3600.0,
              (beta - beta_mean) * 180.0 / M_PI * 3600.0,
              (gamma - gamma_mean) * 180.0 / M_PI * 3600.0);
    }

  fclose(fp);

  return 0;
}

/*
mfield_euler_R3()
  Construct Euler rotation matrix R_3; see
Olsen et al, 2013, eq 2

Inputs: alpha - Euler angle alpha (radians)
        beta  - Euler angle beta (radians)
        gamma - Euler angle gamma (radians)
        R3    - (output) R_3 matrix, 3-by-3

Return: success or error
*/

static int
mfield_euler_R3(const double alpha, const double beta,
                const double gamma, gsl_matrix *R3)
{
  gsl_matrix_set(R3, 0, 0, cos(beta) * cos(gamma));
  gsl_matrix_set(R3, 0, 1, sin(alpha) * sin(beta) * cos(gamma) - cos(alpha) * sin(gamma));
  gsl_matrix_set(R3, 0, 2, cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma));

  gsl_matrix_set(R3, 1, 0, cos(beta) * sin(gamma));
  gsl_matrix_set(R3, 1, 1, cos(alpha) * cos(gamma) + sin(alpha) * sin(beta) * sin(gamma));
  gsl_matrix_set(R3, 1, 2, -sin(alpha) * cos(gamma) + cos(alpha) * sin(beta) * sin(gamma));

  gsl_matrix_set(R3, 2, 0, -sin(beta));
  gsl_matrix_set(R3, 2, 1, sin(alpha) * cos(beta));
  gsl_matrix_set(R3, 2, 2, cos(alpha) * cos(beta));

  return GSL_SUCCESS;
} /* mfield_euler_R3() */

/*
mfield_euler_Rq()
  Construct quaternion rotation matrix R_q; see
Olsen et al, 2013, eq 3

Note: the (3,2) element of the equation 3 is wrong, it should be sign
inverted: 2*(q2*q3 - q1*q4)

Inputs: q  - quaternion vector
        Rq - (output) R_q matrix, 3-by-3

Return: success or error
*/

static int
mfield_euler_Rq(const double *q, gsl_matrix *Rq)
{
  const double q1 = q[0];
  const double q2 = q[1];
  const double q3 = q[2];
  const double q4 = q[3];

  gsl_matrix_set(Rq, 0, 0, 1.0 - 2.0*q2*q2 - 2.0*q3*q3);
  gsl_matrix_set(Rq, 0, 1, 2.0*(q1*q2 + q3*q4));
  gsl_matrix_set(Rq, 0, 2, 2.0*(q1*q3 - q2*q4));

  gsl_matrix_set(Rq, 1, 0, 2.0*(q1*q2 - q3*q4));
  gsl_matrix_set(Rq, 1, 1, 1.0 - 2.0*q1*q1 - 2.0*q3*q3);
  gsl_matrix_set(Rq, 1, 2, 2.0*(q2*q3 + q1*q4));

  gsl_matrix_set(Rq, 2, 0, 2.0*(q1*q3 + q2*q4));
  gsl_matrix_set(Rq, 2, 1, 2.0*(q2*q3 - q1*q4));
  gsl_matrix_set(Rq, 2, 2, 1.0 - 2.0*q1*q1 - 2.0*q2*q2);

  return GSL_SUCCESS;
}

/*
mfield_euler_apply_Rx()
  Apply a rotation matrix around the x axis

Inputs: alpha - rotation angle in radians
        A_in  - input vector
        A_out - (output) rotated vector
*/

static int
mfield_euler_apply_Rx(int deriv, const double alpha, const double A_in[3], double A_out[3])
{
  int s = 0;
  double x = alpha;

  if (deriv)
    {
      /* d/dalpha R(alpha) = R(alpha + pi/2) */
      x += M_PI / 2.0;
      A_out[0] = 0.0;
    }
  else
    A_out[0] = A_in[0];

  A_out[1] = A_in[1] * cos(x) - A_in[2] * sin(x);
  A_out[2] = A_in[1] * sin(x) + A_in[2] * cos(x);

  return s;
} /* mfield_euler_apply_Rx() */

/*
mfield_euler_apply_Ry()
  Apply a rotation matrix around the y axis

Inputs: alpha - rotation angle in radians
        A_in  - input vector
        A_out - (output) rotated vector
*/

static int
mfield_euler_apply_Ry(int deriv, const double alpha, const double A_in[3], double A_out[3])
{
  int s = 0;
  double x = alpha;

  if (deriv)
    {
      /* d/dalpha R(alpha) = R(alpha + pi/2) */
      x += M_PI / 2.0;
      A_out[1] = 0.0;
    }
  else
    A_out[1] = A_in[1];

  A_out[0] = A_in[0] * cos(x) + A_in[2] * sin(x);
  A_out[2] = -A_in[0] * sin(x) + A_in[2] * cos(x);

  return s;
} /* mfield_euler_apply_Ry() */

/*
mfield_euler_apply_Rz()
  Apply a rotation matrix around the z axis

Inputs: deriv - 1 if computing derivative wrt alpha, 0 if not
        alpha - rotation angle in radians
        A_in  - input vector
        A_out - (output) rotated vector
*/

static int
mfield_euler_apply_Rz(int deriv, const double alpha, const double A_in[3], double A_out[3])
{
  int s = 0;
  double x = alpha;

  if (deriv)
    {
      /* d/dalpha R(alpha) = R(alpha + pi/2) */
      x += M_PI / 2.0;
      A_out[2] = 0.0;
    }
  else
    A_out[2] = A_in[2];

  A_out[0] = A_in[0] * cos(x) - A_in[1] * sin(x);
  A_out[1] = A_in[0] * sin(x) + A_in[1] * cos(x);

  return s;
} /* mfield_euler_apply_Rz() */

/*
mfield_euler_bin()
  Return bin number of Euler angles for a given satellite and time

Inputs: sat_idx - satellite index in [0,nsat-1]
        t       - time (CDF_EPOCH)
        w       - workspace

Return: bin number corresponding to Euler angle triplet
*/

static size_t
mfield_euler_bin(const size_t sat_idx, const double t, const mfield_workspace *w)
{
  double t0 = w->data_workspace_p->t0[sat_idx];
  double t1 = w->data_workspace_p->t1[sat_idx];
  double ratio = (t - t0) / (t1 - t0);
  size_t bin = (size_t) (w->nbins_euler[sat_idx] * ratio);

  if (bin >= w->nbins_euler[sat_idx])
    bin = w->nbins_euler[sat_idx] - 1;

  return bin;
}
