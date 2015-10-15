/*
 * euler.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include <satdata/satdata.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>

#include "common.h"

#include "euler.h"

static int euler_apply_Rx(const int deriv, double x,
                          const double A_in[3], double A_out[3]);
static int euler_apply_Ry(const int deriv, double x,
                          const double A_in[3], double A_out[3]);
static int euler_apply_Rz(const int deriv, double x,
                          const double A_in[3], double A_out[3]);
static int euler_find(const double t, size_t *idx, const euler_workspace *w);

/*
euler_alloc()
  Allocate Euler workspace

Inputs: flags - EULER_FLG_xxx for Euler convention

Return: pointer to workspace
*/

euler_workspace *
euler_alloc(const size_t flags)
{
  euler_workspace *w;

  w = calloc(1, sizeof(euler_workspace));
  if (!w)
    return 0;

  w->n = 0;
  w->flags = flags;

  return w;
} /* euler_alloc() */

void
euler_free(euler_workspace *w)
{
  free(w);
} /* euler_free() */

euler_workspace *
euler_read(const char *filename)
{
  FILE *fp;
  char buf[EULER_MAX_BUFFER];
  euler_workspace *w = NULL;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "euler_read: fopen: cannot open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  while (fgets(buf, EULER_MAX_BUFFER, fp) != NULL)
    {
      int c;
      double t, year, alpha, beta, gamma;

      /* search for flags to allocate workspace */
      if (!w)
        {
          size_t flags;

          c = sscanf(buf, "# Flags: %zu", &flags);
          if (c < 1)
            continue;

          w = euler_alloc(flags);
          if (!w)
            return 0;
        }

      if (*buf == '#')
        continue;

      c = sscanf(buf, "%lf %lf %lf %lf %lf",
                 &t,
                 &year,
                 &alpha,
                 &beta,
                 &gamma);
      if (c < 5)
        continue;

      if (w->n > 0)
        {
          /* ensure time series is sorted ascending */
          if (t <= w->t[w->n - 1])
            {
              fprintf(stderr, "euler_read: error: t < t_{-1}: %f, %f\n",
                      t, w->t[w->n - 1]);
              continue;
            }
        }

      w->t[w->n] = t;
      w->alpha[w->n] = alpha * M_PI / 180.0;
      w->beta[w->n] = beta * M_PI / 180.0;
      w->gamma[w->n] = gamma * M_PI / 180.0;

      if (++(w->n) >= EULER_MAX_BINS)
        {
          fprintf(stderr, "euler_read: MAX_BINS too small: %zu\n", w->n);
          break;
        }
    }

  fclose(fp);

  return w;
} /* euler_read() */

int
euler_write(const char *filename, const euler_workspace *w)
{
  FILE *fp;
  size_t i;
  double alpha_mean, beta_mean, gamma_mean;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "euler_write: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  fprintf(fp, "# Flags: %zu\n", w->flags);

  i = 1;
  fprintf(fp, "# Field %zu: timestamp (CDF_EPOCH)\n", i++);
  fprintf(fp, "# Field %zu: time (decimal years)\n", i++);
  fprintf(fp, "# Field %zu: alpha (degrees)\n", i++);
  fprintf(fp, "# Field %zu: beta (degrees)\n", i++);
  fprintf(fp, "# Field %zu: gamma (degrees)\n", i++);
  fprintf(fp, "# Field %zu: alpha - mean (arcseconds)\n", i++);
  fprintf(fp, "# Field %zu: beta - mean (arcseconds)\n", i++);
  fprintf(fp, "# Field %zu: gamma - mean (arcseconds)\n", i++);

  /* compute means of Euler angles */
  alpha_mean = gsl_stats_mean(w->alpha, 1, w->n);
  beta_mean = gsl_stats_mean(w->beta, 1, w->n);
  gamma_mean = gsl_stats_mean(w->gamma, 1, w->n);

  for (i = 0; i < w->n; ++i)
    {
      double t = w->t[i];
      double alpha = w->alpha[i];
      double beta = w->beta[i];
      double gamma = w->gamma[i];

      fprintf(fp, "%f %.4f %.10f %.10f %.10f %f %f %f\n",
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
euler_add()
  Add a set of Euler angles to workspace at a given time. Data must
be added in ascending time order

Inputs: t - timestamp (CDF_EPOCH)
        x - 3-vector of Euler angles in radians
        w - workspace

Return: success/error
*/

int
euler_add(const double t, const gsl_vector *x, euler_workspace *w)
{
  int s = 0;
  const double alpha = gsl_vector_get(x, EULER_IDX_ALPHA);
  const double beta = gsl_vector_get(x, EULER_IDX_BETA);
  const double gamma = gsl_vector_get(x, EULER_IDX_GAMMA);

  w->t[w->n] = t;
  w->alpha[w->n] = alpha;
  w->beta[w->n] = beta;
  w->gamma[w->n] = gamma;

  ++(w->n);

  return s;
} /* euler_add() */

/*
euler_apply()
  Apply Euler angle rotations to satellite dataset

Inputs: data  - satellite data
        w     - workspace

Return: success or error

Notes:
1) data->B is replaced with data->B_VFM vectors rotated
with Euler angles and satellite quaternions
*/

int
euler_apply(satdata_mag *data, const euler_workspace *w)
{
  int s = 0;
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      double *q = &(data->q[4 * i]);
      double *B_VFM = &(data->B_VFM[3 * i]);
      double *B = &(data->B[3 * i]);
      double t = data->t[i];

      s += euler_vfm2nec_t(t, q, B_VFM, B, w);
    }

  return s;
} /* euler_apply() */

/*
euler_nec2vfm_t()
  Convert a vector from NEC to VFM frame at a given time

Inputs: t     - timestamp (CDF_EPOCH)
        q     - quaternions for CRF to NEC
        B_in  - input vector (NEC frame)
        B_out - output vector (VFM frame)

Return: success or error

Notes:
1) In place transform is allowed (ie: B_in == B_out)
*/

int
euler_nec2vfm_t(const double t, const double q[],
                const double B_in[3], double B_out[3],
                const euler_workspace *w)
{
  int s = 0;
  size_t idx;

  s = euler_find(t, &idx, w);
  if (s)
    return s;

  s = euler_nec2vfm(w->flags, w->alpha[idx], w->beta[idx], w->gamma[idx],
                    q, B_in, B_out);

  return s;
} /* euler_nec2vfm_t() */

/*
euler_vfm2nec_t()
  Convert a vector from NEC to VFM frame at a given time

Inputs: t     - timestamp (CDF_EPOCH)
        q     - quaternions for CRF to NEC
        B_in  - input vector (NEC frame)
        B_out - output vector (VFM frame)
        w     - workspace

Return: success or error

Notes:
1) In place transform is allowed (ie: B_in == B_out)
*/

int
euler_vfm2nec_t(const double t, const double q[],
                const double B_in[3], double B_out[3],
                const euler_workspace *w)
{
  int s = 0;
  size_t idx;

  s = euler_find(t, &idx, w);
  if (s)
    return s;

  s = euler_vfm2nec(w->flags, w->alpha[idx], w->beta[idx], w->gamma[idx],
                    q, B_in, B_out);

  return s;
} /* euler_vfm2nec_t() */

/*
euler_vfm2nec()
  Convert a vector from VFM to NEC frame

Inputs: flags - EULER_FLG_xxx
                These flags can indicate derivatives to be taken,
                and also Euler convention
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
euler_vfm2nec(const size_t flags, const double alpha, const double beta,
              const double gamma, const double q[],
              const double B_in[3], double B_out[3])
{
  double Rq_data[9], tmp_data[3];
  gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
  gsl_vector_view out = gsl_vector_view_array(B_out, 3);
  gsl_vector_view tmp = gsl_vector_view_array(tmp_data, 3);
  int deriv_alpha = flags & EULER_FLG_DERIV_ALPHA;
  int deriv_beta = flags & EULER_FLG_DERIV_BETA;
  int deriv_gamma = flags & EULER_FLG_DERIV_GAMMA;

  /* construct quaternion rotation matrix R_q (Olsen et al, 2013, eq 3) */
  euler_Rq(q, &Rq.matrix);

  if (flags & EULER_FLG_RINV)
    {
      /* apply R_inv = [ 0 0 -1; 0 -1 0; -1 0 0 ] matrix to input vector */
      tmp_data[0] = -B_in[2];
      tmp_data[1] = -B_in[1];
      tmp_data[2] = -B_in[0];
    }
  else
    {
      tmp_data[0] = B_in[0];
      tmp_data[1] = B_in[1];
      tmp_data[2] = B_in[2];
    }

  /* apply Euler rotations for VFM to CRF: B_tmp = R_3 B_in */
  if (flags & EULER_FLG_ZYZ)
    {
      /* ZYZ convention */
      euler_apply_Rz(deriv_alpha, alpha, tmp_data, tmp_data);
      euler_apply_Ry(deriv_beta, beta, tmp_data, tmp_data);
      euler_apply_Rz(deriv_gamma, gamma, tmp_data, tmp_data);
    }
  else if (flags & EULER_FLG_ZYX)
    {
      /* ZYX convention */
      euler_apply_Rx(deriv_alpha, alpha, tmp_data, tmp_data);
      euler_apply_Ry(deriv_beta, beta, tmp_data, tmp_data);
      euler_apply_Rz(deriv_gamma, gamma, tmp_data, tmp_data);
    }
  else
    {
      fprintf(stderr, "euler_vfm2nec: error: no Euler convention specified\n");
      return -1;
    }

  /* apply quaternion rotation for CRF to NEC: B_out = R_q B_tmp */
  gsl_blas_dgemv(CblasNoTrans, 1.0, &Rq.matrix, &tmp.vector, 0.0, &out.vector);

  return 0;
} /* euler_vfm2nec() */

/*
euler_nec2vfm()
  Convert a vector from NEC to VFM frame

Inputs: flags - EULER_FLG_xxx
                These flags can indicate derivatives to be taken,
                and also Euler convention
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
euler_nec2vfm(const size_t flags, const double alpha, const double beta,
              const double gamma, const double q[],
              const double B_in[3], double B_out[3])
{
  double Rq_data[9], tmp_data[3];
  gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
  gsl_vector_const_view in = gsl_vector_const_view_array(B_in, 3);
  gsl_vector_view tmp = gsl_vector_view_array(tmp_data, 3);
  int deriv_alpha = flags & EULER_FLG_DERIV_ALPHA;
  int deriv_beta = flags & EULER_FLG_DERIV_BETA;
  int deriv_gamma = flags & EULER_FLG_DERIV_GAMMA;

  /* construct quaternion rotation matrix R_q (Olsen et al, 2013, eq 3) */
  euler_Rq(q, &Rq.matrix);

  /* apply (inverse) quaternion rotation for NEC to CRF: B_tmp = R_q^T B_in */
  gsl_blas_dgemv(CblasTrans, 1.0, &Rq.matrix, &in.vector, 0.0, &tmp.vector);

  /*
   * apply Euler rotations for CRF to VFM: B_tmp = R_3^T B_in
   * Note that matrices are applied in reverse order due to transpose
   */
  if (flags & EULER_FLG_ZYZ)
    {
      /* ZYZ convention */
      euler_apply_Rz(deriv_gamma, -gamma, tmp_data, B_out);
      euler_apply_Ry(deriv_beta, -beta, B_out, B_out);
      euler_apply_Rz(deriv_alpha, -alpha, B_out, B_out);
    }
  else if (flags & EULER_FLG_ZYX)
    {
      /* ZYX convention */
      euler_apply_Rz(deriv_gamma, -gamma, tmp_data, B_out);
      euler_apply_Ry(deriv_beta, -beta, B_out, B_out);
      euler_apply_Rx(deriv_alpha, -alpha, B_out, B_out);
    }
  else
    {
      fprintf(stderr, "euler_nec2vfm: error: no Euler convention specified\n");
      return -1;
    }

  return 0;
} /* euler_nec2vfm() */

/*
euler_Rq()
  Construct quaternion rotation matrix R_q; see
Olsen et al, 2013, eq 3

Note: the (3,2) element of the equation 3 is wrong, it should be sign
inverted: 2*(q2*q3 - q1*q4)

Inputs: q  - quaternion vector
        Rq - (output) R_q matrix, 3-by-3

Return: success or error
*/

int
euler_Rq(const double *q, gsl_matrix *Rq)
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
euler_apply_Rx()
  Apply a rotation matrix around the x axis

Inputs: x     - rotation angle in radians
        A_in  - input vector
        A_out - (output) rotated vector

Notes:
1) In place transform is allowed (A_in == A_out)
*/

static int
euler_apply_Rx(const int deriv, double x, const double A_in[3], double A_out[3])
{
  int s = 0;
  const double A1 = A_in[1]; /* for in-place transform */
  double cn, sn;

  if (deriv)
    {
      /* d/dalpha R(alpha) = R(alpha + pi/2) */
      x += M_PI / 2.0;
      A_out[0] = 0.0;
    }
  else
    A_out[0] = A_in[0];

  cn = cos(x);
  sn = sin(x);

  A_out[1] = A1 * cn - A_in[2] * sn;
  A_out[2] = A1 * sn + A_in[2] * cn;

  return s;
} /* euler_apply_Rx() */

/*
euler_apply_Ry()
  Apply a rotation matrix around the y axis

Inputs: x     - rotation angle in radians
        A_in  - input vector
        A_out - (output) rotated vector

Notes:
1) In place transform allowed (A_in == A_out)
*/

static int
euler_apply_Ry(const int deriv, double x, const double A_in[3], double A_out[3])
{
  int s = 0;
  const double A0 = A_in[0]; /* for in-place transform */
  double cn, sn;

  if (deriv)
    {
      /* d/dalpha R(alpha) = R(alpha + pi/2) */
      x += M_PI / 2.0;
      A_out[1] = 0.0;
    }
  else
    A_out[1] = A_in[1];

  cn = cos(x);
  sn = sin(x);

  A_out[0] = A0 * cn + A_in[2] * sn;
  A_out[2] = -A0 * sn + A_in[2] * cn;

  return s;
} /* euler_apply_Ry() */

/*
euler_apply_Rz()
  Apply a rotation matrix around the z axis

Inputs: deriv - 1 if computing derivative wrt x, 0 if not
        x     - rotation angle in radians
        A_in  - input vector
        A_out - (output) rotated vector

Notes:
1) In place transform allowed (A_in == A_out)
*/

static int
euler_apply_Rz(const int deriv, double x, const double A_in[3], double A_out[3])
{
  int s = 0;
  const double A0 = A_in[0]; /* for in-place transform */
  double cn, sn;

  if (deriv)
    {
      /* d/dalpha R(alpha) = R(alpha + pi/2) */
      x += M_PI / 2.0;
      A_out[2] = 0.0;
    }
  else
    A_out[2] = A_in[2];

  cn = cos(x);
  sn = sin(x);

  A_out[0] = A0 * cn - A_in[1] * sn;
  A_out[1] = A0 * sn + A_in[1] * cn;

  return s;
} /* euler_apply_Rz() */

/*
euler_find()
  Find bin corresponding to angles for time t

Inputs: t   - timestamp (CDF_EPOCH)
        idx - (output) bin index
        w   - workspace

Return: success or error

Notes:
1) Based on gsl_interp_bsearch()

2)
 if ( t == t0 )           then  index == 0
 if ( t > t0 && t <= t1 ) then  index == 0, and sim. for other interior pts
 if ( t == tN )           then  index == N-1
 if ( t > tN )            then  index == N-1
 if ( t < t0 )            then  index == 0 
*/

static int
euler_find(const double t, size_t *idx, const euler_workspace *w)
{
  int s = 0;
  size_t ilo = 0;
  size_t ihi;
  
  if (w->n == 0)
    {
      fprintf(stderr, "euler_find: error: no data in structure\n");
      return 0;
    }

  ihi = w->n - 1;

  /* binary search of w->t */
  while (ihi > ilo + 1)
    {
      size_t i = (ihi + ilo) / 2;

      if (w->t[i] > t) 
        ihi = i;
      else      
        ilo = i;
    }
  
  *idx = ilo;

  return s;
} /* euler_find() */
