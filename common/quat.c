/*
 * quat.c
 * Routines related to quaternions
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "quat.h"

/*
quat_apply()
  Apply quaternion rotation to input vector

output = R_q * input

Inputs: q      - quaternions, length 4
        input  - input vector, length 3
        output - (output) output vector, length 3

Notes:
1) In-place transform is allowed, so input can equal output
*/

int
quat_apply(const double q[4], const double input[3], double output[3])
{
  double Rq_data[9], tmp_data[3];
  gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
  gsl_vector_view out = gsl_vector_view_array(output, 3);
  gsl_vector_view tmp = gsl_vector_view_array(tmp_data, 3);

  tmp_data[0] = input[0];
  tmp_data[1] = input[1];
  tmp_data[2] = input[2];

  /* XXX this is probably inefficient */
  quat_q2R(q, &Rq.matrix);
  gsl_blas_dgemv(CblasNoTrans, 1.0, &Rq.matrix, &tmp.vector, 0.0, &out.vector);

  return GSL_SUCCESS;
}

/*
quat_apply_inverse()
  Apply inverse quaternion rotation to input vector

output = R_q^T * input

Inputs: q      - quaternions, length 4
        input  - input vector, length 3
        output - (output) output vector, length 3

Notes:
1) In-place transform is allowed, so input can equal output
*/

int
quat_apply_inverse(const double q[4], const double input[3], double output[3])
{
  double Rq_data[9], tmp_data[3];
  gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
  gsl_vector_view out = gsl_vector_view_array(output, 3);
  gsl_vector_view tmp = gsl_vector_view_array(tmp_data, 3);

  tmp_data[0] = input[0];
  tmp_data[1] = input[1];
  tmp_data[2] = input[2];

  /* XXX this is probably inefficient */
  quat_q2R(q, &Rq.matrix);
  gsl_blas_dgemv(CblasTrans, 1.0, &Rq.matrix, &tmp.vector, 0.0, &out.vector);

  return GSL_SUCCESS;
}

/*
quat_q2R()
  Construct quaternion rotation matrix R_q; see
Olsen et al, 2013, eq 3

Note: the (3,2) element of the equation 3 is wrong, it should be sign
inverted: 2*(q2*q3 - q1*q4)

Inputs: q - quaternion vector, length 4
        R - (output) R_q matrix, 3-by-3

Return: success or error
*/

int
quat_q2R(const double q[4], gsl_matrix *R)
{
  const double q1 = q[0];
  const double q2 = q[1];
  const double q3 = q[2];
  const double q4 = q[3];

  gsl_matrix_set(R, 0, 0, 1.0 - 2.0*q2*q2 - 2.0*q3*q3);
  gsl_matrix_set(R, 0, 1, 2.0*(q1*q2 + q3*q4));
  gsl_matrix_set(R, 0, 2, 2.0*(q1*q3 - q2*q4));

  gsl_matrix_set(R, 1, 0, 2.0*(q1*q2 - q3*q4));
  gsl_matrix_set(R, 1, 1, 1.0 - 2.0*q1*q1 - 2.0*q3*q3);
  gsl_matrix_set(R, 1, 2, 2.0*(q2*q3 + q1*q4));

  gsl_matrix_set(R, 2, 0, 2.0*(q1*q3 + q2*q4));
  gsl_matrix_set(R, 2, 1, 2.0*(q2*q3 - q1*q4));
  gsl_matrix_set(R, 2, 2, 1.0 - 2.0*q1*q1 - 2.0*q2*q2);

  return GSL_SUCCESS;
}

/*
quat_R2q()
  Convert a rotation matrix to quaternion representation

Inputs: R - rotation matrix, 3-by-3
        q - (output) quaternions, length 4

Return: success/error
*/

int
quat_R2q(const gsl_matrix * R, double q[4])
{
  const double R11 = gsl_matrix_get(R, 0, 0);
  const double R12 = gsl_matrix_get(R, 0, 1);
  const double R13 = gsl_matrix_get(R, 0, 2);
  const double R21 = gsl_matrix_get(R, 1, 0);
  const double R22 = gsl_matrix_get(R, 1, 1);
  const double R23 = gsl_matrix_get(R, 1, 2);
  const double R31 = gsl_matrix_get(R, 2, 0);
  const double R32 = gsl_matrix_get(R, 2, 1);
  const double R33 = gsl_matrix_get(R, 2, 2);
  const double Trace = R11 + R22 + R33;
  const double A = 0.25 * (1.0 - Trace);

  q[0] = 0.5 * R11 + A;
  q[1] = 0.5 * R22 + A;
  q[2] = 0.5 * R33 + A;
  q[3] = 0.25 * (1.0 + Trace);

  if ((q[0] >= q[1]) && (q[0] >= q[2]) && (q[0] >= q[3]))
    {
      q[0] = sqrt(q[0]);
      q[1] = (R21 + R12) / (4.0 * q[0]);
      q[2] = (R31 + R13) / (4.0 * q[0]);
      q[3] = (R23 - R32) / (4.0 * q[0]);
    }

  if ((q[1] >= q[0]) && (q[1] >= q[2]) && (q[1] >= q[3]))
    {
      q[1] = sqrt(q[1]);
      q[2] = (R32 + R23) / (4.0 * q[1]);
      q[3] = (R31 - R13) / (4.0 * q[1]);
      q[0] = (R12 + R21) / (4.0 * q[1]);
    }

  if ((q[2] >= q[0]) && (q[2] >= q[1]) && (q[2] >= q[3]))
    {
      q[2] = sqrt(q[2]);
      q[3] = (R12 - R21) / (4.0 * q[2]);
      q[0] = (R13 + R31) / (4.0 * q[2]);
      q[1] = (R23 + R32) / (4.0 * q[2]);
    }

  if ((q[3] >= q[0]) && (q[3] >= q[1]) && (q[3] >= q[2]))
    {
      q[3] = sqrt(q[3]);
      q[0] = (R23 - R32) / (4.0 * q[3]);
      q[1] = (R31 - R13) / (4.0 * q[3]);
      q[2] = (R12 - R21) / (4.0 * q[3]);
    }

  return GSL_SUCCESS;
}
