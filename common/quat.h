/*
 * quat.h
 */

#ifndef INCLUDED_quat_h
#define INCLUDED_quat_h

/*
 * Prototypes
 */

int quat_apply(const double q[4], const double input[3], double output[3]);
int quat_apply_inverse(const double q[4], const double input[3], double output[3]);
int quat_q2R(const double q[4], gsl_matrix *R);
int quat_R2q(const gsl_matrix * R, double q[4]);

#endif /* INCLUDED_quat_h */
