/*
 * test.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "euler.h"

int
test_quaternion(gsl_rng *r, double q[4])
{
  double u[3];  /* unit vector of quaternion */
  double unorm; /* |u| */
  double theta; /* angle of rotation for quaternion */
  size_t j;

  /* random angle */
  theta = gsl_rng_uniform(r) * 2.0 * M_PI;

  /* random unit vector */
  for (j = 0; j < 3; ++j)
    u[j] = gsl_rng_uniform(r);

  unorm = gsl_hypot3(u[0], u[1], u[2]);

  for (j = 0; j < 3; ++j)
    {
      u[j] /= unorm;
      q[j] = u[j] * sin(theta / 2.0);
    }

  q[3] = cos(theta / 2.0);

  return 0;
}

int
test_euler(euler_workspace *w)
{
  int s = 0;
  const double tol = 1.0e-12;
  const double t0 = w->t[0];
  const double t1 = w->t[w->n - 1];
  const double dt = 5.0 * 86400.0 * 1000.0;
  double B0[3], B1[3], B2[3], q[4];
  double t;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  size_t j;

  for (t = t0; t <= t1; t += dt)
    {
      /* make random quaternion */
      test_quaternion(r, q);

      /* random input vector */
      for (j = 0; j < 3; ++j)
        B0[j] = gsl_rng_uniform(r) * 10.0;

      s += euler_nec2vfm_t(t, q, B0, B1, w);
      s += euler_vfm2nec_t(t, q, B1, B2, w);

      gsl_test_rel(B2[0], B0[0], tol, "t=%f X", t);
      gsl_test_rel(B2[1], B0[1], tol, "t=%f Y", t);
      gsl_test_rel(B2[2], B0[2], tol, "t=%f Z", t);
    }

  gsl_rng_free(r);

  return s;
} /* test_euler() */

int
main(int argc, char *argv[])
{
  euler_workspace *euler_p;
  char *euler_file = "euler.0";

  fprintf(stderr, "main: reading Euler angles from %s...", euler_file);
  euler_p = euler_read(euler_file);
  fprintf(stderr, "done (%zu sets of angles read)\n", euler_p->n);

  test_euler(euler_p);

  euler_free(euler_p);

  exit (gsl_test_summary());
} /* main() */
