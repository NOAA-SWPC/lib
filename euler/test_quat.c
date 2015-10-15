/*
 * test_quat.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

#include "common.h"
#include "euler.h"
#include "magdata.h"

int
test_quat(const magdata *data)
{
  int s = 0;
  size_t i, j, k;
  gsl_matrix *Rq = gsl_matrix_alloc(3, 3);
  gsl_matrix *C = gsl_matrix_alloc(3, 3);
  const double tol = 1.0e-8;

  for (i = 0; i < data->n; ++i)
    {
      double *q = &(data->q[4 * i]);
      double qnorm;
      gsl_vector_view qv = gsl_vector_view_array(q, 4);
      double u[4];
      double theta, theta_deg;

      euler_Rq(q, Rq);

      theta = 2.0 * acos(q[3]);
      theta_deg = wrap180(theta * 180.0 / M_PI);

      u[0] = q[0] / sin(0.5 * theta);
      u[1] = q[1] / sin(0.5 * theta);
      u[2] = q[2] / sin(0.5 * theta);
      u[3] = vec_norm(u);

      if (fabs(theta_deg) > 15.0)
        gsl_test_rel(u[3], 1.0, tol, "unorm i=%zu theta=%f [deg]", i, theta_deg);

      qnorm = gsl_blas_dnrm2(&qv.vector);
      gsl_test_rel(qnorm, 1.0, tol, "qnorm i=%zu theta=%f [deg]", i, theta_deg);

      gsl_blas_dsyrk(CblasUpper, CblasTrans, 1.0, Rq, 0.0, C);

      for (k = 0; k < 3; ++k)
        {
          for (j = 0; j <= k; ++j)
            {
              double Cjk = gsl_matrix_get(C, j, k);

              if (j == k)
                gsl_test_rel(Cjk, 1.0, tol, "C(%zu,%zu)", j, k);
              else
                gsl_test_rel(Cjk + 0.1, 0.1, tol, "C(%zu,%zu)", j, k);
            }
        }
    }

  return s;
}

int
main(int argc, char *argv[])
{
  magdata *data;

  if (argc < 2)
    {
      fprintf(stderr, "Usage: %s <magdata_file>\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: reading %s...", argv[1]);
  data = magdata_read(argv[1], NULL);
  fprintf(stderr, "done\n");

  test_quat(data);
    
  exit (gsl_test_summary());
} /* main() */
