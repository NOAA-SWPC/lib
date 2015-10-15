/*
 * chaos.c
 * Wrapper for CHAOS external field model (links to MATLAB)
 *
 * Calling program must initialize MATLAB as follows:
 *
 * if (!mclInitializeApplication(NULL,0))
 *   {
 *     fprintf(stderr, "main: failed to initialize\n");
 *     exit(1);
 *   }
 *
 * if (!libchaosextInitialize())
 *   {
 *     fprintf(stderr, "main: failed to initialize\n");
 *     exit(1);
 *   }
 *
 * and terminate with:
 *
 * libchaosextTerminate();
 * mclTerminateApplication();
 */

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <satdata/satdata.h>

#include "chaos.h"
#include "libchaosext.h"

int
chaos_init(void)
{
  fprintf(stderr, "chaos_init: initializing MATLAB runtime...");

  if (!mclInitializeApplication(NULL,0))
    {
      fprintf(stderr, "chaos_init: failed to initialize\n");
      exit(1);
    }

  if (!libchaosextInitialize())
    {
      fprintf(stderr, "chaos_init: failed to initialize\n");
      exit(1);
    }

  fprintf(stderr, "done\n");

  return 0;
}

int
chaos_term(void)
{
  libchaosextTerminate();
  mclTerminateApplication();
}

/*
chaos_ext()
  Wrapper for CHAOS external field model

Inputs: tv    - CDF_EPOCH
        altv  - geocentric altitude (km)
        latv  - geocentric latitude (deg)
        phiv  - geocentric longitude (deg)
        RC_e  - external RC index (nT)
        RC_i  - internal RC index (nT)
        B_ext - (output) external field vector in NEC (nT)
                B_ext[0] = X (nT)
                B_ext[1] = Y (nT)
                B_ext[2] = Z (nT)
*/

int
chaos_ext(const gsl_vector *tv, const gsl_vector *altv,
          const gsl_vector *latv, const gsl_vector *phiv,
          const gsl_vector *RC_ev, const gsl_vector *RC_iv,
          gsl_matrix *B_ext)
{
  const size_t n = tv->size;
  int s = 0;

  mxArray *t_a = mxCreateDoubleMatrix(n, 1, mxREAL);
  mxArray *r_a = mxCreateDoubleMatrix(n, 1, mxREAL);
  mxArray *theta_a = mxCreateDoubleMatrix(n, 1, mxREAL);
  mxArray *phi_a = mxCreateDoubleMatrix(n, 1, mxREAL);
  mxArray *RC_e_a = mxCreateDoubleMatrix(n, 1, mxREAL);
  mxArray *RC_i_a = mxCreateDoubleMatrix(n, 1, mxREAL);
  mxArray *B_a = mxCreateDoubleMatrix(n, 3, mxREAL);

  double *t = mxGetPr(t_a);
  double *r = mxGetPr(r_a);
  double *theta = mxGetPr(theta_a);
  double *phi = mxGetPr(phi_a);
  double *RC_e = mxGetPr(RC_e_a);
  double *RC_i = mxGetPr(RC_i_a);

  size_t i;
  gsl_matrix_view Bv;

  /* fill input arrays */
  for (i = 0; i < n; ++i)
    {
      *t++ = satdata_epoch2fday(gsl_vector_get(tv, i));
      *r++ = 6371.2 + gsl_vector_get(altv, i);
      *theta++ = 90.0 - gsl_vector_get(latv, i);
      *phi++ = gsl_vector_get(phiv, i);
      *RC_e++ = gsl_vector_get(RC_ev, i);
      *RC_i++ = gsl_vector_get(RC_iv, i);
    }

  mlfSynth_chaos_ext(1,
                     &(B_a),
                     t_a,
                     r_a,
                     theta_a,
                     phi_a,
                     RC_e_a,
                     RC_i_a);

  /* store the output */
  Bv = gsl_matrix_view_array(mxGetPr(B_a), 3, n);
  gsl_matrix_transpose_memcpy(B_ext, &Bv.matrix);

  for (i = 0; i < n; ++i)
    {
      double Br = gsl_matrix_get(B_ext, i, 0);
      double Bt = gsl_matrix_get(B_ext, i, 1);
      double Bp = gsl_matrix_get(B_ext, i, 2);

      gsl_matrix_set(B_ext, i, 0, -Bt);
      gsl_matrix_set(B_ext, i, 1, Bp);
      gsl_matrix_set(B_ext, i, 2, -Br);
    }

  mxDestroyArray(t_a);
  mxDestroyArray(r_a);
  mxDestroyArray(theta_a);
  mxDestroyArray(phi_a);
  mxDestroyArray(RC_e_a);
  mxDestroyArray(RC_i_a);
  mxDestroyArray(B_a);

  return s;
} /* chaos_ext() */
