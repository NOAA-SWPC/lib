/*
 * test.c
 * Patrick Alken
 *
 * Test MSIS module
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "msis.h"

int
main(void)
{
  size_t i;
  time_t t;
  double theta, phi;
  double altmin, altmax, altstp;
  size_t nalt;
  msis_workspace *msis_workspace_p;

  t = 1111276800;
  theta = M_PI / 2.0 - 10.0 * M_PI / 180.0;
  phi = 15.0 * M_PI / 180.0;

  altmin = 50.0;
  altstp = 2.0;
  nalt = 300;

  altmax = altmin + nalt * altstp;

  msis_workspace_p = msis_alloc(nalt, F107_IDX_FILE);

  /*msis_f107_override(175.0, 175.0, msis_workspace_p);*/

  msis_calc(theta, phi, t, altmin, altstp, nalt, msis_workspace_p);

  /* print results */
  i = 1;
  printf("# Field %zu: altitude (km)\n", i++);
  printf("# Field %zu: He density (m^-3)\n", i++);
  printf("# Field %zu: O density (m^-3)\n", i++);
  printf("# Field %zu: N2 density (m^-3)\n", i++);
  printf("# Field %zu: O2 density (m^-3)\n", i++);
  printf("# Field %zu: Ar density (m^-3)\n", i++);
  printf("# Field %zu: H density (m^-3)\n", i++);
  printf("# Field %zu: N density (m^-3)\n", i++);
  printf("# Field %zu: total density (m^-3)\n", i++);
  printf("# Field %zu: neutral temperature (K)\n", i++);

  for (i = 0; i < nalt; ++i)
    {
      double alt = altmin + i * altstp;
      msis_result *result = msis_get_result(i, msis_workspace_p);

      printf("%f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
             alt,
             result->n_He,
             result->n_O,
             result->n_N2,
             result->n_O2,
             result->n_Ar,
             result->n_H,
             result->n_N,
             result->n_n,
             result->T_n);
    }

  msis_free(msis_workspace_p);

  return 0;
} /* main() */
