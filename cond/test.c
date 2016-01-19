/*
 * test.c
 *
 * Test cond module
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "cond.h"

int
main(void)
{
  size_t i;
  time_t t;
  double theta, phi;
  double altmin, altmax, altstp;
  size_t nalt;
  cond_workspace *cond_workspace_p;

#if 1
  t = 1111276800;
  theta = M_PI / 2.0 - 10.0 * M_PI / 180.0;
  phi = 15.0 * M_PI / 180.0;

  altmin = 50.0;
  altstp = 2.0;
  nalt = 300;
#else

  t = 1048204800;
  theta = M_PI / 2.0 + 1.5271630954950381;
  phi = 170.0 * M_PI / 180.0;

  altmin = 90.0;
  altmax = 773.0;
  nalt = 29;
  altstp = (altmax - altmin) / (nalt - 1);
#endif

  altmax = altmin + nalt * altstp;

  cond_workspace_p = cond_alloc(nalt, F107_IDX_FILE);

  /*cond_f107_override(175.0, 175.0, cond_workspace_p);*/
  /*cond_set_error_scale(4.0, 1.0, cond_workspace_p);*/
  /*cond_set_error_scale(1.0, 4.0, cond_workspace_p);*/

  cond_calc(theta, phi, t, altmin, altstp, nalt, cond_workspace_p);

  /* print results */
  i = 1;
  printf("# Field %zu: altitude (km)\n", i++);
  printf("# Field %zu: electron temperature (K)\n", i++);
  printf("# Field %zu: ion temperature (K)\n", i++);
  printf("# Field %zu: neutral temperature from IRI (K)\n", i++);
  printf("# Field %zu: neutral temperature from MSIS (K)\n", i++);
  printf("# Field %zu: direct conductivity (A/V/m)\n", i++);
  printf("# Field %zu: Pedersen conductivity (A/V/m)\n", i++);
  printf("# Field %zu: Hall conductivity (A/V/m)\n", i++);
  printf("# Field %zu: electron density (cm-3)\n", i++);
  printf("# Field %zu: electron/ion collision frequency (1/s)\n", i++);
  printf("# Field %zu: neutral density (cm-3)\n", i++);

  for (i = 0; i < nalt; ++i)
    {
      double alt = altmin + i * altstp;
      cond_result *result = cond_get_result(i, cond_workspace_p);

      printf("%f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
             alt,
             result->T_e,
             result->T_i,
             result->T_n_iri,
             result->T_n_msis,
             result->sigma_0,
             result->sigma_p,
             result->sigma_h,
             result->n_e,
             result->v_ei,
             result->n_n);
    }

  cond_free(cond_workspace_p);

  return 0;
} /* main() */
