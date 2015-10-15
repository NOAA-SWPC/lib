/*
 * test_omni.c
 *
 * Test the OMNI module by reading all data files, testing a few
 * specific known values, and checking for missing data
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include "omni.h"

int
main(int argc, char *argv[])
{
  int s = 0;
  omni_workspace *omni_workspace_p;
  time_t t0_2011 = 1293840000;
  time_t t0 = 852076800;
  time_t t_now;
  time_t t;
  double B[4];
  double V;

  omni_workspace_p = omni_alloc(OMNI_DATA_DIR);
  time(&t_now);

  /* test some specific values */

  s += omni_get(821984400, B, &V, omni_workspace_p);
  gsl_test_rel(B[0], 1.9, 1.0e-10, "Bx");
  gsl_test_rel(B[1], -3.2, 1.0e-10, "By");
  gsl_test_rel(B[2], -1.0, 1.0e-10, "Bz");
  gsl_test_rel(V, 482.0, 1.0e-10, "V");

  for (t = t0; t < t0_2011; t += 86400)
    {
      s += omni_get(t, B, &V, omni_workspace_p);
      /*printf("%s %f\n", ctime(&t), omni);*/
    }

  omni_free(omni_workspace_p);

  fprintf(stderr, "%s: %d errors\n", argv[0], s);

  return s;
} /* main() */
