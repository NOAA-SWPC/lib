/*
 * test.c
 *
 * Test IRI module
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include "iri.h"

int
main(void)
{
  const char *data_file = "data.dat";
  FILE *fp;
  size_t i;
  time_t t;
  double theta, phi;
  double altmin, altmax, altstp;
  size_t nalt;
  iri_workspace *iri_workspace_p;

#if 0
  /* Mar 20, 2005 */
  t = 1111276800;
  /*t = 1325401200;*/
  t = 1388622455; /* Jan 2 2014 */
  theta = M_PI / 2.0 - 10.0 * M_PI / 180.0;
  phi = 15.0 * M_PI / 180.0;

  altmin = 50.0;
  altstp = 2.0;
  nalt = 300;
#else
  t = 1032048000; /* Sep 15 2002 0UT */
  theta = 75.0 * M_PI / 180.0;
  phi = 225.0 * M_PI / 180.0;

  altmin = 0.0;
  altstp = 1.0;
  nalt = 1000;
#endif

  altmax = altmin + nalt * altstp;

  iri_workspace_p = iri_alloc(nalt, F107_IDX_FILE);

  /*iri_f107_override(30075.0, iri_workspace_p);*/

  iri_calc(theta, phi, t, altmin, altstp, nalt, iri_workspace_p);

  fp = fopen(data_file, "w");

  fprintf(stderr, "main: writing data to %s...", data_file);

  /* print results */
  i = 1;
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: electron density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: neutral temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: ion temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: electron temperature (K)\n", i++);
  fprintf(fp, "# Field %zu: O+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: H+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: HE+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: O2+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: NO+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: N+ density (m^-3)\n", i++);
  fprintf(fp, "# Field %zu: cluster ion density (m^-3)\n", i++);

  for (i = 0; i < nalt; ++i)
    {
      double alt = altmin + i * altstp;
      iri_result *result = iri_get_result(i, iri_workspace_p);

      fprintf(fp, "%f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
              alt,
              result->Ne,
              result->Tn,
              result->Ti,
              result->Te,
              result->n_Op,
              result->n_Hp,
              result->n_HEp,
              result->n_O2p,
              result->n_NOp,
              result->n_Np,
              result->n_clus);
    }

  fprintf(stderr, "done\n");

  iri_free(iri_workspace_p);
  fclose(fp);

  return 0;
} /* main() */
