/*
 * mag2swarm.c
 *
 * Convert a magdata file to Swarm CDF format
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

#include <satdata/satdata.h>

#include "common.h"
#include "magdata.h"
#include "pomme.h"

/*
replace_ephemeris()
  For Swarm crustal field study - subtract 100 km from all
data points, recalculate main field and crustal field model
at new location
*/

int
replace_ephemeris(magdata *mdata)
{
  const size_t n = mdata->n;
  size_t max_threads = (size_t) omp_get_max_threads();
  size_t i;
  size_t *omp_n;
  msynth_workspace *msynth_core[24];
  msynth_workspace *msynth_mf7[24];
  msynth_workspace *msynth_newcrust[24];
  pomme_workspace *pomme_p = pomme_alloc_default();

  pomme_set_radius(R_EARTH_KM, pomme_p);
  omp_n = malloc(max_threads * sizeof(size_t));

  for (i = 0; i < max_threads; ++i)
    {
      omp_n[i] = 0;
      msynth_core[i] = msynth_read(MSYNTH_BOUMME_FILE);
      msynth_mf7[i] = msynth_mf7_read(MSYNTH_MF7_FILE);
      msynth_newcrust[i] = msynth_ngdc720_read("new_coef.txt");

      msynth_set(1, 15, msynth_core[i]);
      msynth_set(16, 133, msynth_mf7[i]);
      msynth_set(16, 300, msynth_newcrust[i]);
    }

  fprintf(stderr, "\n");

#pragma omp parallel for private(i)
  for (i = 0; i < n; ++i)
    {
      int thread_id = omp_get_thread_num();
      msynth_workspace *core_p = msynth_core[thread_id];
      msynth_workspace *mf7_p = msynth_mf7[thread_id];
      msynth_workspace *newcrust_p = msynth_newcrust[thread_id];
      double tyr = satdata_epoch2year(mdata->t[i]);
      time_t unix_time = satdata_epoch2timet(mdata->t[i]);
      double r = mdata->r[i];
      double theta = mdata->theta[i];
      double phi = mdata->phi[i];
      double alt = r - mdata->R;
      double B_core[4], B_crust[4], B_ext[4], B_tot[4];
      int j;

      ++omp_n[thread_id];

      /* evaluate field models at current location */
      msynth_eval(tyr, r, theta, phi, B_core, core_p);
      msynth_eval(tyr, r, theta, phi, B_crust, mf7_p);

      pomme_calc_ext(theta, phi, unix_time, alt, B_ext, pomme_p);
      for (j = 0; j < 3; ++j)
        B_ext[j] *= 1.0e9;

      for (j = 0; j < 3; ++j)
        {
          B_tot[j] = B_core[j] + B_crust[j] + B_ext[j];
        }

      /* subtract field models from current location */
      mdata->Bx_nec[i] -= B_tot[0];
      mdata->By_nec[i] -= B_tot[1];
      mdata->Bz_nec[i] -= B_tot[2];

      /* reduce radius by 100 km */
      r -= 100.0;
      alt -= 100.0;
      mdata->r[i] = r;

      /* evaluate field models at new location (with new crustal field) */
      msynth_eval(tyr, r, theta, phi, B_core, core_p);
      msynth_eval(tyr, r, theta, phi, B_crust, newcrust_p);

      pomme_calc_ext(theta, phi, unix_time, alt, B_ext, pomme_p);
      for (j = 0; j < 3; ++j)
        B_ext[j] *= 1.0e9;

      for (j = 0; j < 3; ++j)
        {
          B_tot[j] = B_core[j] + B_crust[j] + B_ext[j];
        }

      /* add back field models at new location */
      mdata->Bx_nec[i] += B_tot[0];
      mdata->By_nec[i] += B_tot[1];
      mdata->Bz_nec[i] += B_tot[2];

      if (thread_id == 0 && omp_n[0] % 100 == 0)
        {
          double progress = 0.0;

          for (j = 0; j < max_threads; ++j)
            progress += (double) omp_n[j];
          progress /= (double) n;

          progress_bar(stderr, progress, 70);
        }
    }

  for (i = 0; i < max_threads; ++i)
    {
      msynth_free(msynth_core[i]);
      msynth_free(msynth_mf7[i]);
      msynth_free(msynth_newcrust[i]);
    }

  pomme_free(pomme_p);
  free(omp_n);

  return 0;
}

int
main(int argc, char *argv[])
{
  char *swarm_file = NULL;
  magdata *mdata = NULL;
  satdata_mag *data;
  int new_ephemeris = 0;
  int c;

  while ((c = getopt(argc, argv, "i:o:n")) != (-1))
    {
      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            mdata = magdata_read(optarg, NULL);
            fprintf(stderr, "done (%zu data read)\n", mdata->n);
            break;

          case 'o':
            swarm_file = optarg;
            break;

          case 'n':
            new_ephemeris = 1;
            break;

          default:
            break;
        }
    }

  if (mdata == NULL || swarm_file == NULL)
    {
      fprintf(stderr, "Usage: %s <-i magdata_file> <-o swarm_file> [-n]\n", argv[0]);
      exit(1);
    }

  if (new_ephemeris)
    {
      fprintf(stderr, "main: adjusting ephemeris to lower altitude...");
      replace_ephemeris(mdata);
      fprintf(stderr, "done\n");
    }

  data = magdata_mag2sat(mdata);

  fprintf(stderr, "main: writing Swarm CDF file %s...", swarm_file);
  satdata_swarm_write(0, swarm_file, data);
  fprintf(stderr, "done (%zu data written)\n", data->n);

  magdata_free(mdata);
  satdata_mag_free(data);

  return 0;
}
