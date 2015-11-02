/*
 * mfield_plot_ext.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include "mageq.h"
#include "mfield.h"
#include "pomme.h"
#include "wmm.h"

int
doplot(const satdata_mag *data, mfield_workspace *w)
{
  size_t i;
  estist_workspace *est_p = estist_alloc(ESTIST_IDX_FILE);
  double rms_f = 0.0;
  double rms_z = 0.0;
  size_t n = 0;

  for (i = 0; i < data->n; ++i)
    {
      double B_ext[4], B_ext_pomme[4];
      double r = data->altitude[i] + data->R;
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double E_st, I_st;
      time_t unix_time = satdata_epoch2timet(data->t[i]);

#if 0
      estist_get(unix_time, &E_st, &I_st, est_p);
      mfield_eval_ext(data->t[i], r, theta, phi, E_st, I_st, B_ext, w);
#endif

      B_ext_pomme[0] = SATDATA_VEC_X(data->B_ext, i);
      B_ext_pomme[1] = SATDATA_VEC_Y(data->B_ext, i);
      B_ext_pomme[2] = SATDATA_VEC_Z(data->B_ext, i);
      B_ext_pomme[3] = gsl_hypot3(B_ext_pomme[0], B_ext_pomme[1], B_ext_pomme[2]);

      rms_f += pow(B_ext[3] - B_ext_pomme[3], 2.0);
      rms_z += pow(B_ext[2] - B_ext_pomme[2], 2.0);
      ++n;

#if 1
      printf("%f %f %f %f %f %f %f\n",
             data->altitude[i],
             data->latitude[i],
             data->longitude[i],
             B_ext[3],
             B_ext_pomme[3],
             B_ext[2],
             B_ext_pomme[2]);
#endif
    }

  rms_f = sqrt(rms_f / n);
  rms_z = sqrt(rms_z / n);

  fprintf(stderr, "RMS F = %f [nT]\n", rms_f);
  fprintf(stderr, "RMS Z = %f [nT]\n", rms_z);

  estist_free(est_p);

  return 0;
}

int
main(int argc, char *argv[])
{
  int c;
  char *coeffile = NULL;
  char *infile = NULL;
  mfield_workspace *mfield_workspace_p = NULL;
  satdata_mag *data;

  while ((c = getopt(argc, argv, "i:c:")) != (-1))
    {
      switch (c)
        {
          case 'c':
            coeffile = optarg;
            break;

          case 'i':
            infile = optarg;
            break;

          default:
            break;
        }
    }

  if (!coeffile || !infile)
    {
      fprintf(stderr, "usage: %s [-c coef_file] [-i dmsp_idx_file]\n", argv[0]);
      exit(1);
    }

  data = satdata_dmsp_read_idx(infile, 1);

  fprintf(stderr, "main: loading coefficients from %s...", coeffile);
  mfield_workspace_p = mfield_read(coeffile);
  fprintf(stderr, "done\n");

  {
    pomme_workspace *pomme_p = pomme_alloc_default();
    mageq_workspace *mageq_p = mageq_alloc();
    double B_ext[4];
    double alt = 0.0;
    double phi = 20.0 * M_PI / 180.0;
    double theta, lat;
    time_t unix_time = 1262304000;

    lat = mageq_calc(phi, alt, unix_time, pomme_p, mageq_p);
    theta = M_PI / 2.0 - lat;

    pomme_calc_ext(theta, phi, unix_time, alt, B_ext, pomme_p);

    fprintf(stderr, "RC0 = %f [nT]\n", B_ext[3]*1.0e9);

    pomme_free(pomme_p);
    mageq_free(mageq_p);
  }

  {
    double r = 6371.2 + 101.0;
    double lat = -36.7 * M_PI / 180.0;
    double theta = M_PI / 2.0 - lat;
    double phi = 203.2 * M_PI / 180.0;
    time_t unix_time = 1262304000;
    double t = satdata_timet2epoch(unix_time);
    double B_ext[4];
    double Est = -1.0;
    double Ist = -2.0;

#if 0
    mfield_eval_ext(t, r, theta, phi, Est, Ist, B_ext, mfield_workspace_p);
    fprintf(stderr, "BX = %f BY = %f BZ = %f F = %f\n",
            B_ext[0],
            B_ext[1],
            B_ext[2],
            B_ext[3]);
    exit(1);
#endif
  }

  doplot(data, mfield_workspace_p);

  mfield_free(mfield_workspace_p);

  return 0;
}
