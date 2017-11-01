/*
 * magdata_list.c
 *
 * This is a high-level wrapper for multiple magdata sources
 * (multiple satellites, ground observatories, etc)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <assert.h>

#include <satdata/satdata.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <common/common.h>

#include "magdata_list.h"

/*
magdata_list_alloc()
  Allocate magdata_list

Inputs: nsources - number of data sources (satellites, observatories)

Return: pointer to workspace
*/

magdata_list *
magdata_list_alloc(const size_t nsources)
{
  magdata_list *w;

  w = calloc(1, sizeof(magdata_list));
  if (!w)
    return 0;

  w->n_tot = nsources;
  w->n = 0;

  w->mdata = calloc(1, nsources * sizeof(magdata *));
  if (!w->mdata)
    {
      magdata_list_free(w);
      return 0;
    }

  return w;
}

void
magdata_list_free(magdata_list *w)
{
  if (w->mdata)
    {
      size_t i;

      for (i = 0; i < w->n; ++i)
        {
          magdata *mdata = w->mdata[i];

          if (mdata)
            magdata_free(mdata);
        }

      free(w->mdata);
    }

  free(w);
}

/*
magdata_list_add
  Read a magdata file and add it to the list of sources

Inputs: filename - magdata file
        w        - workspace

Return: number of data read
*/

size_t
magdata_list_add(const char *filename, magdata_list *w)
{
  if (w->n >= w->n_tot)
    {
      fprintf(stderr, "magdata_list_add: error: only %zu data sources allocated\n", w->n_tot);
      return 0;
    }
  else
    {
      size_t ndata;

      w->mdata[w->n] = magdata_read(filename, NULL);
      if (w->mdata[w->n] == NULL)
        {
          fprintf(stderr, "magdata_list_add: error reading %s\n", filename);
          return 0;
        }

      ndata = w->mdata[w->n]->n;

      ++(w->n);

      return ndata;
    }
}

/*
magdata_list_filter_time()
  Flag any data points outside of [tmin,tmax] with
MAGDATA_FLG_DISCARD

Inputs: tmin - minimum time (decimal year)
        tmax - maximum time (decimal year)
        w    - workspace

Return: number of data flagged

Notes:
1) tmin/tmax can be set to -1 to exclude them from the test
*/

size_t
magdata_list_filter_time(const double tmin, const double tmax,
                         magdata_list *w)
{
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->n_tot; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
          double t = satdata_epoch2year(mptr->t[j]);

          if ((tmin > 0.0 && t < tmin) ||
              (tmax > 0.0 && t > tmax))
            {
              mptr->flags[j] |= MAGDATA_FLG_DISCARD;
              ++cnt;
            }
        }
    }

  return cnt;
}

/*
magdata_list_filter_euler()
  We are not fitting Euler angles (fit_euler is 0),
so discard any data which is marked Euler only

Inputs: w    - workspace

Return: number of data flagged
*/

size_t
magdata_list_filter_euler(magdata_list *w)
{
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->n_tot; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
          if ((mptr->flags[j] & MAGDATA_FLG_FIT_EULER) &&
              !(mptr->flags[j] & MAGDATA_FLG_FIT_MF))
            {
              mptr->flags[j] |= MAGDATA_FLG_DISCARD;
              ++cnt;
            }
        }
    }

  return cnt;
}

int
magdata_list_map(const char *dir_prefix, const magdata_list *w)
{
  int s = 0;
  size_t i;
  char buf[2048];

  for (i = 0; i < w->n_tot; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, w);

      sprintf(buf, "%s/map%zu", dir_prefix, i);

      fprintf(stderr, "magdata_list_map: printing spatial coverage of satellite %zu to %s...",
              i, buf);
      magdata_map(buf, mptr);
      fprintf(stderr, "done\n");
    }

  return s;
}

/*
magdata_list_print()
  Print out all data which will be used for field modeling
*/

int
magdata_list_print(const char *dir_prefix, const gsl_vector *wts_spatial,
                   const magdata_list *w)
{
  int s = 0;
  size_t i, j;
  char buf[2048];
  FILE *fp[12];
  const size_t n = 12; /* number of components to print */
  const char *fmtstr = "%ld %.8f %.4f %.4f %.4f %.4f %.3f %.4f %.4f\n";
  const char *fmtstr_grad = "%ld %.8f %.4f %.4f %.4f %.4f %.3f %.4f %.4f %.4f %.4f\n";
  size_t idx = 0;

  for (i = 0; i < w->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, w);
      size_t k;

      sprintf(buf, "%s/data%zu_X.dat", dir_prefix, i);
      fp[0] = fopen(buf, "w");

      sprintf(buf, "%s/data%zu_Y.dat", dir_prefix, i);
      fp[1] = fopen(buf, "w");

      sprintf(buf, "%s/data%zu_Z.dat", dir_prefix, i);
      fp[2] = fopen(buf, "w");

      sprintf(buf, "%s/data%zu_F.dat", dir_prefix, i);
      fp[3] = fopen(buf, "w");

      sprintf(buf, "%s/data%zu_DX_NS.dat", dir_prefix, i);
      fp[4] = fopen(buf, "w");

      sprintf(buf, "%s/data%zu_DY_NS.dat", dir_prefix, i);
      fp[5] = fopen(buf, "w");

      sprintf(buf, "%s/data%zu_DZ_NS.dat", dir_prefix, i);
      fp[6] = fopen(buf, "w");

      sprintf(buf, "%s/data%zu_DF_NS.dat", dir_prefix, i);
      fp[7] = fopen(buf, "w");

      sprintf(buf, "%s/data%zu_DX_EW.dat", dir_prefix, i);
      fp[8] = fopen(buf, "w");

      sprintf(buf, "%s/data%zu_DY_EW.dat", dir_prefix, i);
      fp[9] = fopen(buf, "w");

      sprintf(buf, "%s/data%zu_DZ_EW.dat", dir_prefix, i);
      fp[10] = fopen(buf, "w");

      sprintf(buf, "%s/data%zu_DF_EW.dat", dir_prefix, i);
      fp[11] = fopen(buf, "w");

      for (j = 0; j < n; ++j)
        {
          if (fp[j] == NULL)
            {
              fprintf(stderr, "magdata_list_print: fp[%zu] is NULL\n", j);
              return -1;
            }
        }

      /* header line */
      fprintf(fp[0], "# X vector data for field modeling (satellite %zu)\n", i);
      fprintf(fp[1], "# Y vector data for field modeling (satellite %zu)\n", i);
      fprintf(fp[2], "# Z vector data for field modeling (satellite %zu)\n", i);
      fprintf(fp[3], "# F scalar data for field modeling (satellite %zu)\n", i);
      fprintf(fp[4], "# DX gradient (N/S) vector data for field modeling (satellite %zu)\n", i);
      fprintf(fp[5], "# DY gradient (N/S) vector data for field modeling (satellite %zu)\n", i);
      fprintf(fp[6], "# DZ gradient (N/S) vector data for field modeling (satellite %zu)\n", i);
      fprintf(fp[7], "# DF gradient (N/S) scalar data for field modeling (satellite %zu)\n", i);
      fprintf(fp[8], "# DX gradient (E/W) vector data for field modeling (satellite %zu)\n", i);
      fprintf(fp[9], "# DY gradient (E/W) vector data for field modeling (satellite %zu)\n", i);
      fprintf(fp[10], "# DZ gradient (E/W) vector data for field modeling (satellite %zu)\n", i);
      fprintf(fp[11], "# DF gradient (E/W) scalar data for field modeling (satellite %zu)\n", i);

      for (j = 0; j < n; ++j)
        {
          k = 1;
          fprintf(fp[j], "# Field %zu: timestamp (UT seconds since 1970-01-01)\n", k++);
          fprintf(fp[j], "# Field %zu: time (decimal year)\n", k++);
          fprintf(fp[j], "# Field %zu: longitude (degrees)\n", k++);
          fprintf(fp[j], "# Field %zu: geocentric latitude (degrees)\n", k++);
          fprintf(fp[j], "# Field %zu: QD latitude (degrees)\n", k++);
          fprintf(fp[j], "# Field %zu: geocentric radius (km)\n", k++);
          fprintf(fp[j], "# Field %zu: spatial weight factor\n", k++);
        }

      fprintf(fp[0], "# Field %zu: X vector measurement (nT)\n", k);
      fprintf(fp[1], "# Field %zu: Y vector measurement (nT)\n", k);
      fprintf(fp[2], "# Field %zu: Z vector measurement (nT)\n", k);
      fprintf(fp[3], "# Field %zu: F scalar measurement (nT)\n", k);
      fprintf(fp[4], "# Field %zu: X vector measurement (nT)\n", k);
      fprintf(fp[5], "# Field %zu: Y vector measurement (nT)\n", k);
      fprintf(fp[6], "# Field %zu: Z vector measurement (nT)\n", k);
      fprintf(fp[7], "# Field %zu: F scalar measurement (nT)\n", k);
      fprintf(fp[8], "# Field %zu: X vector measurement (nT)\n", k);
      fprintf(fp[9], "# Field %zu: Y vector measurement (nT)\n", k);
      fprintf(fp[10], "# Field %zu: Z vector measurement (nT)\n", k);
      fprintf(fp[11], "# Field %zu: F scalar measurement (nT)\n", k);
      ++k;

      fprintf(fp[0], "# Field %zu: X a priori model (nT)\n", k);
      fprintf(fp[1], "# Field %zu: Y a priori model (nT)\n", k);
      fprintf(fp[2], "# Field %zu: Z a priori model (nT)\n", k);
      fprintf(fp[3], "# Field %zu: F a priori model (nT)\n", k);
      fprintf(fp[4], "# Field %zu: X a priori model (nT)\n", k);
      fprintf(fp[5], "# Field %zu: Y a priori model (nT)\n", k);
      fprintf(fp[6], "# Field %zu: Z a priori model (nT)\n", k);
      fprintf(fp[7], "# Field %zu: F a priori model (nT)\n", k);
      fprintf(fp[8], "# Field %zu: X a priori model (nT)\n", k);
      fprintf(fp[9], "# Field %zu: Y a priori model (nT)\n", k);
      fprintf(fp[10], "# Field %zu: Z a priori model (nT)\n", k);
      fprintf(fp[11], "# Field %zu: F a priori model (nT)\n", k);
      ++k;

      fprintf(fp[4], "# Field %zu: X vector measurement at N/S gradient point (nT)\n", k);
      fprintf(fp[5], "# Field %zu: Y vector measurement at N/S gradient point (nT)\n", k);
      fprintf(fp[6], "# Field %zu: Z vector measurement at N/S gradient point (nT)\n", k);
      fprintf(fp[7], "# Field %zu: F scalar measurement at N/S gradient point (nT)\n", k);
      fprintf(fp[8], "# Field %zu: X vector measurement at E/W gradient point (nT)\n", k);
      fprintf(fp[9], "# Field %zu: Y vector measurement at E/W gradient point (nT)\n", k);
      fprintf(fp[10], "# Field %zu: Z vector measurement at E/W gradient point (nT)\n", k);
      fprintf(fp[11], "# Field %zu: F scalar measurement at E/W gradient point (nT)\n", k);
      ++k;

      fprintf(fp[4], "# Field %zu: X a priori model at N/S gradient point (nT)\n", k);
      fprintf(fp[5], "# Field %zu: Y a priori model at N/S gradient point (nT)\n", k);
      fprintf(fp[6], "# Field %zu: Z a priori model at N/S gradient point (nT)\n", k);
      fprintf(fp[7], "# Field %zu: F a priori model at N/S gradient point (nT)\n", k);
      fprintf(fp[8], "# Field %zu: X a priori model at E/W gradient point (nT)\n", k);
      fprintf(fp[9], "# Field %zu: Y a priori model at E/W gradient point (nT)\n", k);
      fprintf(fp[10], "# Field %zu: Z a priori model at E/W gradient point (nT)\n", k);
      fprintf(fp[11], "# Field %zu: F a priori model at E/W gradient point (nT)\n", k);
      ++k;

      for (j = 0; j < mptr->n; ++j)
        {
          double t = satdata_epoch2year(mptr->t[j]);
          time_t unix_time = satdata_epoch2timet(mptr->t[j]);
          double phi = wrap180(mptr->phi[j] * 180.0 / M_PI);
          double lat = 90.0 - mptr->theta[j] * 180.0 / M_PI;
          double qdlat = mptr->qdlat[j];
          double r = mptr->r[j];
          double B[4], B_grad[4];
          double B_model[4], B_grad_model[4];

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          B[0] = mptr->Bx_nec[j];
          B[1] = mptr->By_nec[j];
          B[2] = mptr->Bz_nec[j];
          B[3] = mptr->F[j];

          B_model[0] = mptr->Bx_model[j];
          B_model[1] = mptr->By_model[j];
          B_model[2] = mptr->Bz_model[j];
          B_model[3] = gsl_hypot3(B_model[0], B_model[1], B_model[2]);

          B_grad[0] = mptr->Bx_nec_ns[j];
          B_grad[1] = mptr->By_nec_ns[j];
          B_grad[2] = mptr->Bz_nec_ns[j];
          B_grad[3] = mptr->F_ns[j];

          B_grad_model[0] = mptr->Bx_model_ns[j];
          B_grad_model[1] = mptr->By_model_ns[j];
          B_grad_model[2] = mptr->Bz_model_ns[j];
          B_grad_model[3] = gsl_hypot3(B_grad_model[0], B_grad_model[1], B_grad_model[2]);

          if ((j > 0) && (mptr->flags[j] & MAGDATA_FLG_TRACK_START))
            {
              size_t k;

              for (k = 0; k < n; ++k)
                fprintf(fp[k], "\n\n");
            }

          if (MAGDATA_ExistX(mptr->flags[j]))
            {
              double wj = gsl_vector_get(wts_spatial, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                fprintf(fp[0], fmtstr, unix_time, t, phi, lat, qdlat, r, wj, B[0], B_model[0]);
            }

          if (MAGDATA_ExistY(mptr->flags[j]))
            {
              double wj = gsl_vector_get(wts_spatial, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                fprintf(fp[1], fmtstr, unix_time, t, phi, lat, qdlat, r, wj, B[1], B_model[1]);
            }

          if (MAGDATA_ExistZ(mptr->flags[j]))
            {
              double wj = gsl_vector_get(wts_spatial, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                fprintf(fp[2], fmtstr, unix_time, t, phi, lat, qdlat, r, wj, B[2], B_model[2]);
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            {
              double wj = gsl_vector_get(wts_spatial, idx++);
              fprintf(fp[3], fmtstr, unix_time, t, phi, lat, qdlat, r, wj, B[3], B_model[3]);
            }

          if (MAGDATA_ExistDX_NS(mptr->flags[j]))
            {
              double wj = gsl_vector_get(wts_spatial, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                fprintf(fp[4], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[0], B_model[0], B_grad[0], B_grad_model[0]);
            }

          if (MAGDATA_ExistDY_NS(mptr->flags[j]))
            {
              double wj = gsl_vector_get(wts_spatial, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                fprintf(fp[5], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[1], B_model[1], B_grad[1], B_grad_model[1]);
            }

          if (MAGDATA_ExistDZ_NS(mptr->flags[j]))
            {
              double wj = gsl_vector_get(wts_spatial, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                fprintf(fp[6], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[2], B_model[2], B_grad[2], B_grad_model[2]);
            }

          if (MAGDATA_ExistDF_NS(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            {
              double wj = gsl_vector_get(wts_spatial, idx++);
              fprintf(fp[7], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[3], B_model[3], B_grad[3], B_grad_model[3]);
            }

          if (MAGDATA_ExistDX_EW(mptr->flags[j]))
            {
              double wj = gsl_vector_get(wts_spatial, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                fprintf(fp[8], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[0], B_model[0], B_grad[0], B_grad_model[0]);
            }

          if (MAGDATA_ExistDY_EW(mptr->flags[j]))
            {
              double wj = gsl_vector_get(wts_spatial, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                fprintf(fp[9], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[1], B_model[1], B_grad[1], B_grad_model[1]);
            }

          if (MAGDATA_ExistDZ_EW(mptr->flags[j]))
            {
              double wj = gsl_vector_get(wts_spatial, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                fprintf(fp[10], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[2], B_model[2], B_grad[2], B_grad_model[2]);
            }

          if (MAGDATA_ExistDF_EW(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            {
              double wj = gsl_vector_get(wts_spatial, idx++);
              fprintf(fp[11], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[3], B_model[3], B_grad[3], B_grad_model[3]);
            }
        }

      for (j = 0; j < n; ++j)
        fclose(fp[j]);
    }

  assert(idx == wts_spatial->size);

  return s;
}

magdata *
magdata_list_ptr(const size_t idx, const magdata_list *w)
{
  if (idx >= w->n)
    {
      fprintf(stderr, "magdata_list_ptr: invalid index: %zu\n", idx);
      return 0;
    }

  return w->mdata[idx];
}
