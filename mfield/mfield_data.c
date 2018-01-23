/*
 * mfield_data.c
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

#include "mfield_data.h"

/*
mfield_data_alloc()
  Allocate mfield_data_workspace

Inputs: nsources - number of data sources (satellites,observatories)
        params   - data parameters

Return: pointer to workspace
*/

mfield_data_workspace *
mfield_data_alloc(const size_t nsources, const mfield_data_parameters *params)
{
  mfield_data_workspace *w;

  w = calloc(1, sizeof(mfield_data_workspace));
  if (!w)
    return 0;

  w->nsources = nsources;
  w->params = *params;

  w->t0 = malloc(nsources * sizeof(double));
  w->t1 = malloc(nsources * sizeof(double));
  if (!w->t0 || !w->t1)
    {
      mfield_data_free(w);
      return 0;
    }

  w->mdata = calloc(1, nsources * sizeof(magdata *));
  if (!w->mdata)
    {
      mfield_data_free(w);
      return 0;
    }

  w->rstat_workspace_p = gsl_rstat_alloc();
  if (!w->rstat_workspace_p)
    {
      mfield_data_free(w);
      return 0;
    }

  w->t_mu = -1.0;
  w->t_sigma = -1.0;
  w->t0_data = -1.0;
  w->t1_data = -1.0;

  return w;
} /* mfield_data_alloc() */

void
mfield_data_free(mfield_data_workspace *w)
{
  if (w->t0)
    free(w->t0);

  if (w->t1)
    free(w->t1);

  if (w->rstat_workspace_p)
    gsl_rstat_free(w->rstat_workspace_p);

  if (w->mdata)
    {
      size_t i;

      for (i = 0; i < w->nsources; ++i)
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
mfield_data_filter_time()
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
mfield_data_filter_time(const double tmin, const double tmax,
                        mfield_data_workspace *w)
{
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

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
} /* mfield_data_filter_time() */

/*
mfield_data_filter_euler()
  We are not fitting Euler angles (fit_euler is 0),
so discard any data which is marked Euler only

Inputs: w    - workspace

Return: number of data flagged
*/

size_t
mfield_data_filter_euler(mfield_data_workspace *w)
{
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

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
} /* mfield_data_filter_euler() */

/*
mfield_data_filter_comp()
  Discard any data according to config file flags

Inputs: w    - workspace

Return: number of data flagged
*/

size_t
mfield_data_filter_comp(mfield_data_workspace *w)
{
  const mfield_data_parameters *params = &(w->params);
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
          double qdlat = mptr->qdlat[j];

          if (fabs(qdlat) <= params->qdlat_fit_cutoff)
            {
              /* select components for mid/low latitudes */

              if (!params->fit_X)
                mptr->flags[j] &= ~MAGDATA_FLG_X;

              if (!params->fit_Y)
                mptr->flags[j] &= ~MAGDATA_FLG_Y;

              if (!params->fit_Z)
                mptr->flags[j] &= ~MAGDATA_FLG_Z;

              if (!params->fit_F)
                mptr->flags[j] &= ~MAGDATA_FLG_F;

              if (!params->fit_DX_NS)
                mptr->flags[j] &= ~MAGDATA_FLG_DX_NS;

              if (!params->fit_DY_NS)
                mptr->flags[j] &= ~MAGDATA_FLG_DY_NS;

              if (!params->fit_DZ_NS)
                mptr->flags[j] &= ~MAGDATA_FLG_DZ_NS;

              if (!params->fit_DF_NS)
                mptr->flags[j] &= ~MAGDATA_FLG_DF_NS;

              if (!params->fit_DX_EW)
                mptr->flags[j] &= ~MAGDATA_FLG_DX_EW;

              if (!params->fit_DY_EW)
                mptr->flags[j] &= ~MAGDATA_FLG_DY_EW;

              if (!params->fit_DZ_EW)
                mptr->flags[j] &= ~MAGDATA_FLG_DZ_EW;

              if (!params->fit_DF_EW)
                mptr->flags[j] &= ~MAGDATA_FLG_DF_EW;
            }
          else
            {
              /* select components for high-latitudes */

              /* don't fit X/Y at high-latitudes, including gradients */
              mptr->flags[j] &= ~(MAGDATA_FLG_X | MAGDATA_FLG_Y);
              mptr->flags[j] &= ~(MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS);
              mptr->flags[j] &= ~(MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW);

              if (!params->fit_Z_highlat)
                mptr->flags[j] &= ~MAGDATA_FLG_Z;

              if (!params->fit_F_highlat)
                mptr->flags[j] &= ~MAGDATA_FLG_F;

              if (!params->fit_DZ_NS_highlat)
                mptr->flags[j] &= ~MAGDATA_FLG_DZ_NS;

              if (!params->fit_DF_NS_highlat)
                mptr->flags[j] &= ~MAGDATA_FLG_DF_NS;

              if (!params->fit_DZ_EW_highlat)
                mptr->flags[j] &= ~MAGDATA_FLG_DZ_EW;

              if (!params->fit_DF_EW_highlat)
                mptr->flags[j] &= ~MAGDATA_FLG_DF_EW;
            }
        }
    }

  return cnt;
}

/*
mfield_data_init()
  Compute mean and stddev of timestamps minus epoch for later time scaling

w_i = t_i - epoch
t_mu = mean(w_i)
t_sigma = stddev(w_i)

Inputs: w - workspace

Return: success/error

Notes:
1) w->t_mu and w->t_sigma are updated with timestamp mean/stddev in years

2) w->t0_data is initialized to the timestamp of the first data point (CDF_EPOCH)

3) w->t1_data is initialized to the timestamp of the last data point (CDF_EPOCH)

3) w->t0 and w->t1 are initialized to the first/last timestamps of each satellite
*/

int
mfield_data_init(mfield_data_workspace *w)
{
  int s = 0;
  size_t i, j;

  gsl_rstat_reset(w->rstat_workspace_p);

  w->t0_data = 1.0e15;
  w->t1_data = -1.0e15;
  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      magdata_t(&(w->t0[i]), &(w->t1[i]), mptr);

      if (mptr->n != 0)
        {
          if (w->t0[i] > 0.0)
            w->t0_data = GSL_MIN(w->t0_data, w->t0[i]);

          if (w->t1[i] > 0.0)
            w->t1_data = GSL_MAX(w->t1_data, w->t1[i]);
        }

      for (j = 0; j < mptr->n; ++j)
        {
          double t;

          if (mptr->flags[j] & MAGDATA_FLG_DISCARD)
            continue;

          t = satdata_epoch2year(mptr->t[j]) - w->params.epoch;
          gsl_rstat_add(t, w->rstat_workspace_p);
        }
    }

  w->t_mu = gsl_rstat_mean(w->rstat_workspace_p);
  w->t_sigma = gsl_rstat_sd(w->rstat_workspace_p);

  if (w->t_sigma == 0.0)
    {
      /* this can happen for a fixed time grid like EMAG2 */
      w->t_mu = 0.0;
      w->t_sigma = 1.0;
    }

  return s;
} /* mfield_data_init() */

/*
mfield_data_epoch()
  Compute epoch of input data by averaging all timestamps
*/

double
mfield_data_epoch(mfield_data_workspace *w)
{
  /* initialize t_mu and t_sigma */
  mfield_data_init(w);

  return w->t_mu + w->params.epoch;
} /* mfield_data_epoch() */

int
mfield_data_map(const char *dir_prefix, const mfield_data_workspace *w)
{
  int s = 0;
  size_t i;
  char buf[2048];

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      sprintf(buf, "%s/datamap%zu", dir_prefix, i);

      fprintf(stderr, "mfield_data_map: printing spatial coverage of satellite %zu to %s...",
              i, buf);
      magdata_map(buf, mptr);
      fprintf(stderr, "done\n");
    }

  return s;
}

/*
mfield_data_print()
  Print out all data which will be used for main field modeling
*/

int
mfield_data_print(const char *dir_prefix, const gsl_vector *wts_spatial,
                  const mfield_data_workspace *w)
{
  int s = 0;
  size_t i, j;
  char buf[2048];
  FILE *fp[12];
  const size_t n = 12; /* number of components to print */
  const char *fmtstr = "%ld %.8f %.4f %.4f %.4f %.4f %.3f %.4f %.4f\n";
  const char *fmtstr_grad = "%ld %.8f %.4f %.4f %.4f %.4f %.3f %.4f %.4f %.4f %.4f\n";
  size_t idx = 0;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);
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
              fprintf(stderr, "mfield_data_print: fp[%zu] is NULL\n", j);
              return -1;
            }
        }

      /* header line */
      fprintf(fp[0], "# X vector data for MF modeling (satellite %zu)\n", i);
      fprintf(fp[1], "# Y vector data for MF modeling (satellite %zu)\n", i);
      fprintf(fp[2], "# Z vector data for MF modeling (satellite %zu)\n", i);
      fprintf(fp[3], "# F scalar data for MF modeling (satellite %zu)\n", i);
      fprintf(fp[4], "# DX gradient (N/S) vector data for MF modeling (satellite %zu)\n", i);
      fprintf(fp[5], "# DY gradient (N/S) vector data for MF modeling (satellite %zu)\n", i);
      fprintf(fp[6], "# DZ gradient (N/S) vector data for MF modeling (satellite %zu)\n", i);
      fprintf(fp[7], "# DF gradient (N/S) scalar data for MF modeling (satellite %zu)\n", i);
      fprintf(fp[8], "# DX gradient (E/W) vector data for MF modeling (satellite %zu)\n", i);
      fprintf(fp[9], "# DY gradient (E/W) vector data for MF modeling (satellite %zu)\n", i);
      fprintf(fp[10], "# DZ gradient (E/W) vector data for MF modeling (satellite %zu)\n", i);
      fprintf(fp[11], "# DF gradient (E/W) scalar data for MF modeling (satellite %zu)\n", i);

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
mfield_data_ptr(const size_t idx, const mfield_data_workspace *w)
{
  if (idx >= w->nsources)
    {
      fprintf(stderr, "mfield_data_ptr: invalid index: %zu\n", idx);
      return 0;
    }

  return w->mdata[idx];
}

/*
mfield_data_add_noise()
  For WMM simulation studies - add random Gaussian noise to vector
data

Inputs: sigma - standard deviation of noise (nT);
                if negative, not used
        bias  - bias for vector components in VFM frame (nT)
        w     - workspace
*/

int
mfield_data_add_noise(const double sigma, const double bias, mfield_data_workspace * w)
{
  int s = 0;
  size_t i, j;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* if no vector measurement, do nothing */
          if (MAGDATA_ExistVector(mptr->flags[j]))
            {
              if (sigma > 0.0)
                {
                  mptr->Bx_vfm[j] += gsl_ran_gaussian(r, sigma);
                  mptr->By_vfm[j] += gsl_ran_gaussian(r, sigma);
                  mptr->Bz_vfm[j] += gsl_ran_gaussian(r, sigma);
                }

              mptr->Bx_vfm[j] += bias;
              mptr->By_vfm[j] += bias;
              mptr->Bz_vfm[j] += bias;

              /* recompute scalar field measurement */
              mptr->F[j] = gsl_hypot3(mptr->Bx_vfm[j], mptr->By_vfm[j], mptr->Bz_vfm[j]);
            }
          else if (MAGDATA_ExistScalar(mptr->flags[j]))
            {
              /* scalar only measurement (high latitudes) */
              mptr->F[j] += bias + gsl_ran_gaussian(r, sigma);
            }
        }
    }

  gsl_rng_free(r);

  return s;
}
