/*
 * secs1d.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <satdata/satdata.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>

#include "common.h"
#include "secs1d.h"
#include "track.h"

static int build_matrix_row(const double r, const double theta, const double phi,
                            gsl_vector *X, gsl_vector *Y, gsl_vector *Z,
                            secs1d_workspace *w);

/*
secs1d_alloc()
  Allocate secs 1d workspace

Inputs: lmax   - maximum degree for Legendre functions in expansion
        R_iono - radius of ionosphere (km)
        npoles - number of SECS poles

Return: pointer to workspace
*/

secs1d_workspace *
secs1d_alloc(const size_t lmax, const double R_iono, const size_t npoles)
{
  secs1d_workspace *w;

  w = calloc(1, sizeof(secs1d_workspace));
  if (!w)
    return 0;

  w->nmax = 30000;
  w->n = 0;
  w->p = npoles;
  w->R_iono = R_iono;
  w->lmax = lmax;

  w->X = gsl_matrix_alloc(w->nmax, w->p);
  w->c = gsl_vector_alloc(w->p);
  w->rhs = gsl_vector_alloc(w->nmax);
  w->wts = gsl_vector_alloc(w->nmax);

  w->theta0 = malloc(w->p * sizeof(double));
  w->Pltheta0 = malloc((w->lmax + 1) * sizeof(double));
  w->Pltheta = malloc((w->lmax + 1) * sizeof(double));

  return w;
}

void
secs1d_free(secs1d_workspace *w)
{
  if (w->X)
    gsl_matrix_free(w->X);

  if (w->c)
    gsl_vector_free(w->c);

  if (w->rhs)
    gsl_vector_free(w->rhs);

  if (w->wts)
    gsl_vector_free(w->wts);

  if (w->theta0)
    free(w->theta0);

  if (w->Pltheta0)
    free(w->Pltheta0);

  if (w->Pltheta)
    free(w->Pltheta);

  free(w);
}

int
secs1d_add_track(const track_data *tptr, const satdata_mag *data,
                 track_workspace *track_p, secs1d_workspace *w)
{
  size_t rowidx = w->n;
  size_t i;

  for (i = 0; i < tptr->n; ++i)
    {
      size_t didx = i + tptr->start_idx;
      double r = data->altitude[didx] + data->R;
      double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
      double phi = data->longitude[didx] * M_PI / 180.0;
      double B_res[3];
      double wi = 1.0;
      gsl_vector_view vx = gsl_matrix_row(w->X, rowidx);
      gsl_vector_view vy = gsl_matrix_row(w->X, rowidx + 1);
      gsl_vector_view vz = gsl_matrix_row(w->X, rowidx + 2);

      if (SATDATA_BadData(data->flags[didx]))
        continue;

      /* fit only low-latitude data */
      if (fabs(data->qdlat[didx]) > 50.0)
        continue;

      /* compute magnetic residual */
      track_residual(i, didx, B_res, data, track_p);

      /* set rhs vector */
      gsl_vector_set(w->rhs, rowidx, B_res[0]);
      gsl_vector_set(w->rhs, rowidx + 1, B_res[1]);
      gsl_vector_set(w->rhs, rowidx + 2, B_res[2]);

      /* set weight vector */
      gsl_vector_set(w->wts, rowidx, wi);
      gsl_vector_set(w->wts, rowidx + 1, wi);
      gsl_vector_set(w->wts, rowidx + 2, wi);

      /* build 3 rows of the LS matrix */
      build_matrix_row(r, theta, phi, &vx.vector, &vy.vector, &vz.vector, w);

      rowidx += 3;
    }

  w->n = rowidx;

  return 0;
}

/*
secs1d_green()
  Compute magnetic field Green's function for a single divergence-free
1D SECS

Inputs: r      - radius (km)
        theta  - colatitude (radians)
        theta0 - pole position (radians)
        B      - (output) magnetic field Green's function (X,Y,Z)
        w      - workspace

Return: success/error
*/

int
secs1d_df_green(const double r, const double theta, const double theta0,
                double B[3], secs1d_workspace *w)
{
  const double costheta = cos(theta);
  B[0] = 0.0;
  B[1] = 0.0;
  B[2] = 0.0;

  /* compute P_l(cos(theta)) */
  gsl_sf_legendre_Pl_array(w->lmax, costheta, w->Pltheta);

  /* compute P_l(cos(theta0)) */
  gsl_sf_legendre_Pl_array(w->lmax, cos(theta0), w->Pltheta0);

  if (r > w->R_iono)
    {
      double ratio = w->R_iono / r;
      double rterm = ratio; /* (R / r)^{l+1} */
      size_t l;

      for (l = 1; l <= w->lmax; ++l)
        {
          double Pl1 = gsl_sf_legendre_Plm(l, 1, costheta);

          /* (R / r)^{l+1} */
          rterm *= ratio;

          B[0] += rterm / (l + 1.0) * w->Pltheta0[l] * Pl1;
          B[2] -= rterm * w->Pltheta0[l] * w->Pltheta[l];
        }

      B[0] *= SECS1D_MU_0 / (2.0 * r);
      B[2] *= SECS1D_MU_0 / (2.0 * r);
    }
  else
    {
    }

  return GSL_SUCCESS;
}

/*
secs1d_df_green_J()
  Compute current density Green's function for a single divergence-free
1D SECS

Inputs: theta  - colatitude (radians)
        theta0 - pole position (radians)
        K      - (output) current density Green's function (X,Y,Z)
        w      - workspace

Return: success/error
*/

int
secs1d_df_green_J(const double theta, const double theta0, double K[3], secs1d_workspace *w)
{
  K[0] = 0.0;
  K[2] = 0.0;

  if (theta < theta0)
    {
      K[1] = -tan(0.5 * theta) / (2.0 * w->R_iono);
    }
  else
    {
      K[1] = 1.0 / tan(0.5 * theta) / (2.0 * w->R_iono);
    }

  return GSL_SUCCESS;
}

static int
build_matrix_row(const double r, const double theta, const double phi,
                 gsl_vector *X, gsl_vector *Y, gsl_vector *Z,
                 secs1d_workspace *w)
{
  const size_t p = X->size;
  size_t i;
  double B[3] = { 0.0, 0.0, 0.0 };

  /* compute P_l(cos(theta)) */
  gsl_sf_legendre_Pl_array(w->lmax, cos(theta), w->Pltheta);

#if 0
  if (r > w->R_iono)
    {
      double ratio = w->R_iono / r;

      for (i = 0; i < p; ++i)
        {
          double theta0 = w->theta0[i];
          double rterm = ratio; /* (R / r)^{l+1} */
          size_t l;

          /* compute P_l(cos(theta0)) */
          gsl_sf_legendre_Pl_array(w->lmax, cos(theta0), w->Pltheta0);

          for (l = 1; l <= w->lmax; ++l)
            {
              /* (R / r)^{l+1} */
              rterm *= ratio;

              B[1] += rterm / (l + 1.0) * w->Pltheta0[l] * Pl1;
              B[2] -= rterm * w->Pltheta0[l] * w->Pltheta[l];
            }

          B[2] *= SECS1D_MU_0 / (2.0 * r);

          /* field from curl-free current */
          if (theta > theta0)
            B[1] = -SECS1D_MU_0 / (2.0 * r) * tan(0.5 * (M_PI - theta));
          else
            B[1] = SECS1D_MU_0 / (2.0 * r) * tan(0.5 * theta);

          gsl_vector_set(X, i, B[0]);
          gsl_vector_set(Y, i, B[1]);
          gsl_vector_set(Z, i, B[2]);
        }
    }
  else
    {
    }
#endif

  return 0;
}
