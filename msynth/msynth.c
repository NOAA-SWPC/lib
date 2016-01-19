/*
 * msynth.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>

#include "msynth.h"

static int msynth_eval_g(const double t, const double r, const double theta,
                         const double phi, const double epoch, const double *g,
                         const double *dg, const double *ddg, double B[4],
                         msynth_workspace *w);
static int msynth_eval_dgdt(const double t, const double r, const double theta,
                            const double phi, const double epoch,
                            const double *dg, const double *ddg,
                            double dBdt[4], msynth_workspace *w);
static int msynth_compare(const void *a, const void *b);

/*
msynth_alloc()
  Allocate an msynth workspace

Inputs: nmax      - maximum spherical harmonic degree
        nsnapshot - number of snapshot models
        epochs    - array of size nsnapshot containing model epochs
                    (may be NULL if user wishes to fill in later)
*/

msynth_workspace *
msynth_alloc(const size_t nmax, const size_t nsnapshot, const double *epochs)
{
  msynth_workspace *w;
  size_t nnm_tot;
  size_t plm_array_size = gsl_sf_legendre_array_n(nmax);
  size_t i;

  w = calloc(1, sizeof(msynth_workspace));
  if (!w)
    return 0;

  w->nmax = nmax;
  w->nmax_ext = 2;
  w->R = 6371.2;

  /* exclude (0,0) coefficient */
  nnm_tot = (nmax + 1) * (nmax + 1);
  w->nnm = nnm_tot - 1;
  w->p = 3 * w->nnm; /* include SV,SA coefficients */

  nnm_tot = (w->nmax_ext + 1) * (w->nmax_ext + 1);
  w->nnm_ext = nnm_tot - 1;

  w->sv_offset = w->nnm;
  w->sa_offset = 2 * w->nnm;

  w->n_snapshot = nsnapshot;

  w->c = malloc(w->n_snapshot * w->p * sizeof(double));
  w->epochs = malloc(w->n_snapshot * sizeof(double));

  /* fill in epochs */
  if (epochs)
    {
      for (i = 0; i < nsnapshot; ++i)
        w->epochs[i] = epochs[i];

      w->n_epochs = w->n_snapshot;
    }

  w->cosmphi = malloc((nmax + 1) * sizeof(double));
  w->sinmphi = malloc((nmax + 1) * sizeof(double));

  w->Plm = malloc(plm_array_size * sizeof(double));
  w->dPlm = malloc(plm_array_size * sizeof(double));
  if (!w->Plm || !w->dPlm)
    {
      msynth_free(w);
      return 0;
    }

  w->dX = malloc(w->nnm * sizeof(double));
  w->dY = malloc(w->nnm * sizeof(double));
  w->dZ = malloc(w->nnm * sizeof(double));

  w->dX_ext = malloc(w->nnm_ext * sizeof(double));
  w->dY_ext = malloc(w->nnm_ext * sizeof(double));
  w->dZ_ext = malloc(w->nnm_ext * sizeof(double));

  w->eval_nmin = 1;
  w->eval_nmax = nmax;

  w->data_start = -1.0;
  w->data_end = -1.0;

  /* initialize coefficient array to 0 */
  for (i = 0; i < w->n_snapshot * w->p; ++i)
    w->c[i] = 0.0;

  return w;
} /* msynth_alloc() */

void
msynth_free(msynth_workspace *w)
{
  if (w->c)
    free(w->c);

  if (w->epochs)
    free(w->epochs);

  if (w->cosmphi)
    free(w->cosmphi);

  if (w->sinmphi)
    free(w->sinmphi);

  if (w->Plm)
    free(w->Plm);

  if (w->dPlm)
    free(w->dPlm);

  if (w->dX)
    free(w->dX);

  if (w->dY)
    free(w->dY);

  if (w->dZ)
    free(w->dZ);

  if (w->dX_ext)
    free(w->dX_ext);

  if (w->dY_ext)
    free(w->dY_ext);

  if (w->dZ_ext)
    free(w->dZ_ext);

  free(w);
} /* msynth_free() */

int
msynth_set(const size_t nmin, const size_t nmax, msynth_workspace *w)
{
  w->eval_nmin = GSL_MAX(nmin, 1);
  w->eval_nmax = GSL_MIN(nmax, w->nmax);

  return 0;
}

int
msynth_set_data_interval(const double data_start, const double data_end,
                         msynth_workspace *w)
{
  w->data_start = data_start;
  w->data_end = data_end;
  return 0;
}

/*
msynth_eval()
  Evaluate magnetic field model at given point

Inputs: t     - timestamp (decimal years)
        r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - geocentric longitude (radians)
        B     - (output) magnetic field (nT)
                B[0] = B_x
                B[1] = B_y
                B[2] = B_z
                B[3] = |B|
        w     - workspace
*/

int
msynth_eval(const double t, const double r, const double theta,
            const double phi, double B[4], msynth_workspace *w)
{
  int s = 0;
  const size_t epoch_idx = msynth_epoch_idx(t, w);
  const double epoch = w->epochs[epoch_idx];
  const double *g = w->c + epoch_idx * w->p;
  const double *dg = g + w->sv_offset;
  const double *ddg = g + w->sa_offset;

  s = msynth_eval_g(t, r, theta, phi, epoch, g, dg, ddg, B, w);

  return s;
} /* msynth_eval() */

/*
msynth_eval_dBdt()
  Evaluate magnetic field model at given point

Inputs: t     - timestamp (decimal years)
        r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - geocentric longitude (radians)
        dBdt  - (output) magnetic field SV (nT/year)
                dBdt[0] = d/dt B_x
                dBdt[1] = d/dt B_y
                dBdt[2] = d/dt B_z
                dBdt[3] = |dBdt|
        w     - workspace
*/

int
msynth_eval_dBdt(const double t, const double r, const double theta,
                 const double phi, double dBdt[4], msynth_workspace *w)
{
  int s = 0;
  const size_t epoch_idx = msynth_epoch_idx(t, w);
  const double epoch = w->epochs[epoch_idx];
  const double *g = w->c + epoch_idx * w->p;
  const double *dg = g + w->sv_offset;
  const double *ddg = g + w->sa_offset;

  s = msynth_eval_dgdt(t, r, theta, phi, epoch, dg, ddg, dBdt, w);

  return s;
} /* msynth_eval_dBdt() */

/*
msynth_copy()
  Make a copy of an msynth workspace that can be used
interchangeably with msynth_eval(). This is useful for
OpenMP applications, since msynth_eval() uses internal
arrays (dX,dY,dZ) so we need separate workspaces for each
thread for evaluation
*/

msynth_workspace *
msynth_copy(const msynth_workspace *w)
{
  msynth_workspace *w_copy;

  w_copy = msynth_alloc(w->nmax, w->n_snapshot, w->epochs);
  if (!w_copy)
    return 0;

  /* copy coefficients */
  memcpy(w_copy->c, w->c, w->n_snapshot * w->p * sizeof(double));

  w_copy->eval_nmin = w->eval_nmin;
  w_copy->eval_nmax = w->eval_nmax;

  return w_copy;
} /* msynth_copy() */

/*
msynth_diff()
  Create an msynth workspace which is a difference of
two input workspaces, so that:

g = g1 - g2
dg = dg1 - dg2
ddg = ddg1 - ddg2

Inputs: epoch   - epoch at which to compute difference
        msynth1 - workspace 1
        msynth2 - workspace 2
*/

msynth_workspace *
msynth_diff(const double epoch, const msynth_workspace *msynth1,
            const msynth_workspace *msynth2)
{
  msynth_workspace *w1 = msynth_copy(msynth1);
  msynth_workspace *w2 = msynth_copy(msynth2);
  msynth_workspace *w;
  const size_t nmax = GSL_MIN(w1->nmax, w2->nmax);
  size_t n;

  w = msynth_alloc(nmax, 1, &epoch);
  if (!w)
    return 0;

  /* extrapolate both workspaces to epoch */
  msynth_extrapolate_g(epoch, w1);
  msynth_extrapolate_g(epoch, w2);

  for (n = 1; n <= w->nmax; ++n)
    {
      int ni = (int) n;
      int m;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = msynth_nmidx(n, m, w);
          double gnm = msynth_get_mf(epoch, cidx, w1) -
                       msynth_get_mf(epoch, cidx, w2);
          double dgnm = msynth_get_sv(epoch, cidx, w1) -
                        msynth_get_sv(epoch, cidx, w2);
          double ddgnm = msynth_get_sa(epoch, cidx, w1) -
                         msynth_get_sa(epoch, cidx, w2);

          w->c[cidx] = gnm;
          w->c[cidx + w->sv_offset] = dgnm;
          w->c[cidx + w->sa_offset] = ddgnm;
        }
    }

  msynth_free(w1);
  msynth_free(w2);

  return w;
} /* msynth_diff() */

/*
msynth_extrapolate_g()
  Extrapolate Gauss coefficients to a different time

Inputs: t    - new epoch (decimal years)
        w    - workspace

Return: success/error

Notes:
1) w->c is updated with new coefficients
2) w->epochs is updated with new epoch
3) if 't' falls outside previously specified data
   interval, no SA is used, and the SV is set to
   the instantaneous SV at the edge of the interval
*/

int
msynth_extrapolate_g(const double t, msynth_workspace *w)
{
  int s = 0;
  const size_t epoch_idx = msynth_epoch_idx(t, w);
  double *g = w->c + epoch_idx * w->p;
  double *dg = g + w->sv_offset;
  double *ddg = g + w->sa_offset;
  const double t0 = w->epochs[epoch_idx];
  const double t1 = t - t0;
  size_t n;

  if ((w->data_end > 0.0) && (t > w->data_end))
    {
      const double t2 = w->data_end - t0;
      const double t3 = t - w->data_end;

      /*
       * t is outside (and after) the data interval, use:
       * gnm(t) = gnm(t_end) + d/dt gnm(t_end) (t - t_end)
       */
      for (n = 1; n <= w->nmax; ++n)
        {
          int ni = (int) n;
          int m;

          for (m = -ni; m <= ni; ++m)
            {
              size_t cidx = msynth_nmidx(n, m, w);
              double gnm = g[cidx] + t1 * dg[cidx] +
                           ddg[cidx] * t2 * (0.5*t2 + t3);
              double dgnm = dg[cidx] + t2 * ddg[cidx];

              g[cidx] = gnm;
              dg[cidx] = dgnm;
              ddg[cidx] = 0.0;
            }
        }
    }
  else if ((w->data_start > 0.0) && (t < w->data_start))
    {
      const double t2 = w->data_start - t0;
      const double t3 = w->data_start - t;

      /*
       * t is outside (and before) the data interval, use:
       * gnm(t) = gnm(t_start) + d/dt gnm(t_start) (t_start - t)
       */
      for (n = 1; n <= w->nmax; ++n)
        {
          int ni = (int) n;
          int m;

          for (m = -ni; m <= ni; ++m)
            {
              size_t cidx = msynth_nmidx(n, m, w);
              double gnm = g[cidx] + dg[cidx] * (t2 + t3) +
                           ddg[cidx] * t2 * (0.5*t2 + t3);
              double dgnm = -(dg[cidx] + t2 * ddg[cidx]);

              g[cidx] = gnm;
              dg[cidx] = dgnm;
              ddg[cidx] = 0.0;
            }
        }
    }
  else
    {
      const double t2 = 0.5 * t1 * t1;

      /* t is in the data interval, use normal 2nd order taylor series */
      for (n = 1; n <= w->nmax; ++n)
        {
          int ni = (int) n;
          int m;

          for (m = -ni; m <= ni; ++m)
            {
              size_t cidx = msynth_nmidx(n, m, w);
              double gnm = g[cidx] + t1 * dg[cidx] + t2 * ddg[cidx];
              double dgnm = dg[cidx] + t1 * ddg[cidx];
              double ddgnm = ddg[cidx];

              g[cidx] = gnm;
              dg[cidx] = dgnm;
              ddg[cidx] = ddgnm;
            }
        }
    }

  /* update epoch */
  w->epochs[epoch_idx] = t;

  return s;
}

int
msynth_print_spectrum(const char *filename, const double t,
                      const msynth_workspace *w)
{
  size_t n;
  FILE *fp;
  
  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "msynth_print_spectrum: cannot open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  for (n = w->eval_nmin; n <= w->eval_nmax; ++n)
    {
      double gn = msynth_spectrum(t, n, w);
      double dgn = msynth_spectrum_sv(t, n, w);
      double ddgn = msynth_spectrum_sa(t, n, w);
      fprintf(fp, "%zu %.12e %.12e %.12e\n", n, gn, dgn, ddgn);
    }

  fclose(fp);

  return 0;
} /* msynth_print_spectrum() */

int
msynth_print_spectrum_m(const char *filename, const msynth_workspace *w)
{
  size_t n;
  FILE *fp;
  
  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "msynth_print_spectrum_m: cannot open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  n = 1;
  fprintf(fp, "# Field %zu: m/n\n", n++);
  fprintf(fp, "# Field %zu: n + 0.25*(m/n)\n", n++);
  fprintf(fp, "# Field %zu: gnm^2 (nT^2)\n", n++);
  fprintf(fp, "# Field %zu: dgnm^2 (nT/year)^2\n", n++);
  fprintf(fp, "# Field %zu: ddgnm^2 (nT/year^2)^2\n", n++);

  for (n = 1; n <= w->eval_nmax; ++n)
    {
      int ni = (int) n;
      int m;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = msynth_nmidx(n, m, w);
          double gnm = w->c[cidx];
          double dgnm = w->c[cidx + w->sv_offset];
          double ddgnm = w->c[cidx + w->sa_offset];
          double ratio = (double)m / (double)ni;

          fprintf(fp, "%f %f %.12e %.12e %.12e\n",
                  ratio,
                  n + 0.25 * ratio,
                  gnm * gnm,
                  dgnm * dgnm,
                  ddgnm * ddgnm);
        }
    }

  fclose(fp);

  return 0;
} /* msynth_print_spectrum_m() */

/*
msynth_calc_sv()
  Recompute SV coefficients as forward differences between
successive snapshot models. The latest snapshot model is
not modified. Also, set SA coefficients to 0.
*/

int
msynth_calc_sv(msynth_workspace *w)
{
  int s = 0;
  size_t i;

  for (i = 0; i < w->n_snapshot - 1; ++i)
    {
      double *cptr = w->c + i * w->p;
      double *cptr_next = w->c + (i + 1) * w->p;
      double dt = w->epochs[i + 1] - w->epochs[i];
      size_t n;

      for (n = 1; n <= w->nmax; ++n)
        {
          int m, ni = (int) n;

          for (m = -ni; m <= ni; ++m)
            {
              size_t cidx = msynth_nmidx(n, m, w);
              double gnm = cptr[cidx];
              double gnm_next = cptr_next[cidx];
              double dgnm = (gnm_next - gnm) / dt;

              /* store new SV value */
              cptr[cidx + w->sv_offset] = dgnm;

              /* set SA value to 0 */
              cptr[cidx + w->sa_offset] = 0.0;
            }
        }
    }

  return s;
} /* msynth_calc_sv() */

/*******************************************************
 *             INTERNAL ROUTINES                       *
 *******************************************************/

/*
msynth_eval_g()
  Evaluate magnetic field model for given g,dg,ddg coefficients

Inputs: t     - timestamp (decimal years)
        r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        epoch - epoch of closest snapshot model
        g     - main field coefficients (nT)
        dg    - secular variation coefficients (nT/year)
        ddg   - secular acceleration coefficients (nT/year^2)
        B     - (output) magnetic field
                B[0] = B_x
                B[1] = B_y
                B[2] = B_z
                B[3] = |B|
        w     - workspace
*/

static int
msynth_eval_g(const double t, const double r, const double theta,
              const double phi, const double epoch, const double *g,
              const double *dg, const double *ddg, double B[4],
              msynth_workspace *w)
{
  int s = 0;
  size_t n;
  int m;

  /* subtract epoch */
  const double t0 = epoch;
  const double t1 = t - t0;        /* SV term (years) */
  const double t2 = 0.5 * t1 * t1; /* SA term (years^2) */

  /* after this date, use a linear model for gnm */
  const double tend = 2013.5;
  const double t3 = tend - t0;
  const double t4 = 0.5 * t3 * t3;

  s += msynth_green(r, theta, phi, w);

  B[0] = B[1] = B[2] = 0.0;

  for (n = w->eval_nmin; n <= w->eval_nmax; ++n)
    {
      int ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = msynth_nmidx(n, m, w);
          double gnm;

          /*
           * use the parabola taylor series for the period when
           * data was available, and after that use the tangent
           * line to the parabola at the point tend
           */
          if (t <= tend)
            {
              gnm = g[cidx] + dg[cidx] * t1 + ddg[cidx] * t2;
            }
          else
            {
              /* set gnm to the tangent line to the parabola at tend */
              gnm = g[cidx] + dg[cidx] * t3 + ddg[cidx] * t4 +
                    (dg[cidx] + ddg[cidx] * t3) * (t - tend);
            }

          B[0] += gnm * w->dX[cidx];
          B[1] += gnm * w->dY[cidx];
          B[2] += gnm * w->dZ[cidx];
        }
    }

  B[3] = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

  return s;
} /* msynth_eval_g() */

/*
msynth_eval_dgdt()
  Evaluate magnetic field model for given g,dg,ddg coefficients

Inputs: t     - timestamp (decimal years)
        r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        epoch - epoch of closest snapshot model (decimal years)
        g     - main field coefficients (nT)
        dg    - secular variation coefficients (nT/year)
        ddg   - secular acceleration coefficients (nT/year^2)
        dBdt  - (output) magnetic field SV
                dBdt[0] = d/dt B_x
                dBdt[1] = d/dt B_y
                dBdt[2] = d/dt B_z
                dBdt[3] = |dBdt|
        w     - workspace
*/

static int
msynth_eval_dgdt(const double t, const double r, const double theta,
                 const double phi, const double epoch, const double *dg,
                 const double *ddg, double dBdt[4], msynth_workspace *w)
{
  int s = 0;
  size_t n;
  int m;

  /* subtract epoch */
  const double t0 = epoch;
  const double t1 = t - t0;        /* SV term (years) */

  s += msynth_green(r, theta, phi, w);

  dBdt[0] = dBdt[1] = dBdt[2] = 0.0;

  for (n = w->eval_nmin; n <= w->eval_nmax; ++n)
    {
      int ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = msynth_nmidx(n, m, w);
          double gnm = dg[cidx] + ddg[cidx] * t1;

          dBdt[0] += gnm * w->dX[cidx];
          dBdt[1] += gnm * w->dY[cidx];
          dBdt[2] += gnm * w->dZ[cidx];
        }
    }

  dBdt[3] = sqrt(dBdt[0]*dBdt[0] + dBdt[1]*dBdt[1] + dBdt[2]*dBdt[2]);

  return s;
} /* msynth_eval_dgdt() */

/*
msynth_read()
  Read ASCII coefficient file and return a workspace pointer
*/

msynth_workspace *
msynth_read(const char *filename)
{
  size_t nmax;
  double epoch;
  int got_epoch = 0;
  int got_nmax = 0;
  msynth_workspace *w = NULL;
  FILE *fp;
  char buf[2048];

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "msynth_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  while (fgets(buf, 2048, fp) != NULL)
    {
      int c;
      size_t n, cidx;
      int m;
      double gnm, dgnm, ddgnm;

      if (!got_nmax)
        {
          c = sscanf(buf, "%% nmax: %zu", &nmax);
          if (c == 1)
            got_nmax = 1;
        }

      if (!got_epoch)
        {
          c = sscanf(buf, "%% epoch: %lf", &epoch);
          if (c == 1)
            got_epoch = 1;
        }

      if (*buf == '%' || *buf == '#')
        continue;

      if (w == NULL)
        {
          assert(got_nmax == 1 && got_epoch == 1);
          w = msynth_alloc(nmax, 1, &epoch);
        }

      c = sscanf(buf, "%zu %d %lf %lf %lf",
                 &n,
                 &m,
                 &gnm,
                 &dgnm,
                 &ddgnm);
      if (c < 5)
        continue;

      if (n > nmax)
        {
          fprintf(stderr, "msynth_read: error: n = %zu\n", n);
          return w;
        }
      else if (abs(m) > (int) n)
        {
          fprintf(stderr, "msynth_read: error: m(%d) > n(%zu)\n", m, n);
          return w;
        }

      cidx = msynth_nmidx(n, m, w);
      w->c[cidx] = gnm;
      w->c[cidx + w->sv_offset] = dgnm;
      w->c[cidx + w->sa_offset] = ddgnm;
    }

  fclose(fp);

  return w;
} /* msynth_read() */

/*
msynth_read2()
  Read ASCII coefficient file and store in given workspace pointer

Inputs: filename - coefficient file
        w        - (input/output) if NULL, a new workspace is allocated
                   for the new coefficients. If non-NULL, coefficients
                   for new epoch are stored in existing workspace
*/

msynth_workspace *
msynth_read2(const char *filename, msynth_workspace *w)
{
  size_t nmax;
  double epoch;
  int got_epoch = 0;
  int got_nmax = 0;
  FILE *fp;
  char buf[2048];
  double *cptr; /* pointer to coefficients for this epoch */

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "msynth_read2: unable to open %s: %s\n",
              filename, strerror(errno));
      return w;
    }

  while (fgets(buf, 2048, fp) != NULL)
    {
      int c;
      size_t n, cidx;
      int m;
      double gnm, dgnm, ddgnm;

      if (!got_nmax)
        {
          c = sscanf(buf, "%% nmax: %zu", &nmax);
          if (c == 1)
            got_nmax = 1;
        }

      if (!got_epoch)
        {
          c = sscanf(buf, "%% epoch: %lf", &epoch);
          if (c == 1)
            {
              got_epoch = 1;

              /* store new epoch */
              if (w != NULL)
                w->epochs[w->n_epochs++] = epoch;
            }
        }

      if (*buf == '%' || *buf == '#')
        continue;

      assert(got_nmax == 1 && got_epoch == 1);

      if (w == NULL)
        {
          /* allocate new workspace */
          w = msynth_alloc(nmax, 1, &epoch);
        }

      assert(w->n_epochs >= 1);
      cptr = w->c + (w->n_epochs - 1) * w->p;

      c = sscanf(buf, "%zu %d %lf %lf %lf",
                 &n,
                 &m,
                 &gnm,
                 &dgnm,
                 &ddgnm);
      if (c < 5)
        continue;

      if (n > nmax)
        {
          fprintf(stderr, "msynth_read2: error: n = %zu\n", n);
          return w;
        }
      else if (abs(m) > (int) n)
        {
          fprintf(stderr, "msynth_read2: error: m(%d) > n(%zu)\n", m, n);
          return w;
        }

      cidx = msynth_nmidx(n, m, w);
      cptr[cidx] = gnm;
      cptr[cidx + w->sv_offset] = dgnm;
      cptr[cidx + w->sa_offset] = ddgnm;
    }

  fclose(fp);

  return w;
} /* msynth_read2() */

/*
msynth_write()
  Write ascii coefficient file

Inputs: filename - output file
        t        - time (decimal years); model coefficients are
                   extrapolated to time t from nearest model epoch
        w        - workspace
*/

int
msynth_write(const char *filename, const double t,
             const msynth_workspace *w)
{
  int s = 0;
  FILE *fp;
  size_t n;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "msynth_write: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  /* print header information */
  fprintf(fp, "%% Magnetic field model coefficients\n");
  fprintf(fp, "%% nmax:  %zu\n", w->nmax);
  fprintf(fp, "%% epoch: %.4f\n", t);
  fprintf(fp, "%% radius: %.1f\n", w->R);
  fprintf(fp, "%% %3s %5s %20s %20s %20s\n",
          "n",
          "m",
          "MF gnm (nT)",
          "SV gnm (nT/year)",
          "SA gnm (nT/year^2)");

  for (n = 1; n <= w->nmax; ++n)
    {
      int m, ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          double gnm = msynth_get_gnm(t, n, m, w); 
          double dgnm = msynth_get_dgnm(t, n, m, w); 
          double ddgnm = msynth_get_ddgnm(t, n, m, w); 

          fprintf(fp, "%5zu %5d %20.4f %20.4f %20.4f\n",
                  n,
                  m,
                  gnm,
                  dgnm,
                  ddgnm);
        }
    }

  fclose(fp);

  return s;
} /* msynth_write() */

double
msynth_spectrum(const double t, const size_t n, const msynth_workspace *w)
{
  const size_t epoch_idx = msynth_epoch_idx(t, w);
  const double *g = w->c + epoch_idx * w->p;
  int m, ni = (int) n;
  double sum = 0.0;

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = msynth_nmidx(n, m, w);
      double gnm = g[cidx];

      sum += gnm * gnm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
} /* msynth_spectrum() */

double
msynth_spectrum_sv(const double t, const size_t n, const msynth_workspace *w)
{
  const size_t epoch_idx = msynth_epoch_idx(t, w);
  const double *g = w->c + epoch_idx * w->p;
  int m, ni = (int) n;
  double sum = 0.0;

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = msynth_nmidx(n, m, w);
      double dgnm = g[cidx + w->sv_offset];

      sum += dgnm * dgnm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
} /* msynth_spectrum_sv() */

double
msynth_spectrum_sa(const double t, const size_t n, const msynth_workspace *w)
{
  const size_t epoch_idx = msynth_epoch_idx(t, w);
  const double *g = w->c + epoch_idx * w->p;
  int m, ni = (int) n;
  double sum = 0.0;

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = msynth_nmidx(n, m, w);
      double ddgnm = g[cidx + w->sa_offset];

      sum += ddgnm * ddgnm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
} /* msynth_spectrum_sa() */

/*
msynth_nmidx()
  This function returns a unique index in [0,w->p-1] corresponding
to a given (l,m) pair. The array will look like:

[(1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) ...]

(the (0,0) coefficient is not solved for)

Inputs: l - SH degree (> 0)
        m - SH order (-l <= m <= l)

Return: index in [0,nnm-1]
*/

size_t
msynth_nmidx(const size_t n, const int m, const msynth_workspace *w)
{
  size_t base = n * n; /* index of block for this n */
  int offset = m + n;  /* offset within block for this m */
  size_t nmidx;

  if (n == 0)
    {
      fprintf(stderr, "msynth_nmidx: error: n = 0\n");
      return 0;
    }
  else if (n > w->nmax)
    {
      fprintf(stderr, "msynth_nmidx: error: n(%zu) larger than nmax(%zu)\n",
              n, w->nmax);
      return 0;
    }

  nmidx = base + offset;

  /* subtract 1 to exclude (0,0) coefficient */
  return nmidx - 1;
} /* msynth_nmidx() */

/*
msynth_get_epoch()
  Return model epoch corresponding to given time

Inputs: t   - time (decimal years)
        w   - workspace
*/

double
msynth_get_epoch(const double t, const msynth_workspace *w)
{
  const size_t epoch_idx = msynth_epoch_idx(t, w);
  return w->epochs[epoch_idx];
}

/*
msynth_get_mf()
  Return main field coefficient for a given epoch and index

Inputs: t   - epoch time (decimal years)
        idx - coefficient index (n,m) from msynth_nmidx()
        w   - workspace
*/

double
msynth_get_mf(const double t, const size_t idx,
              const msynth_workspace *w)
{
  const size_t epoch_idx = msynth_epoch_idx(t, w);
  const double *g = w->c + epoch_idx * w->p;

  return g[idx];
}

/*
msynth_get_sv()
  Return secular variation coefficient for a given epoch and index

Inputs: t   - epoch time (decimal years)
        idx - coefficient index (n,m) from msynth_nmidx()
        w   - workspace
*/

double
msynth_get_sv(const double t, const size_t idx,
              const msynth_workspace *w)
{
  const size_t epoch_idx = msynth_epoch_idx(t, w);
  const double *g = w->c + epoch_idx * w->p;

  return g[idx + w->sv_offset];
}

/*
msynth_get_sa()
  Return main field coefficient for a given epoch and index

Inputs: t   - epoch time (decimal years)
        idx - coefficient index (n,m) from msynth_nmidx()
        w   - workspace
*/

double
msynth_get_sa(const double t, const size_t idx,
              const msynth_workspace *w)
{
  const size_t epoch_idx = msynth_epoch_idx(t, w);
  const double *g = w->c + epoch_idx * w->p;

  return g[idx + w->sa_offset];
}

/*
msynth_get_gnm()
  Return main field coefficient for a given time and (n,m)

Inputs: t   - time (decimal years)
        n   - spherical harmonic degree
        m   - spherical harmonic order
        w   - workspace
*/

double
msynth_get_gnm(const double t, const size_t n, int m,
               const msynth_workspace *w)
{
  const size_t epoch_idx = msynth_epoch_idx(t, w);
  const double epoch = w->epochs[epoch_idx];
  const size_t cidx = msynth_nmidx(n, m, w);
  const double *g = w->c + epoch_idx * w->p;
  const double *dg = g + w->sv_offset;
  const double *ddg = g + w->sa_offset;
  const double t1 = t - epoch;
  const double t2 = 0.5 * t1 * t1;
  double gnm;

  gnm = g[cidx] + dg[cidx] * t1 + ddg[cidx] * t2;

  return gnm;
} /* msynth_get_gnm() */

/*
msynth_get_dgnm()
  Return SV coefficient for a given time and (n,m)

Inputs: t   - time (decimal years)
        n   - spherical harmonic degree
        m   - spherical harmonic order
        w   - workspace
*/

double
msynth_get_dgnm(const double t, const size_t n, int m,
                const msynth_workspace *w)
{
  const size_t epoch_idx = msynth_epoch_idx(t, w);
  const double epoch = w->epochs[epoch_idx];
  const size_t cidx = msynth_nmidx(n, m, w);
  const double *g = w->c + epoch_idx * w->p;
  const double *dg = g + w->sv_offset;
  const double *ddg = g + w->sa_offset;
  const double t1 = t - epoch;
  double dgnm;

  dgnm = dg[cidx] + ddg[cidx] * t1;

  return dgnm;
} /* msynth_get_dgnm() */

/*
msynth_get_ddgnm()
  Return SA coefficient for a given time and (n,m)

Inputs: t   - time (decimal years)
        n   - spherical harmonic degree
        m   - spherical harmonic order
        w   - workspace
*/

double
msynth_get_ddgnm(const double t, const size_t n, int m,
                 const msynth_workspace *w)
{
  const size_t epoch_idx = msynth_epoch_idx(t, w);
  const size_t cidx = msynth_nmidx(n, m, w);
  const double *g = w->c + epoch_idx * w->p;
  const double *ddg = g + w->sa_offset;

  return ddg[cidx];
} /* msynth_get_ddgnm() */

size_t
msynth_epoch_idx(const double t, const msynth_workspace *w)
{
  if (w->n_snapshot == 1)
    return 0; /* only 1 snapshot model available */

  if (t >= w->epochs[w->n_snapshot - 1])
    return w->n_snapshot - 1;
  else if (t <= w->epochs[0])
    return 0;
  else
    {
      double *ptr = bsearch(&t, w->epochs, w->n_snapshot - 1, sizeof(double), msynth_compare);

      if (ptr == NULL)
        {
          fprintf(stderr, "msynth_epoch_idx: epoch not found for time %f\n", t);
          return 0;
        }

      return ptr - w->epochs;
    }
}

static int
msynth_compare(const void *a, const void *b)
{
  double t = *(double *) a;
  double *db = (double *) b;
  double epoch = *db;
  double epoch_next = *(db + 1);

  if (t < epoch)
    return -1;
  else if (t >= epoch_next)
    return 1;
  else
    return 0;
}
