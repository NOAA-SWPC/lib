/*
 * mfield_eval.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <gsl/gsl_sf_legendre.h>

#include "mfield_eval.h"

static int mfield_eval_g(const double t, const double r, const double theta,
                         const double phi, const double *g, const double *dg,
                         const double *ddg, double B[4], mfield_eval_workspace *w);
static int mfield_eval_dgdt(const double t, const double r,
                            const double theta, const double phi,
                            const double *dg, const double *ddg,
                            double dBdt[4], mfield_eval_workspace *w);

mfield_eval_workspace *
mfield_eval_alloc(const size_t nmax, const double epoch)
{
  mfield_eval_workspace *w;
  size_t nnm_tot;
  size_t plm_array_size = gsl_sf_legendre_array_n(nmax);

  w = calloc(1, sizeof(mfield_eval_workspace));
  if (!w)
    return 0;

  w->nmax = nmax;
  w->R = 6371.2;
  w->epoch = epoch;

  nnm_tot = (nmax + 1) * (nmax + 1);

  /* exclude (0,0) coefficient */
  w->nnm = nnm_tot - 1;
  w->p = 3 * w->nnm; /* include SV,SA coefficients */

  w->sv_offset = w->nnm;
  w->sa_offset = 2 * w->nnm;

  w->c = malloc(w->p * sizeof(double));

  w->cosmphi = malloc((nmax + 1) * sizeof(double));
  w->sinmphi = malloc((nmax + 1) * sizeof(double));

  w->Plm = malloc(plm_array_size * sizeof(double));
  w->dPlm = malloc(plm_array_size * sizeof(double));
  if (!w->Plm || !w->dPlm)
    {
      mfield_eval_free(w);
      return 0;
    }

  w->dX = malloc(w->nnm * sizeof(double));
  w->dY = malloc(w->nnm * sizeof(double));
  w->dZ = malloc(w->nnm * sizeof(double));

  return w;
} /* mfield_eval_alloc() */

void
mfield_eval_free(mfield_eval_workspace *w)
{
  if (w->c)
    free(w->c);

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

  free(w);
} /* mfield_eval_free() */

/*
mfield_eval()
  Evaluate magnetic field model at given point

Inputs: t     - timestamp (decimal years)
        r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - geocentric longitude (radians)
        B     - (output) magnetic field
                B[0] = B_x
                B[1] = B_y
                B[2] = B_z
                B[3] = |B|
        w     - workspace
*/

int
mfield_eval(const double t, const double r, const double theta,
            const double phi, double B[4], mfield_eval_workspace *w)
{
  int s = 0;
  double *g = &(w->c[0]);
  double *dg = &(w->c[w->sv_offset]);
  double *ddg = &(w->c[w->sa_offset]);

  s = mfield_eval_g(t, r, theta, phi, g, dg, ddg, B, w);

  return s;
} /* mfield_eval() */

/*
mfield_eval_dBdt()
  Evaluate magnetic field model at given point

Inputs: t     - timestamp (decimal years)
        r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - geocentric longitude (radians)
        dBdt  - (output) magnetic field SV
                dBdt[0] = d/dt B_x
                dBdt[1] = d/dt B_y
                dBdt[2] = d/dt B_z
                dBdt[3] = |dBdt|
        w     - workspace
*/

int
mfield_eval_dBdt(const double t, const double r, const double theta,
                 const double phi, double dBdt[4], mfield_eval_workspace *w)
{
  int s = 0;
  double *dg = &(w->c[w->sv_offset]);
  double *ddg = &(w->c[w->sa_offset]);

  s = mfield_eval_dgdt(t, r, theta, phi, dg, ddg, dBdt, w);

  return s;
} /* mfield_eval_dBdt() */

/*******************************************************
 *             INTERNAL ROUTINES                       *
 *******************************************************/

/*
mfield_eval_g()
  Evaluate magnetic field model for given g,dg,ddg coefficients

Inputs: t     - timestamp (decimal years)
        r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
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
mfield_eval_g(const double t, const double r, const double theta,
              const double phi, const double *g, const double *dg,
              const double *ddg, double B[4], mfield_eval_workspace *w)
{
  int s = 0;
  size_t n;
  int m;

  /* subtract epoch */
  const double t0 = w->epoch;
  const double t1 = t - t0;        /* SV term (years) */
  const double t2 = 0.5 * t1 * t1; /* SA term (years^2) */

  /* after this date, use a linear model for gnm */
  const double tend = 2013.5;
  const double t3 = tend - t0;
  const double t4 = 0.5 * t3 * t3;

  s += mfield_eval_green(r, theta, phi, w);

  B[0] = B[1] = B[2] = 0.0;

  for (n = 1; n <= w->nmax; ++n)
    {
      int ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_eval_nmidx(n, m);
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
} /* mfield_eval_g() */

/*
mfield_eval_dgdt()
  Evaluate magnetic field model for given g,dg,ddg coefficients

Inputs: t     - timestamp (decimal years)
        r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
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
mfield_eval_dgdt(const double t, const double r, const double theta,
                 const double phi, const double *dg, const double *ddg,
                 double dBdt[4], mfield_eval_workspace *w)
{
  int s = 0;
  size_t n;
  int m;

  /* subtract epoch */
  const double t0 = w->epoch;
  const double t1 = t - t0;        /* SV term (years) */

  s += mfield_eval_green(r, theta, phi, w);

  dBdt[0] = dBdt[1] = dBdt[2] = 0.0;

  for (n = 1; n <= w->nmax; ++n)
    {
      int ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_eval_nmidx(n, m);
          double gnm = dg[cidx] + ddg[cidx] * t1;

          dBdt[0] += gnm * w->dX[cidx];
          dBdt[1] += gnm * w->dY[cidx];
          dBdt[2] += gnm * w->dZ[cidx];
        }
    }

  dBdt[3] = sqrt(dBdt[0]*dBdt[0] + dBdt[1]*dBdt[1] + dBdt[2]*dBdt[2]);

  return s;
} /* mfield_eval_g() */

/*
mfield_eval_ext()
  Evaluate external magnetic field model at given point

Inputs: t     - timestamp (decimal years)
        r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - geocentric longitude (radians)
        E_st  - external Dst component (nT)
        I_st  - induced Dst component (nT)
        B     - (output) magnetic field
                B[0] = B_x
                B[1] = B_y
                B[2] = B_z
                B[3] = |B|
        w     - workspace
*/

int
mfield_eval_ext(const double t, const double r, const double theta, const double phi,
                const double E_st, const double I_st, double B[4], mfield_eval_workspace *w)
{
  int s = 0;

  /* these are the RC0 MF and SV coefficients from HDGM */
  const double rc_g10 = 23.6956;
  const double rc_g11 = 1.2411;
  const double rc_h11 = -3.8961;
  const double rc_dg10 = 3.6339;
  const double rc_dg11 = 0.1903;
  const double rc_dh11 = -0.5975;

  size_t i;
  double g10, g11, h11, g1;
  double g_ext[3], dg_ext[3];
  double B_ext[4], B_rc[4];

  g10 = w->c[mfield_eval_nmidx(1, 0)];
  g11 = w->c[mfield_eval_nmidx(1, 1)];
  h11 = w->c[mfield_eval_nmidx(1, -1)];
  g1 = sqrt(g10*g10 + g11*g11 + h11*h11);

  /* scale to a unit vector */
  g_ext[mfield_eval_nmidx(1, 0)] = g10 / g1;
  g_ext[mfield_eval_nmidx(1, 1)] = g11 / g1;
  g_ext[mfield_eval_nmidx(1, -1)] = h11 / g1;

  for (i = 0; i < 3; ++i)
    dg_ext[i] = 0.0;

  mfield_eval_g_ext(t, r, theta, phi, E_st, I_st,
                    g_ext, dg_ext, B_ext, w);

  /* now compute the steady ring current contribution from HDGM coefficients */
  g_ext[mfield_eval_nmidx(1, 0)] = rc_g10;
  g_ext[mfield_eval_nmidx(1, 1)] = rc_g11;
  g_ext[mfield_eval_nmidx(1, -1)] = rc_h11;

  dg_ext[mfield_eval_nmidx(1, 0)] = rc_dg10;
  dg_ext[mfield_eval_nmidx(1, 1)] = rc_dg11;
  dg_ext[mfield_eval_nmidx(1, -1)] = rc_dh11;

  mfield_eval_g_ext(t, r, theta, phi, 1.0, 0.0,
                    g_ext, dg_ext, B_rc, w);

  B[0] = B_ext[0] + B_rc[0];
  B[1] = B_ext[1] + B_rc[1];
  B[2] = B_ext[2] + B_rc[2];
  B[3] = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

  return s;
} /* mfield_eval_ext() */

/*
mfield_eval_g_ext()
  Evaluate external magnetic field model at given point

Inputs: t     - timestamp (decimal years)
        r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - geocentric longitude (radians)
        E_st  - external Dst component (nT)
        I_st  - induced Dst component (nT)
        B     - (output) magnetic field
                B[0] = B_x
                B[1] = B_y
                B[2] = B_z
                B[3] = |B|
        w     - workspace
*/

int
mfield_eval_g_ext(const double t, const double r, const double theta, const double phi,
                  const double E_st, const double I_st,
                  const double *g, const double *dg,
                  double B[4], mfield_eval_workspace *w)
{
  int s = 0;
  int m;
  double sint = sin(theta);
  double rterm;
  double dt = t - w->epoch;

  /* no radial term for n = 1 external, only internal */
  rterm = pow(w->R / r, 3.0);

  /* compute associated legendres */
  gsl_sf_legendre_deriv_array(GSL_SF_LEGENDRE_SCHMIDT,
                              w->nmax, cos(theta), w->Plm, w->dPlm);

  B[0] = B[1] = B[2] = 0.0;

  for (m = 0; m <= 1; ++m)
    {
      size_t cidx = mfield_eval_nmidx(1, m);
      size_t pidx = gsl_sf_legendre_array_index(1, m);
      double g1m = g[cidx] + dt * dg[cidx];
      double h1m = 0.0;

      if (m != 0)
        {
          cidx = mfield_eval_nmidx(1, -m);
          h1m = g[cidx] + dt * dg[cidx];
        }

      /* external contribution */
      B[0] += E_st * (g1m * cos(m * phi) + h1m * sin(m * phi)) * w->dPlm[pidx] * (-sint);
      B[1] += E_st * m / sint * (g1m * sin(m * phi) - h1m * cos(m * phi)) * w->Plm[pidx];
      B[2] += E_st * (g1m * cos(m * phi) + h1m * sin(m * phi)) * w->Plm[pidx];

      /* internal contribution */
      B[0] += I_st * rterm * (g1m * cos(m * phi) + h1m * sin(m * phi)) * w->dPlm[pidx] * (-sint);
      B[1] += I_st * rterm * m / sint * (g1m * sin(m * phi) - h1m * cos(m * phi)) * w->Plm[pidx];
      B[2] -= I_st * 2.0 * rterm * (g1m * cos(m * phi) + h1m * sin(m * phi)) * w->Plm[pidx];
    }

  B[3] = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

  return s;
} /* mfield_eval_g_ext() */

/*
mfield_eval_read()
  Read ASCII coefficient file and return a workspace pointer
*/

mfield_eval_workspace *
mfield_eval_read(const char *filename)
{
  size_t nmax;
  double epoch;
  int got_epoch = 0;
  int got_nmax = 0;
  mfield_eval_workspace *w = NULL;
  FILE *fp;
  char buf[2048];

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "mfield_eval_read: unable to open %s: %s\n",
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
          w = mfield_eval_alloc(nmax, epoch);
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
          fprintf(stderr, "mfield_eval_read: error: n = %zu\n", n);
          return w;
        }
      else if (abs(m) > (int) n)
        {
          fprintf(stderr, "mfield_eval_read: error: m(%d) > n(%zu)\n", m, n);
          return w;
        }

      cidx = mfield_eval_nmidx(n, m);
      w->c[cidx] = gnm;
      w->c[cidx + w->sv_offset] = dgnm;
      w->c[cidx + w->sa_offset] = ddgnm;
    }

  fclose(fp);

  return w;
} /* mfield_eval_read() */

/*
mfield_eval_green()
  Compute Green's functions for X,Y,Z spherical harmonic expansion. These
are simply the basis functions multiplying the g_{nm} and h_{nm} coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Plm
w->dPlm
w->sinmphi
w->cosmphi

2) The output Green's functions are stored in w->dX, w->dY, w->dZ
*/

int
mfield_eval_green(const double r, const double theta, const double phi,
                  mfield_eval_workspace *w)
{
  int s = 0;
  size_t n;
  int m;
  const double sint = sin(theta);
  const double cost = cos(theta);
  double ratio = w->R / r;
  double term = ratio * ratio;     /* (a/r)^{n+2} */

  /* precompute cos(m phi) and sin(m phi) */
  for (n = 0; n <= w->nmax; ++n)
    {
      w->cosmphi[n] = cos(n * phi);
      w->sinmphi[n] = sin(n * phi);
    }

  /* compute associated legendres */
  gsl_sf_legendre_deriv_array(GSL_SF_LEGENDRE_SCHMIDT,
                              w->nmax, cost, w->Plm, w->dPlm);

  for (n = 1; n <= w->nmax; ++n)
    {
      int ni = (int) n;

      /* (a/r)^{n+2} */
      term *= ratio;

      for (m = -ni; m <= ni; ++m)
        {
          int mabs = abs(m);
          size_t cidx = mfield_eval_nmidx(n, m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);

          if (m < 0)
            {
              /* h_{nm} */
              w->dX[cidx] = term * w->sinmphi[mabs] * w->dPlm[pidx] * (-sint);
              w->dY[cidx] = -term / sint * mabs * w->cosmphi[mabs] * w->Plm[pidx];
              w->dZ[cidx] = -(n + 1.0) * term * w->sinmphi[mabs] * w->Plm[pidx];
            }
          else
            {
              /* g_{nm} */
              w->dX[cidx] = term * w->cosmphi[mabs] * w->dPlm[pidx] * (-sint);
              w->dY[cidx] = term / sint * mabs * w->sinmphi[mabs] * w->Plm[pidx];
              w->dZ[cidx] = -(n + 1.0) * term * w->cosmphi[mabs] * w->Plm[pidx];
            }
        }
    }

  return s;
} /* mfield_eval_green() */

/*
mfield_eval_nmidx()
  This function returns a unique index in [0,w->p-1] corresponding
to a given (l,m) pair. The array will look like:

[(1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) ...]

(the (0,0) coefficient is not solved for)

Inputs: l - SH degree (> 0)
        m - SH order (-l <= m <= l)

Return: index in [0,nnm-1]
*/

size_t
mfield_eval_nmidx(const size_t n, const int m)
{
  size_t base = n * n; /* index of block for this n */
  int offset = m + n;  /* offset within block for this m */
  size_t nmidx;

  if (n == 0)
    {
      fprintf(stderr, "mfield_eval_nmidx: error: n = 0\n");
      return 0;
    }

  nmidx = base + offset;

  /* subtract 1 to exclude (0,0) coefficient */
  return nmidx - 1;
} /* mfield_eval_nmidx() */
