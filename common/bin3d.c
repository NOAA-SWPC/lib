/*
 * bin3d.c
 *
 * Bin a dataset in 3 dimensions
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <gsl/gsl_rstat.h>

#include "bin3d.h"

static size_t bin3d_findbin(const double val, const double a, const double b, const size_t n);
static double bin3d_interp1d(const double a, const double b,
                             const double fa, const double fb, const double x);

bin3d_workspace *
bin3d_alloc(const double xmin, const double xmax, const size_t nx,
            const double ymin, const double ymax, const size_t ny,
            const double zmin, const double zmax, const size_t nz)
{
  bin3d_workspace *w;
  size_t i;

  w = calloc(1, sizeof(bin3d_workspace));
  if (!w)
    return 0;

  w->nbins = nx * ny * nz;
  
  w->bins = malloc(w->nbins * sizeof(gsl_rstat_workspace *));
  if (!w->bins)
    {
      fprintf(stderr, "bin3d_alloc: malloc failed: %s\n", strerror(errno));
      bin3d_free(w);
      return 0;
    }

  for (i = 0; i < w->nbins; ++i)
    {
      w->bins[i] = gsl_rstat_alloc();
      if (!w->bins[i])
        {
          fprintf(stderr, "bin3d_alloc: gsl_rstat_alloc failed\n");
          bin3d_free(w);
          return 0;
        }
    }

  w->nx = nx;
  w->ny = ny;
  w->nz = nz;

  w->xmin = xmin;
  w->xmax = xmax;
  w->ymin = ymin;
  w->ymax = ymax;
  w->zmin = zmin;
  w->zmax = zmax;

  return w;
} /* bin3d_alloc() */

void
bin3d_free(bin3d_workspace *w)
{
  size_t i;

  if (w->bins)
    {
      for (i = 0; i < w->nbins; ++i)
        {
          if (w->bins[i])
            gsl_rstat_free(w->bins[i]);
        }

      free(w->bins);
    }

  free(w);
} /* bin3d_free() */

int
bin3d_add_element(const double x, const double y, const double z,
                  const double data, bin3d_workspace *w)
{
  int s = 0;
  size_t binx, biny, binz;

  if (x < w->xmin || x > w->xmax)
    {
      fprintf(stderr, "bin_add_element: error: x outside allowed range: %f\n", x);
      return -1;
    }
  else if (y < w->ymin || y > w->ymax)
    {
      fprintf(stderr, "bin_add_element: error: y outside allowed range: %f\n", y);
      return -1;
    }
  else if (z < w->zmin || z > w->zmax)
    {
      fprintf(stderr, "bin_add_element: error: z outside allowed range: %f\n", z);
      return -1;
    }

  binx = bin3d_findbin(x, w->xmin, w->xmax, w->nx);
  biny = bin3d_findbin(y, w->ymin, w->ymax, w->ny);
  binz = bin3d_findbin(z, w->zmin, w->zmax, w->nz);

  gsl_rstat_add(data, w->bins[BIN3D_IDX(binx, biny, binz, w)]);

  return s;
} /* bin3d_add_element() */

double
bin3d_mean(const double x, const double y, const double z, const bin3d_workspace *w)
{
  size_t binx = bin3d_findbin(x, w->xmin, w->xmax, w->nx);
  size_t biny = bin3d_findbin(y, w->ymin, w->ymax, w->ny);
  size_t binz = bin3d_findbin(z, w->zmin, w->zmax, w->nz);
  double mean;

  mean = gsl_rstat_mean(w->bins[BIN3D_IDX(binx, biny, binz, w)]);

  return mean;
} /* bin3d_mean() */

double
bin3d_sd(const double x, const double y, const double z, const bin3d_workspace *w)
{
  size_t binx = bin3d_findbin(x, w->xmin, w->xmax, w->nx);
  size_t biny = bin3d_findbin(y, w->ymin, w->ymax, w->ny);
  size_t binz = bin3d_findbin(z, w->zmin, w->zmax, w->nz);
  double sd;

  sd = gsl_rstat_sd(w->bins[BIN3D_IDX(binx, biny, binz, w)]);

  return sd;
} /* bin3d_sd() */

double
bin3d_median(const double x, const double y, const double z, const bin3d_workspace *w)
{
  size_t binx = bin3d_findbin(x, w->xmin, w->xmax, w->nx);
  size_t biny = bin3d_findbin(y, w->ymin, w->ymax, w->ny);
  size_t binz = bin3d_findbin(z, w->zmin, w->zmax, w->nz);
  double median;

  median = gsl_rstat_median(w->bins[BIN3D_IDX(binx, biny, binz, w)]);

  return median;
} /* bin3d_median() */

size_t
bin3d_n(const double x, const double y, const double z, const bin3d_workspace *w)
{
  size_t binx = bin3d_findbin(x, w->xmin, w->xmax, w->nx);
  size_t biny = bin3d_findbin(y, w->ymin, w->ymax, w->ny);
  size_t binz = bin3d_findbin(z, w->zmin, w->zmax, w->nz);
  size_t n;

  n = gsl_rstat_n(w->bins[BIN3D_IDX(binx, biny, binz, w)]);

  return n;
} /* bin3d_n() */

/*
bin3d_xyzval()
  Determine x, y and z bin values for a given bin

Inputs: i - x bin (i \in [0,nx-1])
        j - y bin (j \in [0,ny-1])
        k - z bin (k \in [0,nz-1])
        x - (output) x = xmin + i * dx
        y - (output) y = ymin + j * dy
        z - (output) z = zmin + k * dz
        w - workspace

Return: success or error
*/

int
bin3d_xyzval(const size_t i, const size_t j, const size_t k,
            double *x, double *y, double *z, const bin3d_workspace *w)
{
  int s = 0;
  const double dx = (w->xmax - w->xmin) / (w->nx - 1.0);
  const double dy = (w->ymax - w->ymin) / (w->ny - 1.0);
  const double dz = (w->zmax - w->zmin) / (w->nz - 1.0);

  *x = w->xmin + i * dx;
  *y = w->ymin + j * dy;
  *z = w->zmin + k * dz;

  return s;
} /* bin3d_xyzval() */

/* print 2d grid at fixed z to file */
int
bin3d_print_z(const char *filename, const double z, const bin3d_workspace *w)
{
  int s = 0;
  size_t i, j;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "bin3d_print_z: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Grid at z = %f\n", z);
  fprintf(fp, "# Field %zu: x\n", i++);
  fprintf(fp, "# Field %zu: y\n", i++);
  fprintf(fp, "# Field %zu: mean\n", i++);
  fprintf(fp, "# Field %zu: median\n", i++);
  fprintf(fp, "# Field %zu: stddev\n", i++);
  fprintf(fp, "# Field %zu: number of data points in bin\n", i++);

  for (i = 0; i < w->nx; ++i)
    {
      for (j = 0; j < w->ny; ++j)
        {
          double x, y, dummy;
          double mean, sd, median;
          size_t n;

          bin3d_xyzval(i, j, 0, &x, &y, &dummy, w);

          mean = bin3d_mean(x, y, z, w);
          sd = bin3d_sd(x, y, z, w);
          median = bin3d_median(x, y, z, w);
          n = bin3d_n(x, y, z, w);

          fprintf(fp, "%f %f %.12e %.12e %.12e %zu\n", x, y, mean, median, sd, n);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

/*
bin3d_findbin()
  Determine bin number corresponding to val \in [a, b] where n is
the number of bins using a linear map separately in x, y, z
directions
*/

static size_t
bin3d_findbin(const double val, const double a, const double b, const size_t n)
{
  size_t binno = (size_t) round(bin3d_interp1d(a, b, 0.0, n - 1.0, val));

  assert(binno < n);

  return binno;
} /* bin3d_findbin() */

/*
bin3d_interp1d()
  1D interpolation: find f(x) using endpoints f(a) and f(b)
*/

static double
bin3d_interp1d(const double a, const double b,
               const double fa, const double fb, const double x)
{
  double fx;

  fx = (b - x) / (b - a) * fa + (x - a) / (b - a) * fb;

  return fx;
} /* bin3d_interp1d() */
