/*
 * bin3d.h
 */

#ifndef INCLUDED_bin3d_h
#define INCLUDED_bin3d_h

#include <gsl/gsl_rstat.h>

typedef struct
{
  size_t nx;
  size_t ny;
  size_t nz;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;

  gsl_rstat_workspace **bins;
  size_t nbins;
} bin3d_workspace;

#define CIDX3(i,Ni,j,Nj,k,Nk) ((k) + (Nk)*((j) + (Nj)*(i)))
#define BIN3D_IDX(x, y, z, w)    (CIDX3((x), (w)->nx, (y), (w)->ny, (z), (w)->nz))

/*
 * Prototypes
 */

bin3d_workspace *bin3d_alloc(const double xmin, const double xmax, const size_t nx,
                             const double ymin, const double ymax, const size_t ny,
                             const double zmin, const double zmax, const size_t nz);
void bin3d_free(bin3d_workspace *w);
int bin3d_add_element(const double x, const double y, const double z,
                      const double data, bin3d_workspace *w);
double bin3d_mean(const double x, const double y, const double z, const bin3d_workspace *w);
double bin3d_median(const double x, const double y, const double z, const bin3d_workspace *w);
double bin3d_sd(const double x, const double y, const double z, const bin3d_workspace *w);
size_t bin3d_n(const double x, const double y, const double z, const bin3d_workspace *w);
int bin3d_xyval(const size_t i, const size_t j, const size_t k,
                double *x, double *y, double *z, const bin3d_workspace *w);
int bin3d_print_z(const char *filename, const double z, const bin3d_workspace *w);

#endif /* INCLUDED_bin3d_h */
