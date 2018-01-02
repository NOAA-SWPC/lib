/*
 * tiegcm3d.h
 */

#ifndef INCLUDED_tiegcm3d_h
#define INCLUDED_tiegcm3d_h

#include <time.h>
#include <netcdf.h>

#define TIEGCM3D_IDX(it,ir,ilat,ilon,data)    (CIDX4((it),(data->nt),(ir),(data->nr),(ilon),(data->nlon),(ilat),(data->nlat)))

typedef struct
{
  time_t *t;    /* timestamps, size nt */
  double *year; /* year, size nt */
  double *doy;  /* day of year, size nt */
  double *ut;   /* UT (hours), size nt */
  double *r;    /* radius (km), size nr */
  double *glon; /* geodetic longitude (degrees), size nlon */
  double *glat; /* geodetic latitude (degrees), size nlat */

  double *Jr;   /* J_r grid, nt-by-nr-by-nlon-by-nlat */
  double *Jt;   /* J_t grid, nt-by-nr-by-nlon-by-nlat */
  double *Jp;   /* J_p grid, nt-by-nr-by-nlon-by-nlat */

  double *workx; /* temporary B_x grid, nt_max-by-nlon-by-nlat */
  double *worky; /* temporary B_y grid, nt_max-by-nlon-by-nlat */
  double *workz; /* temporary B_z grid, nt_max-by-nlon-by-nlat */
  size_t nt_max; /* maximum time snapshots in 1 TIEGCM file */

  size_t nt;    /* number of time snapshots stored */
  size_t nr;    /* length of 'r' */
  size_t nlon;  /* length of 'glon' */
  size_t nlat;  /* length of 'glat' */
  size_t ntot;  /* total data allocated */
} tiegcm3d_data;

/* tiegcm3d_alloc.c */
tiegcm3d_data *tiegcm3d_alloc(const size_t nt, const size_t nr, const size_t nlon, const size_t nlat);
tiegcm3d_data *tiegcm3d_realloc(const size_t nt, const size_t nr, const size_t nlon, const size_t nlat,
                                tiegcm3d_data *data);
void tiegcm3d_free(tiegcm3d_data *data);

/* tiegcm3d_read.c */
tiegcm3d_data *tiegcm3d_read(const char *filename, tiegcm3d_data *data);

#endif /* INCLUDED_tiegcm3d_h */
