/*
 * tiegcm.h
 */

#ifndef INCLUDED_tiegcm_h
#define INCLUDED_tiegcm_h

#include <time.h>
#include <netcdf.h>

#define TIEGCM_BIDX(it,ilat,ilon,data)    (CIDX3((it),(data->nt),(ilat),(data->nlat),(ilon),(data->nlon)))

typedef struct
{
  time_t *t;    /* timestamps, size nt */
  double *year; /* year, size nt */
  double *doy;  /* day of year, size nt */
  double *ut;   /* UT (hours), size nt */
  double *glon; /* geodetic longitude (degrees), size nlon */
  double *glat; /* geodetic latitude (degrees), size nlat */

  double *Bx;   /* B_x grid, nt-by-nlon-by-nlat */
  double *By;   /* B_y grid, nt-by-nlon-by-nlat */
  double *Bz;   /* B_z grid, nt-by-nlon-by-nlat */

  double *workx; /* temporary B_x grid, nt_max-by-nlon-by-nlat */
  double *worky; /* temporary B_y grid, nt_max-by-nlon-by-nlat */
  double *workz; /* temporary B_z grid, nt_max-by-nlon-by-nlat */
  size_t nt_max; /* maximum time snapshots in 1 TIEGCM file */

  size_t nt;    /* number of time snapshots stored */
  size_t nlon;  /* length of 'glon' */
  size_t nlat;  /* length of 'glat' */
  size_t ntot;  /* total data allocated */
} tiegcm_data;

/* tiegcm_alloc.c */
tiegcm_data *tiegcm_alloc(const size_t nt, const size_t nlon, const size_t nlat);
tiegcm_data *tiegcm_realloc(const size_t nt, const size_t nlon, const size_t nlat,
                            tiegcm_data *data);
void tiegcm_free(tiegcm_data *data);

/* tiegcm_read.c */
tiegcm_data *tiegcm_read(const char *filename, tiegcm_data *data);

#endif /* INCLUDED_tiegcm_h */
