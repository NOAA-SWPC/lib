/*
 * tiegcm.h
 */

#ifndef INCLUDED_tiegcm_h
#define INCLUDED_tiegcm_h

#include <time.h>
#include <netcdf.h>

typedef struct
{
  time_t *t;    /* timestamps */
  double *year; /* year */
  double *doy;  /* day of year */
  double *ut;   /* UT (hours) */
  size_t n;     /* number of data points stored */
  size_t ntot;  /* total data allocated */
} tiegcm_data;

/* tiegcm_alloc.c */
tiegcm_data *tiegcm_alloc(const size_t n);
tiegcm_data *tiegcm_realloc(const size_t n, tiegcm_data *data);
void tiegcm_free(tiegcm_data *data);

/* tiegcm_read.c */
tiegcm_data *tiegcm_read(const char *filename, tiegcm_data *data);

#endif /* INCLUDED_tiegcm_h */
