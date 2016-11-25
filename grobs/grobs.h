/*
 * grobs.h
 *
 * Data structures for ground observatory data
 */

#ifndef INCLUDED_grobs_h
#define INCLUDED_grobs_h

#include <time.h>

#define GROBS_MAX_BUFFER  4096

/* maximum years of data to allocate at once */
#define GROBS_MAX_YEAR       5

/* 1-minute data */
typedef struct
{
  time_t *t;     /* timestamp */
  double *X;     /* X component (nT) */
  double *Y;     /* Y component (nT) */
  double *Z;     /* Z component (nT) */
  double *H;     /* H component (nT) */
  double *D;     /* declination (deg) */
  double *I;     /* inclination (deg) */

  double glat;   /* geodetic latitude (deg) */
  double glon;   /* geodetic longitude (deg) */
  char name[16]; /* station name */

  size_t ntot;   /* total data allocated */
  size_t n;      /* number of data stored */
} grobs_data;

/*
 * Prototypes
 */

grobs_data *grobs_alloc(const size_t ntot);
grobs_data *grobs_realloc(const size_t ntot, grobs_data *data);
void grobs_free(grobs_data *data);

/* iaga.c */
grobs_data *grobs_iaga_read(const char *filename, grobs_data *data);

/* wamnet.c */
grobs_data *grobs_wamnet_read(const char *filename, grobs_data *data);

#endif /* INCLUDED_grobs_h */
