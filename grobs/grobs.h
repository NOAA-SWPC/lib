/*
 * grobs.h
 *
 * Data structures for ground observatory data
 */

#ifndef INCLUDED_grobs_h
#define INCLUDED_grobs_h

#include <time.h>

/* maximum number of years of data to read */
#define GROBS_MAX_YEAR    5

#define GROBS_MAX_BUFFER  4096

/* 1-minute data */
typedef struct
{
  time_t t[GROBS_MAX_YEAR * 527040]; /* timestamp */
  double X[GROBS_MAX_YEAR * 527040]; /* X component (nT) */
  double Y[GROBS_MAX_YEAR * 527040]; /* Y component (nT) */
  double Z[GROBS_MAX_YEAR * 527040]; /* Z component (nT) */
  double D[GROBS_MAX_YEAR * 527040]; /* declination (deg) */
  double I[GROBS_MAX_YEAR * 527040]; /* inclination (deg) */
  size_t n;                          /* number of data stored */
} grobs_data;

/*
 * Prototypes
 */

grobs_data *iaga_read_HDZF(const char *filename, grobs_data *data);
grobs_data *iaga_read_XYZF(const char *filename, grobs_data *data);

#endif /* INCLUDED_grobs_h */
