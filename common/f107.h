/*
 * f107.h
 * Patrick Alken
 */

#ifndef INCLUDED_f107_h
#define INCLUDED_f107_h

#include <time.h>

#define F107_MAX_BUFFER    512
#define F107_MAX_YEAR      500

/* first year to store data */
#define F107_FIRST_YEAR    1900

typedef struct
{
  /* data for a whole year */
  double data[12][31];
} f107_data;

typedef struct
{
  f107_data f107[F107_MAX_YEAR];
} f107_workspace;

#define F107_DATAFILE      "/home/palken/data/F107/F107.txt"

/*
 * Prototypes
 */

f107_workspace *f107_alloc(const char *filename);
void f107_free(f107_workspace *w);
int f107_get(time_t t, double *result, f107_workspace *w);
int f107a_get(time_t t, double *result, f107_workspace *w);
int f107a_get2(time_t t, double *result, f107_workspace *w);

#endif /* INCLUDED_f107_h */
