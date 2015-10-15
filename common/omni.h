/*
 * omni.h
 * Patrick Alken
 */

#ifndef INCLUDED_omni_h
#define INCLUDED_omni_h

#include <time.h>

#define OMNI_MAX_BUFFER    512
#define OMNI_MAX_YEAR      500

/* first year to store data */
#define OMNI_FIRST_YEAR    1900

typedef struct
{
  double Bx; /* IMF Bx (nT) */
  double By; /* IMF By (nT) */
  double Bz; /* IMF Bz (nT) */
  double V;  /* solar wind velocity (km/s) */
} omni_data;

typedef struct
{
  /* data is every hour of each day */
  omni_data data[OMNI_MAX_YEAR][366][24];
} omni_workspace;

#define OMNI_DATA_DIR     "/home/palken/data/OMNI"
#define OMNI_MISSING      (1.0e9)

/*
 * Prototypes
 */

omni_workspace *omni_alloc(const char *datadir);
void omni_free(omni_workspace *w);
int omni_get(time_t t, double B[4], double *V, omni_workspace *w);

#endif /* INCLUDED_omni_h */
