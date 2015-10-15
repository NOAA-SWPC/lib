/*
 * dst.h
 * Patrick Alken
 */

#ifndef INCLUDED_dst_h
#define INCLUDED_dst_h

#include <time.h>

#define DST_MAX_BUFFER    2048
#define DST_MAX_YEAR      500
#define DST_FIRST_YEAR    1900

/* maximum number of Dst values per day */
#define DST_MAX_DATA      100

typedef struct
{
  double *dst;
  double *est;
  double *ist;
  time_t *t;
  size_t n; /* number of dst data on this day */
} dst_data;

typedef struct
{
  dst_data data[DST_MAX_YEAR][12][31];
  time_t t0;         /* timestamp of first date */
  size_t n;          /* total number of dst data points */
} dst_workspace;


#define DST_DATA_FILE     "/nfs/satmag/data/Indices/Dst/Est_Ist_index.pli"

#define DST_MIN(a,b)      ((a) < (b) ? (a) : (b))

/*
 * Prototypes
 */

dst_workspace *dst_alloc(const char *filename);
void dst_free(dst_workspace *w);
int dst_get(time_t t, double *result, dst_workspace *w);
int est_get(time_t t, double *result, dst_workspace *w);
int ist_get(time_t t, double *result, dst_workspace *w);

extern int putenv(char *string);
double round(double x);

#endif /* INCLUDED_dst_h */
