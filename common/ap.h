/*
 * ap.h
 * Patrick Alken
 */

#ifndef INCLUDED_ap_h
#define INCLUDED_ap_h

#include <time.h>

typedef struct
{
  double *ap_data;
  time_t t0;         /* timestamp of first date */
  size_t n;          /* number of ap data points */
} ap_workspace;

#define AP_MAX_DATA      50000
#define AP_MAX_BUFFER    512

#define AP_DATA_FILE     "/home/palken/data/SPIDR/AP/AP.txt"

#define AP_MIN(a,b)      ((a) < (b) ? (a) : (b))

/*
 * Prototypes
 */

ap_workspace *ap_alloc(const char *filename);
void ap_free(ap_workspace *w);
int ap_get(time_t t, double *result, ap_workspace *w);

extern int putenv(char *string);
double round(double x);

#endif /* INCLUDED_ap_h */
