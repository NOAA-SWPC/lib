/*
 * kp.h
 * Patrick Alken
 */

#ifndef INCLUDED_kp_h
#define INCLUDED_kp_h

#include <time.h>

typedef struct
{
  double *kp_data;
  time_t t0;         /* timestamp of first date */
  size_t n;          /* number of kp data points */
} kp_workspace;

#define KP_MAX_DATA      50000
#define KP_MAX_BUFFER    512

#define KP_DATA_FILE     "/home/palken/data/SPIDR/KP/KP.txt"

#define KP_MIN(a,b)      ((a) < (b) ? (a) : (b))

/*
 * Prototypes
 */

kp_workspace *kp_alloc(const char *filename);
void kp_free(kp_workspace *w);
int kp_get(time_t t, double *result, kp_workspace *w);

extern int putenv(char *string);
double round(double x);

#endif /* INCLUDED_kp_h */
