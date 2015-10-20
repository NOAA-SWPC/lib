/*
 * swarmeef.h
 */

#ifndef INCLUDED_swarmeef_h
#define INCLUDED_swarmeef_h

#include <time.h>
#include <cdf.h>

#define SWARMEEF_MAX_BUF          8192

typedef struct
{
  double *t;         /* timestamp (CDF_EPOCH) */

  double *latitude;  /* latitude in deg */
  double *longitude; /* longitude in deg */
  double *EEF;       /* EEF (V/m) */
  double *RelErr;    /* relative error */

  unsigned short *Flags; /* data flags */

  size_t n;          /* number of data points stored */
  size_t ntot;       /* total data allocated */

  char cdf_error[SATDATA_MAX_BUF];
} swarm_eef;

/*
 * Prototypes
 */

int putenv(char *string);

swarm_eef *swarm_eef_alloc(size_t n);
swarm_eef *swarm_eef_realloc(size_t n, swarm_eef *data);
void swarm_eef_free(swarm_eef *data);
char *cdf_error(CDFstatus status);

swarm_eef *swarm_eef_read(const char *filename, swarm_eef *data);
swarm_eef *swarm_eef_read_idx(const char *filename, const int verbose);

#endif /* INCLUDED_swarmeef_h */
