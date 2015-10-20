/*
 * hwm.h
 * Patrick Alken
 */

#ifndef INCLUDED_hwm_h
#define INCLUDED_hwm_h

#include <indices/indices.h>

/* define to use new HWM07 model instead of old HWM93 model */
#define USE_HWM07

typedef struct
{
  double f107_override;   /* override values of F10.7 and F10.7A */
  double f107a_override;
  kp_workspace *kp_workspace_p;
  f107_workspace *f107_workspace_p;

  double scale;           /* scale factor for error analysis */
} hwm_workspace;

/*
 * Prototypes
 */

hwm_workspace *hwm_alloc(const char *f107_datafile);
void hwm_free(hwm_workspace *w);
void hwm_f107_override(double f107, double f107a, hwm_workspace *w);
int hwm_set_error_scale(double scale, hwm_workspace *w);
size_t hwm_call(double theta, double longitude, time_t t,
             double hstart, double hstep, size_t nalt,
             double *merid, double *zonal, hwm_workspace *w);

/* HWM model prototypes */
void hwm14_(int *iyd, float *sec, float *alt, float *glat, float *glon,
            float *stl, float *f107a, float *f107, float *ap, float *w);

#endif /* INCLUDED_hwm_h */
