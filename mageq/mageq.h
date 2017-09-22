/*
 * $Rev: 4004 $
 * $Date: 2012-01-26 11:51:43 +0100 (Thu, 26 Jan 2012) $
 */
/*
 * mageq.h
 * Patrick Alken
 */

#ifndef INCLUDED_mageq_h
#define INCLUDED_mageq_h

#include <msynth/msynth.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

typedef struct
{
  gsl_min_fminimizer *s;
  msynth_workspace *msynth_workspace_p;
} mageq_workspace;

typedef struct
{
  double r;    /* geocentric radius in km */
  double phi;  /* longitude in radians */
  double t;    /* decimal year */
  mageq_workspace *w;
} mageq_params;

/*
 * Prototypes
 */

mageq_workspace *mageq_alloc();
void mageq_free(mageq_workspace *w);
double mageq_calc(double longitude, double altitude, double t,
                  mageq_workspace *w);
double mageq_angle(double longitude, double altitude, double t,
                   mageq_workspace *w);

#endif /* INCLUDED_mageq_h */
