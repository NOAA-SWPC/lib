/*
 * magpole.h
 * Patrick Alken
 */

#ifndef INCLUDED_magpole_h
#define INCLUDED_magpole_h

#include <msynth/msynth.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

typedef struct
{
  double theta_pole; /* pole co-latitude (radians) */
  double phi_pole;   /* pole longitude (radians) */
  size_t N;          /* number of samples for Monte-Carlo */
  msynth_workspace *igrf_workspace_p;
  gsl_rng *rng_p;
} magpole_workspace;

typedef struct
{
  double r;    /* geocentric radius in km */
  double t;    /* decimal year */
  magpole_workspace *w;
} magpole_params;

/*
 * Prototypes
 */

magpole_workspace *magpole_alloc(const size_t N);
void magpole_free(magpole_workspace *w);
int magpole_calc(const double t, const double r, magpole_workspace *w);

#endif /* INCLUDED_magpole_h */
