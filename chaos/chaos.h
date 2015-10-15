/*
 * chaos.h
 */

#ifndef INCLUDED_chaos_h
#define INCLUDED_chaos_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "libchaosext.h"

/*
 * Prototypes
 */

int chaos_init(void);
int chaos_term(void);
int chaos_ext(const gsl_vector *tv, const gsl_vector *altv,
              const gsl_vector *latv, const gsl_vector *phiv,
              const gsl_vector *RC_ev, const gsl_vector *RC_iv,
              gsl_matrix *B_ext);

#endif
