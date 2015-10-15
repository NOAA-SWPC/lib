/*
 * rm.c - Running mean routines
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include "rm.h"

gsl_statistics_rm_workspace *
gsl_statistics_rm_alloc()
{
  gsl_statistics_rm_workspace *w;

  w = calloc(1, sizeof(gsl_statistics_rm_workspace));

  w->mean = 0.0;
  w->n = 0;

  return w;
} /* gsl_statistics_rm_alloc() */

void
gsl_statistics_rm_free(gsl_statistics_rm_workspace *w)
{
  free(w);
} /* gsl_statistics_rm_free() */

size_t
gsl_statistics_rm_n(gsl_statistics_rm_workspace *w)
{
  return w->n;
} /* gsl_statistics_rm_n() */

/* add a data point to the running totals */
int
gsl_statistics_rm_add(const double x, gsl_statistics_rm_workspace *w)
{
  double delta = x - w->mean;

  w->mean += delta / (w->n + 1);
  w->ssq += delta * (x - w->mean);

  w->n += 1;

  return GSL_SUCCESS;
} /* gsl_statistics_rm_add() */

double
gsl_statistics_rm_mean(gsl_statistics_rm_workspace *w)
{
  return w->mean;
} /* gsl_statistics_rm_mean() */

double
gsl_statistics_rm_variance(gsl_statistics_rm_workspace *w)
{
  return (w->ssq / (w->n - 1));
} /* gsl_statistics_rm_variance() */

double
gsl_statistics_rm_sd(gsl_statistics_rm_workspace *w)
{
  double var = gsl_statistics_rm_variance(w);

  return (sqrt(var));
} /* gsl_statistics_rm_sd() */

/* standard deviation of the mean: sigma / sqrt(n) */
double
gsl_statistics_rm_sd_mean(gsl_statistics_rm_workspace *w)
{
  double sd = gsl_statistics_rm_sd(w);

  return (sd / sqrt((double) w->n));
} /* gsl_statistics_rm_sd_mean() */

int
gsl_statistics_rm_reset(gsl_statistics_rm_workspace *w)
{
  w->mean = 0.0;
  w->ssq = 0.0;
  w->n = 0;

  return GSL_SUCCESS;
} /* gsl_statistics_rm_reset() */
