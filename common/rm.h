/*
 * rm.h
 */

#ifndef INCLUDED_rm_h
#define INCLUDED_rm_h

typedef struct
{
  double mean; /* current mean */
  double ssq;  /* current sum of squares: ssq = sum_{i=1..n} [x_i - mean]^2 */
  size_t n;
} gsl_statistics_rm_workspace;

gsl_statistics_rm_workspace *gsl_statistics_rm_alloc(void);
void gsl_statistics_rm_free(gsl_statistics_rm_workspace *w);
size_t gsl_statistics_rm_n(gsl_statistics_rm_workspace *w);
int gsl_statistics_rm_add(const double x, gsl_statistics_rm_workspace *w);
double gsl_statistics_rm_mean(gsl_statistics_rm_workspace *w);
double gsl_statistics_rm_variance(gsl_statistics_rm_workspace *w);
double gsl_statistics_rm_sd(gsl_statistics_rm_workspace *w);
double gsl_statistics_rm_sd_mean(gsl_statistics_rm_workspace *w);
int gsl_statistics_rm_reset(gsl_statistics_rm_workspace *w);

#endif /* INCLUDED_rm_h */
