/*
 * poltor_shell.c
 *
 * Routines related to poloidal field in satellite shell
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "poltor.h"

static double func_An(double s, void *params);
static double func_Bn(double s, void *params);

typedef struct
{
  double n;      /* spherical harmonic degree */
  double j;      /* term in power series */
  double ratio;  /* rmid / r */
} poltor_shell_int_params;

/*
func_An()
  Integrand for A_n^{(j)}(r) integrals

  = s^{n+2} (s - rmid/r)^j
*/

static double
func_An(double s, void *params)
{
  poltor_shell_int_params *p = (poltor_shell_int_params *) params;
  double f;

  f = pow(s, p->n + 2.0) * pow(s - p->ratio, p->j);

  return f;
} /* func_An() */

/*
func_Bn()
  Integrand for B_n^{(j)}(r) integrals

  = s^{1-n} (s - rmid/r)^j
*/

static double
func_Bn(double s, void *params)
{
  poltor_shell_int_params *p = (poltor_shell_int_params *) params;
  double f;

  f = pow(s, 1.0 - p->n) * pow(s - p->ratio, p->j);

  return f;
} /* func_Bn() */

/*
poltor_shell_An()
  Compute \tilde{A}_n^{(j)}(r) integral
*/

int
poltor_shell_An(const size_t n, const size_t j, const double r,
                double *result, poltor_workspace *w)
{
  int s = 0;
  const double ratio1 = r / w->R;
  const double ratio2 = w->rmin / r;
  const double fac = 1.0 / (2.0 * n + 1.0) * pow(ratio1, j + 2.0);
  double intval;

  if (j == 0)
    {
      /* trivial case j = 0 */
      intval = (1.0 - pow(ratio2, n + 3.0)) / (n + 3.0);
    }
  else if (j == 1)
    {
      /* special case j = 1 */
      const double ratio3 = w->rmid / r;
      intval = ((3.0 + n) * (1.0 - pow(ratio2, 4.0 + n)) -
                ratio3*(4.0 + n)*(1.0 - pow(ratio2, 3.0 + n))) /
               ((3.0 + n) * (4.0 + n));
    }
  else
    {
      gsl_function F;
      poltor_shell_int_params params;

      params.n = (double) n;
      params.j = (double) j;
      params.ratio = w->rmid / r;

      F.function = &func_An;
      F.params = &params;

      s = gsl_integration_cquad(&F, ratio2, 1.0, 0.0, 1.0e-8,
                                w->cquad_workspace_p, &intval, NULL, NULL);
    }

  *result = fac * intval;

  return s;
} /* poltor_shell_An() */

/*
poltor_shell_Bn()
  Compute \tilde{B}_n^{(j)}(r) integral
*/

int
poltor_shell_Bn(const size_t n, const size_t j, const double r,
                double *result, poltor_workspace *w)
{
  int s = 0;
  const double ratio1 = r / w->R;
  const double ratio2 = w->rmax / r;
  const double fac = 1.0 / (2.0 * n + 1.0) * pow(ratio1, j + 2.0);
  double intval;

  if (j == 0)
    {
      /* trivial case j = 0 */
      if (n == 2)
        {
          intval = log(ratio2);
        }
      else
        {
          const double y = 2.0 - (double)n;
          intval = (pow(ratio2, y) - 1.0) / y;
        }
    }
  else if (j == 1)
    {
      /* special case j = 1 */
      const double ratio3 = w->rmid / r;

      if (n == 2)
        {
          intval = -1.0 + ratio2 - ratio3 * log(ratio2);
        }
      else if (n == 3)
        {
          intval = ratio3 * (-1.0 + 1.0/ratio2) + log(ratio2);
        }
      else
        {
          intval = -ratio3 * (1.0 - pow(ratio2, 2.0 - n)) * (n - 3.0) +
                   (1.0 - pow(ratio2, 3.0 - n)) * (n - 2.0);
          intval /= (n - 3.0) * (n - 2.0);
        }
    }
  else
    {
      gsl_function F;
      poltor_shell_int_params params;

      params.n = (double) n;
      params.j = (double) j;
      params.ratio = w->rmid / r;

      F.function = &func_Bn;
      F.params = &params;

      s = gsl_integration_cquad(&F, 1.0, ratio2, 0.0, 1.0e-8,
                                w->cquad_workspace_p, &intval, NULL, NULL);
    }

  *result = fac * intval;

  return s;
} /* poltor_shell_Bn() */
