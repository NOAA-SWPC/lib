/*
 * ema.c
 *
 * Exponential moving average
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int
ema(const double alpha, double y[], const size_t n)
{
  size_t i;

  for (i = 1; i < n; ++i)
    y[i] = alpha * y[i] + (1.0 - alpha) * y[i - 1];

  return 0;
} /* ema() */

/*
ema_reverse()
  Compute EMA of array, starting from end of array and working
forward
*/

int
ema_reverse(const double alpha, double y[], const size_t n)
{
  size_t i;

  for (i = 2; i <= n; ++i)
    {
      size_t j = n - i;

      y[j] = alpha * y[j] + (1.0 - alpha) * y[j + 1];
    }

  return 0;
} /* ema() */
