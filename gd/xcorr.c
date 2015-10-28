#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static double
xcorr(double x[], double y[], int n, int maxdelay)
{
   int i,j;
   int delay;
   double mx,my,sx,sy,sxy,denom,r;
   double rmax2 = 0.0;
   double rmax = -1.0;
   double rmin = 1.0;
   int imax = 0;
   int imin = 0;
   int sgn = 1;
   
   /* Calculate the mean of the two series x[], y[] */
   mx = 0;
   my = 0;   
   for (i=0;i<n;i++) {
      mx += x[i];
      my += y[i];
   }
   mx /= n;
   my /= n;

   /* Calculate the denominator */
   sx = 0;
   sy = 0;
   for (i=0;i<n;i++) {
      sx += (x[i] - mx) * (x[i] - mx);
      sy += (y[i] - my) * (y[i] - my);
   }
   denom = sqrt(sx*sy);

   /* Calculate the correlation series */
   for (delay=-maxdelay;delay<=maxdelay;delay++) {
      sxy = 0;
      for (i=0;i<n;i++)
        {
          j = i + delay;
          if (j < 0 || j >= n)
            continue;
          else
            sxy += (x[i] - mx) * (y[j] - my);
        }

      /* r is the correlation coefficient at "delay" */
      r = sxy / denom;

      if (delay == 0)
        sgn = GSL_SIGN(r);

      if (r > rmax)
        {
          rmax = r;
          imax = delay;
        }
      if (r < rmin)
        {
          rmin = r;
          imin = delay;
        }
      if (fabs(r) > fabs(rmax2))
        rmax2 = r;

      /*printf("%d %f\n", delay, r);*/
   }

  printf("rmax = %f lag = %d, rmin = %f lag = %d\n", rmax, imax, rmin, imin);

#if 0
  if (sgn == 1)
    return rmax;
  else
    return rmin;
#endif
  return rmax2;
}
