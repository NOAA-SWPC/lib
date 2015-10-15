
/* 
   magfdz.f -- translated by f2c (version 20000121).

   2004-10-25
   modified to get rid of the 'f2c.h' include. rother.
*/

#include <math.h>

int magfdz_
(
   double *p, 
   double *dp, 
   double *theta, 
   double *phi, 
   double *r, 
   int    *lmax, 
   double *g, 
   double *dx, 
   double *dy, 
   double *dz, 
   double *dpot, 
   double *x, 
   double *y, 
   double *z, 
   double *h, 
   double *f, 
   double *ainc, 
   double *d, 
   double *pot
)
{
/*
=========================================================================
   
   j bloxham  8 nov 1982 & 11 oct 1983 
   
   modified version of dg13.sv.fort:magfd & jb62.sv.progs:magfds 
   
   gives field components at radius r 
   
   if igeo=1 then the field elements are rotated to 
   the local geodetic frame at the point specified by 
   (r,theta,phi) in the spherical coordinate system.
   
=========================================================================
   
   this version 16 jan 87 
  
   saves dx dy dz in computation 
   
   stefan: 
   modified to do the potential and thrown out igeo stuff,
   disabled 'ainc'.
=========================================================================
*/
   int    i1, i2;
   double cost, sint, b;
   int k, l, m;
   double t, dpotd, sinth;
   int k1, l1;
   double bb, dxd, dzd, dxy;
   
   /* Parameter adjustments */
   --dpot;
   --dz;
   --dy;
   --dx;
   --g;
   --dp;
   --p;
   
   /* Function Body */

   b     = (double)6371.2 / *r;
   *x    = (double)0.;
   *y    = (double)0.;
   *z    = (double)0.;
   *f    = (double)0.;
   *h    = (double)0.;
   *pot  = (double)0.;

   sinth = sin(*theta);

   if (fabs(sinth) < ((double)1e-10)) 
   {
      sinth = (double)1e-10;
   }
   i1 = *lmax;

   for (l = 1; l <= i1; ++l) 
   {
      l1 = l + 1;
      i2 = l + 2;
      bb = pow(b, (double)i2);/* printf("---- i b = %d %f\n",i2,b);*/
      k  = l * l;
      k1 = (l * l1) / 2 + 1;

      /*
      printf("---- k k1 g[k] dp[k1] p[k1] bb = %4d %4d %f %f %f %f\n", 
             k, k1, g[k], dp[k1], p[k1], bb);
      */

      dx[k] = dp[k1] * bb;
      dy[k] = (double)0.;
      dz[k] = -p[k1] * l1 * bb;

      /*
      printf("----      dx(k) dy(k) dz(k) = %f %f %f\n\n",
             dx[k], dy[k], dz[k]);
      */

      dpot[k]  = p[k1] * *r * bb;
      *x      += g[k]  * dx[k];
      *z      += g[k]  * dz[k];
      *pot    += g[k]  * dpot[k];

      i2 = l;

      for (m = 1; m <= i2; ++m) 
      {
         t = (double) m * *phi;
         k = l * l + (m << 1) - 1;
         k1 = l * l1 / 2 + m + 1;
         sint = sin(t);
         cost = cos(t);
         dxd = dp[k1] * bb;
         dx[k] = dxd * cost;
         dx[k + 1] = dxd * sint;

         *x = *x + g[k] * dx[k] + g[k + 1] * dx[k + 1];

         dxy = m * p[k1] * bb / sinth;
         dy[k] = dxy * sint;
         dy[k + 1] = -dxy * cost;

         *y = *y + g[k] * dy[k] + g[k + 1] * dy[k + 1];

         dzd = -l1 * p[k1] * bb;
         dz[k] = dzd * cost;
         dz[k + 1] = dzd * sint;

         *z = *z + g[k] * dz[k] + g[k + 1] * dz[k + 1];

         dpotd = *r * p[k1] * bb;
         dpot[k] = dpotd * cost;
         dpot[k + 1] = dpotd * sint;
         *pot = *pot + g[k] * dpot[k] + g[k + 1] * dpot[k + 1];
         /*
         printf("------2: x,y,z = %14.4f %14.4f %14.4f\n",*x, *y, *z);
         */
      }
   }
   *h    = sqrt((*x) * (*x) + (*y)* (*y));
   *f    = sqrt((*h) * (*h) + (*z) * (*z));
   *ainc = 0.0; /* not used */
   *d    = atan2(*y, *x);
   
   /*
   printf("-----3: x,y,z,f = %14.4f %14.4f %14.4f %14.4f\n",*x, *y, *z, *f);
   */

   return 0;
}
/* magfdz */

int magfdz_ext_
(
   double *p, 
   double *dp, 
   double *theta, 
   double *phi, 
   double *r, 
   int    *lmax, 
   double *g, 
   double *dx, 
   double *dy, 
   double *dz,
   double *dpot, 
   double *x, 
   double *y, 
   double *z, 
   double *h, 
   double *f, 
   double *ainc, 
   double *d, 
   double *pot
)
{
/*
=========================================================================
   
   same thing for external field.... 
   
   j bloxham  8 nov 1982 & 11 oct 1983 
   
   modified version of dg13.sv.fort:magfd & jb62.sv.progs:magfds 
   
   gives field components at radius r 
   
   if igeo=1 then the field elements are rotated to 
   the local geodetic frame at the point specified by 
   (r,theta,phi) in the spherical coordinate system. 

========================================================================   
   
   this version 16 jan 87 
   
   saves dx dy dz in computation 
   
   stefan:
   modified to do the potential and thrown out igeo stuff,
   disabled 'ainc'.
=======================================================================   
*/
   int    i1, i2;

   double cost, sint, b;
   int    k, l, m;
   double t, dpotd, sinth;
   int    k1, l1;
   double bb, dxd, dzd, dxy;
   
   /* Parameter adjustments */
   --dpot;
   --dz;
   --dy;
   --dx;
   --g;
   --dp;
   --p;
   
   /* Function Body */

   b     = (double)6371.2 / *r;
   *x    = (double)0.;
   *y    = (double)0.;
   *z    = (double)0.;
   *pot  = (double)0.;

   sinth = sin(*theta);

   if (fabs(sinth) < (double)1e-10) 
   {
      sinth = (double)1e-10;
   }
   i1 = *lmax;

   for (l = 1; l <= i1; ++l) 
   {
      l1 = l + 1;
      i2 = -l + 1;
      bb = pow(b, (double)i2);
      k  = l * l;
      k1 = l * l1 / 2 + 1;

      dx[k] = dp[k1] * bb;
      dy[k] = (double)0.;
      dz[k] = l * p[k1] * bb;

      dpot[k]  = p[k1] * *r * bb;
      *x      += g[k]  * dx[k];
      *z      += g[k]  * dz[k];
      *pot    += g[k]  * dpot[k];

      i2 = l;

      for (m = 1; m <= i2; ++m) 
      {
         t = (double) m * *phi;
         k = l * l + (m << 1) - 1;
         k1 = l * l1 / 2 + m + 1;
         sint = sin(t);
         cost = cos(t);
         dxd = dp[k1] * bb;
         dx[k] = dxd * cost;
         dx[k + 1] = dxd * sint;

         *x = *x + g[k] * dx[k] + g[k + 1] * dx[k + 1];

         dxy = m * p[k1] * bb / sinth;
         dy[k] = dxy * sint;
         dy[k + 1] = -dxy * cost;

         *y = *y + g[k] * dy[k] + g[k + 1] * dy[k + 1];

         dzd = l * p[k1] * bb;
         dz[k] = dzd * cost;
         dz[k + 1] = dzd * sint;

         *z = *z + g[k] * dz[k] + g[k + 1] * dz[k + 1];

         dpotd = *r * p[k1] * bb;
         dpot[k] = dpotd * cost;
         dpot[k + 1] = dpotd * sint;
         *pot = *pot + g[k] * dpot[k] + g[k + 1] * dpot[k + 1];
      }
   }

   *h    = sqrt(*x * *x + *y * *y);
   *f    = sqrt(*h * *h + *z * *z);
   *ainc = 0.0; /*not used*/
   *d    = atan2(*y, *x);

   return 0;

} 
/* magfdz_ext */
