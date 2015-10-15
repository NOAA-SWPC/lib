/* 
   plmbar.f -- translated by f2c (version 20000121).
   
   2004-10-25. 
   changed intendation; modified to get rid of f2c specials and the 
   'f2c.h' include. m.rother (GFZ Potsdam).
*/

#include <math.h>
#include <stdio.h>

int plmbar_
(
   double *p, 
   double *dp, 
   double *z, 
   int    *lmax, 
   int    *inorm
)
{
   /*

   evaluates normalized associated legendre function p(l,m) as function of 
   z=cos(colatitude) using recurrence relation starting with p(l,l) 
   and then increasing l keeping m fixed.  normalization is: 
   integral(y(l,m)*y(l,m))=4.*pi, where y(l,m) = p(l,m)*exp(i*m*longitude), 
   which is incorporated into the recurrence relation. p(k) contains p(l,m) 
   with k=(l+1)*l/2+m+1; i.e. m increments through range 0 to l before 
   incrementing l. routine is stable in single and real*8 precision to 
   l,m = 511 at least; timing proportional to lmax**2 
   r.j.o'connell 7 sept. 1989 

   a.jackson 19 october 1989  code added at end: 
   (1) choice of normalisation added --- 
       if inorm = 1 schmidt normalisation is chosen, where 
       p[schmidt](l,m) = 1/sqrt(2*l+1)*p[normalised](l,m)
       inorm = 2 for fully normalised harmonics, otherwise error 
   (2) derivatives added and stored in dp(k) 
       using same arrangement as for p(k) 

   -- dimension of p, dp must be (lmax+1)*(lmax+2)/2 in calling program 
   */

   /* Format strings */

  char fmt_99[] = "("\
    "\002plmbar.f: inorm incorrect: "\
    "use inorm=1 for schmidt norm\002,/,\002 "\
    " inorm=2 for full norm.\002"\
    ")";
  
    /* System generated locals */
  int i1, i2;

    /* Local variables */
    double fden, fnum;
    int    k, l, m;
    double f1, f2;
    int    kstart;
    double pm1, pm2, sintsq, fac, plm, pmm, fac1, fac2;

    /* Parameter adjustments */
    --dp;
    --p;

    /* Function Body */
    if ((*lmax < 0) || (fabs(*z) > 1.)) 
    {
       printf("plmbar: bad arguments. lmax = %d, z = %14.4f\n",
              *lmax, *z);
       return(-1);
    }
    if ((*inorm != 1) && (*inorm != 2))
    {
       printf("plmbar: %s\n",fmt_99);
       return(-2);
    }
    /*       --case for p(l,0) */
    pm2 = 1.;
    p[1] = 1.;
    dp[1] = 0.;
    if (*lmax == 0) 
    {
       return 0;
    }
    pm1 = *z;
    p[2] = sqrt(3.) * pm1;
    k = 2;
    i1 = *lmax;
    for (l = 2; l <= i1; ++l) 
    {
       k += l;
       plm = ((double) ((l << 1) - 1) * *z * pm1 - (double) (l - 1)
              * pm2) / (double) l;
       p[k] = sqrt((double) ((l << 1) + 1)) * plm;
       pm2 = pm1;
       /* L4: */
       pm1 = plm;
    }
    /*     --case for m > 0 */
    pmm = 1.;
    sintsq = (1. - *z) * (*z + 1.);
    fnum = -1.;
    fden = 0.;
    kstart = 1;
    i1 = *lmax;
    for (m = 1; m <= i1; ++m) 
    {
       /*     --case for p(m,m) */
       kstart = kstart + m + 1;
       fnum += 2.;
       fden += 2.;
       pmm = pmm * sintsq * fnum / fden;
       pm2 = sqrt((double) ((m << 2) + 2) * pmm);
       p[kstart] = pm2;
       if (m == *lmax) 
       {
          goto L100;
       }
       /*     --case for p(m+1,m) */
       pm1 = *z * sqrt((double) ((m << 1) + 3)) * pm2;
       k = kstart + m + 1;
       p[k] = pm1;
       /*     --case for p(l,m) with l > m+1 */
       if (m < *lmax - 1) 
       {
          i2 = *lmax;
          for (l = m + 2; l <= i2; ++l) 
          {
             k += l;
             f1 = sqrt((double) (((l << 1) + 1) * ((l << 1) - 1)) 
                       / 
                       (
                          double) ((l + m) * (l - m)));
             f2 = sqrt((double) (((l << 1) + 1) * (l - m - 1) * 
                                     (l + m - 1)) 
                       / 
                       (double) (((l << 1) - 3) * (l + m) * 
                                     (l - m)));
             plm = *z * f1 * pm1 - f2 * pm2;
             p[k] = plm;
             pm2 = pm1;
             /* L10: */
             pm1 = plm;
          }
       }
       /* L20: */
    }
  L100:
    /*     choice of normalisation: */
    if (*inorm == 1) 
    {
       k = 1;
       i1 = *lmax;
       for (l = 1; l <= i1; ++l) 
       {
          fac = 1. / sqrt((double) ((l << 1) + 1));
          i2 = l;
          for (m = 0; m <= i2; ++m) 
          {
             ++k;
             p[k] *= fac;
             /* L30: */
          }
       }
    }
    /*     now find derivatives of p(z) wrt theta, where z=cos(theta) */
    /*     recurrence is same regardless of normalisation, since l=constant */
    dp[2] = -p[3];
    dp[3] = p[2];
    k = 3;
    i2 = *lmax;
    for (l = 2; l <= i2; ++l) 
    {
       ++k;
       /*     treat m=0 and m=l separately */
       dp[k] = -sqrt((double) (l * (l + 1)) / 2.) * p[k + 1];
       dp[k + l] = sqrt((double) l / 2.) * p[k + l - 1];
       i1 = l - 1;
       for (m = 1; m <= i1; ++m) 
       {
          ++k;
          fac1 = sqrt((double) ((l - m) * (l + m + 1)));
          fac2 = sqrt((double) ((l + m) * (l - m + 1)));
          if (m == 1) 
          {
             fac2 *= sqrt(2.);
          }
          dp[k] = (fac2 * p[k - 1] - fac1 * p[k + 1]) * .5;
          /* L300: */
       }
       ++k;
       /* L200: */
    }
    return 0;

} 
/* plmbar_ */
