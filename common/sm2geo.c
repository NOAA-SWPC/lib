/* 
   Procedures for transforming a constant spherical harmonic expansion 
   given in SM coordinates to a time-varying spherical harmonic expansion 
   in GEO coordinates. The underlying equations are described in
   Maus and Luehr, Geophysical Journal International, 2005. These
   routines were developed using the Intel C-compiler, version 8.0

   This file contains 3 functions:

   sm2geo_init()  has to be called once to read in the coefficient file

   sm2geo()       transforms a SH expansion to GEO for a given time.

   sm2geo_green() this routine is required when induction effects are 
                  to be modeled. 
                  It returns the real and imaginary parts of the external
                  fields, as seen in GEO for the given time.

   Stefan Maus, August 2005
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>

#define WD              (2.0 * M_PI)             /* diurnal frequency */
#define WA              (2.0 * M_PI / 365.25)    /* annual frequency */

static FILE *cfile;


/* sub-function converting indices n,m to a single index i */
static int indx(int n,int m)
{
  int i;

  if (m == 0) {
    i = (n - 1) * (n + 1) + 1;
    return i-1;
  }
  if (m >= 0)
    i = (n - 1) * (n + 1) + m * 2;
  else
    i = (n - 1) * (n + 1) - m * 2 + 1;
  return i-1;
}  /*indx*/


/* sub-function giving the order m for the single index i */
static int i2m(int i)
{
  int m;
  m = (i+1) / 2;
  if ((i & 1) == 0)
    m = -m;
  return m;
}  /*i2m*/


/* 
   function sm2geo_init

   reads the file containing the transform coefficients and returns them
   in the array 'coeff' of type 'double'. 

   input:
   char fname[]    name of file containing the sm2geo coefficients, including path
   int maxcoeff    dimension of given array coeff[] (to check whether it
                   is large enough)

   output:
   int *deg_ext    maximum spatial degree of the tranform
   int *deg_day    maximum temporal degree of the daily variation expansion

   returns:
   int ncoeff      actual number of coefficients in transform matrix
 */

int sm2geo_init(char fname[],int *deg_ext,int *deg_day,int maxcoeff,double coeff[])
{
   int i,ism,igeo,iday,ncoeff;
   int n,mgeo,msm,i1,i2;
   int ncoeff_ext,ncoeff_day,ie0,ie1,id,in0,im0,im1;

   cfile = fopen(fname, "r");
   if (cfile == NULL) {printf("%s%s%s\n","file ",fname," could not be opened"); exit(EXIT_FAILURE);}
   fscanf(cfile,"%*[^\n]");
   getc(cfile);
   fscanf(cfile,"%d%d%*[^\n]",deg_ext,deg_day);
   ncoeff_ext = *deg_ext*(*deg_ext+2);
   ncoeff_day = *deg_day*2+1;
   ncoeff = ncoeff_ext * ncoeff_ext * ncoeff_day;
   if (ncoeff > maxcoeff) {printf("array size for coefficients is too small"); exit(EXIT_FAILURE);}
   getc(cfile);
   fscanf(cfile,"%*[^\n]");
   getc(cfile);
   fscanf(cfile,"%*[^\n]");
   getc(cfile);

   i = 0;
   for (n = 1;n<=*deg_ext;n++)
     for (i1 = 0; i1<=2*n; i1++)
       {
         msm = i2m(i1);
         ism = indx(n,msm);
         for (i2 = 0; i2<=2*n; i2++)
           {
             mgeo = i2m(i2);
             igeo = indx(n,mgeo);
             for (iday = 0;iday<ncoeff_day;iday++)
               {
                 fscanf(cfile,"%d%d%d%d",&in0,&im0,&im1,&id);
/*                   printf("sm2geo read: n %d m %d m %d iday %d\n",in0,im0,im1,id); */
                 ie0 = indx(in0,im0);
                 ie1 = indx(in0,im1);
                 if (ie0 != ism || ie1 != igeo || id != iday-*deg_day) 
                   {
                     printf("loop: n %d m %d m %d iday %d\n",n,msm,mgeo,id);
                     printf("read: n %d m %d m %d iday %d\n",in0,im0,im1,id);
                     printf("sm2geo: something is wrong with input file! %d <> %d or %d <> %d or %d <> %d\n",
                            ie0,ism,ie1,igeo,id,iday-*deg_day);
                     exit(EXIT_FAILURE);
                   }
                 fscanf(cfile,"%lg",&coeff[i]);
                 /*                 printf("%12.5f\n",coeff[i]); */
                 i++;
                 fscanf(cfile,"%*[^\n]");
                 getc(cfile);
               } /* for iday */
           } /* for mgeo */
       } /* for msm */
   fclose(cfile);
   return ncoeff;
}

/* 
   function sm2geo

   transforms a SH expansion from SM to GEO for a given time.
   Use routine sm2geo_green for the Earth-induced field. 

   input:
   double fday      time in floating point days since 00:00:00 on 1.1.2000
   int deg_want     spatial degree of the expansion to transformed
   int deg_ext      as returned by sm2geo_init (no choice)
   int deg_day      as returned by sm2geo_init (no choice)
   double coeff[]   as returned by sm2geo_init
   double smcoeff[]  spherical harmonic expansion in SM

   output:
   double geocoeff[]  corresponding spherical harmonic expansion in GEO
                      for the specified time
 */

void sm2geo(double fday,int deg_want,                        
             int deg_ext,int deg_day,double coeff[], 
             double smcoeff[],double geocoeff[])
{
  int ncoeff_want = deg_want*(deg_want+2);
  int ncoeff_ext = deg_ext*(deg_ext+2);
  int ncoeff_day = deg_day*2+1;
  int ncoeff = ncoeff_ext * ncoeff_ext * ncoeff_day;
  int igeo,ism,iday,md,ma,ic;
  int n,mgeo,msm,i1,i2;
  double harm_a,harm_d;

  for (igeo = 0; igeo < ncoeff_ext; igeo++) 
    {
      geocoeff[igeo] = 0.0;
    }

  ic = 0;
  for (n = 1;n<=deg_want;n++)
    for (i1 = 0; i1<=2*n; i1++)
      {
        msm = i2m(i1);
        ism = indx(n,msm);
        for (i2 = 0; i2<=2*n; i2++)
          {
            mgeo = i2m(i2);
            igeo = indx(n,mgeo);
            for (md = -deg_day; md <= deg_day; md++) 
              {
                if (md == 0) harm_d = 1.0;
                if (md > 0)  harm_d = cos(md * WD * fday);
                if (md < 0)  harm_d = sin(md * WD * fday);
                geocoeff[igeo]  += smcoeff[ism] * harm_d * coeff[ic];
                ic++;
                /*                    printf("%d %d %d %d %d %12.5f\n",ism,igeo,md,ic,coeff[ic]); */
              }  /*for md*/
          } /* for mgeo */
      } /* for msm */
  return;
} /*sm2geo*/

/* 
   function sm2geo_green

   transforms a SH expansion from GSM to GEO for a given time.
   If the Earth is assumed to be an ideal conductor, then
   this routine is sufficient. Otherwise, use routine gsm2geo_green

   input:
   double fday      time in floating point days since 00:00:00 on 1.1.2000
   int deg_want     spatial degree of the expansion to transformed
   int deg_ext      as returned by sm2geo_init (no choice)
   int deg_day      as returned by sm2geo_init (no choice)
   double coeff[]   as returned by sm2geo_init

   output:
   array green_day_r    real part of diurnal variation due to unit GM field
   array green_day_i    imaginary part of diurnal variation due to unit GM field
   Note that these arrays have to be defined by the calling program in the same
   dimension as by 'typedef tsm2geo_green[MAX_S_COEFF][MAX_S_COEFF][MAXPERIOD]', below

   The following lines of code show how to then afterwards get the induced field:
 
   for (igeo=0;igeo<MAX_S_COEFF;igeo++) induc[igeo] = 0;       initialize induction SH array
          
   for (n = 1; n<=deg_want; n++)                               loop over spatial degree
     for (msm = -n; msm<=n; msm++)                             loop over order in SM
        {
          ism = indx(n,msm);
          for (mgeo = -n;mgeo<=n;mgeo++)                       loop over order in GEO
            {
              igeo = indx(n,mgeo);
              for (iperiod = 1; iperiod <= deg_day; iperiod++)	    
                induc[igeo] += qfac_day_r[n][iperiod] * sm[ism] * green_day_r[ism][igeo][iperiod]
                             - qfac_day_i[n][iperiod] * sm[ism] * green_day_i[ism][igeo][iperiod];
            }
        } 

    where the qfactors can be taken, e.g. from the table in Maus and Luehr, GJI, 2005.
*/

#define MAX_S_COEFF 15
#define MAXPERIOD 7             /* const, 24h, 12h, 8h */
typedef double tsm2geo_green[MAX_S_COEFF][MAX_S_COEFF][MAXPERIOD];

void sm2geo_green(double fday,int deg_want,                     /* chosen */
                  int deg_ext,int deg_day,double coeff[], /* given */
                  tsm2geo_green green_day_r,tsm2geo_green green_day_i)

{
  int ncoeff_want = deg_want*(deg_want+2);
  int ncoeff_ext = deg_ext*(deg_ext+2);
  int ncoeff_day = deg_day*2+1;
  int ncoeff = ncoeff_ext * ncoeff_ext * ncoeff_day;
  int igeo,ism,iday,md,ic;
  int n,mgeo,msm,i1,i2;
  double harm_d_r,harm_d_i;

  for (ism = 0; ism < ncoeff_ext; ism++) 
    for (igeo = 0; igeo < ncoeff_ext; igeo++)
      { 
        for (iday = 0;iday <= deg_day;iday++) green_day_r[ism][igeo][iday] = 0.0;
        for (iday = 0;iday <= deg_day;iday++) green_day_i[ism][igeo][iday] = 0.0;
      }
    
  ic = 0;
  for (n = 1;n<=deg_want;n++)
    for (i1 = 0; i1<=2*n; i1++)
      {
        msm = i2m(i1);
        ism = indx(n,msm);
        for (i2 = 0; i2<=2*n; i2++)
          {
            mgeo = i2m(i2);
            igeo = indx(n,mgeo);
            for (md = -deg_day; md <= deg_day; md++) 
              {
                if (md == 0)
                  {
                    harm_d_r = 1.0;
                    harm_d_i = 0.0;
                  }
                if (md > 0) 
                  {
                    harm_d_r = cos(md * WD * fday);
                    harm_d_i = sin(md * WD * fday);
                  }
                if (md < 0)
                  {
                    harm_d_r = sin(md * WD * fday);
                    harm_d_i = cos(md * WD * fday);
                  }
                ic++;
                if  (md < 0) iday = -md; else iday = md; 
                green_day_r[ism][igeo][iday] += harm_d_r * coeff[ic-1];
                green_day_i[ism][igeo][iday] += harm_d_i * coeff[ic-1];
                        /*                        printf("%d %d %d %d %d %12.5f\n",ism,igeo,md,ic,coeff[ic-1]); */
              }  /*for md*/
          } /* for mgeo */
      } /* for msm */
  return;
} /* sm2geo_green */
