/* 
   Procedures for transforming a constant spherical harmonic expansion 
   given in GSM coordinates to a time-varying spherical harmonic expansion 
   in GEO coordinates. The underlying equations are described in
   Maus and Luehr, Geophysical Journal International, 2005. These
   routines were developed using the Intel C-compiler, version 8.0

   This file contains 3 functions:

   gsm2geo_init()  has to be called once to read in the coefficient file

   gsm2geo()       transforms a SH expansion to GEO for a given time.

   gsm2geo_green() this routine is required when induction effects
                   are to be modeled. 
                   It returns the real and imaginary parts of the external
                   fields, as seen in GEO for the given time.

   Stefan Maus, August 2005
*/

#include <stdio.h>     
#include <stdlib.h>
#include <math.h>

#define WD              (2.0 * M_PI)             /* diurnal frequency */
#define WA              (2.0 * M_PI / 365.25)    /* annual frequency */

#define DEBUG 0                 /* set to 1 to get messages to help debug */

static FILE *cfile;             /* transform coefficient file (read) */


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
   function gsm2geo_init

   reads the file containing the transform coefficients and returns them
   in the array 'coeff' of type 'double'. 

   input:
   char fname[]    name of file containing the gsm2geo coefficients, including path
   int maxcoeff    dimension of given array coeff[] (to check whether it
                   is large enough)

   output:
   int *deg_ext    maximum spatial degree of the tranform
   int *deg_day    maximum temporal degree of the daily variation expansion
   int *deg_ann    maximum temporal degree of the annual variation expansion

   returns:
   int ncoeff      actual number of coefficients in transform matrix
 */

int gsm2geo_init(char fname[],int *deg_ext,int *deg_day,int *deg_ann,int maxcoeff,double coeff[])
{
   int i,igsm,igeo,iday,iann,ncoeff;
   int n,mgeo,mgsm,i1,i2;
   int ncoeff_ext,ncoeff_day,ncoeff_ann,ie0,ie1,id,in0,im0,in1,im1;


   cfile = fopen(fname, "r");
   if (cfile == NULL) {printf("%s%s%s\n","file ",fname," could not be opened"); exit(EXIT_FAILURE);}
   fscanf(cfile,"%*[^\n]");
   getc(cfile);
   fscanf(cfile,"%d%d%d%*[^\n]",deg_ext,deg_day,deg_ann);
   if (DEBUG) printf("gsm2geo: deg_ext = %1d, deg_day = %1d, deg_ann = %1d\n",*deg_ext,*deg_day,*deg_ann);
   ncoeff_ext = *deg_ext*(*deg_ext+2);
   ncoeff_day = *deg_day*2+1;
   ncoeff_ann = *deg_ann*2+1;
   ncoeff = ncoeff_ext * ncoeff_ext * ncoeff_day * ncoeff_ann;
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
         mgsm = i2m(i1);
         igsm = indx(n,mgsm);
         for (i2 = 0; i2<=2*n; i2++)
           {
             mgeo = i2m(i2);
             igeo = indx(n,mgeo);
             for (iday = 0;iday<ncoeff_day;iday++)
               {
                 fscanf(cfile,"%d%d%d%d",&in0,&im0,&im1,&id);
                 if (DEBUG) printf("gsm read: n %d m %d m %d iday %d\n",in0,im0,im1,id);
                 ie0 = indx(in0,im0);
                 ie1 = indx(in0,im1);
                 if (ie0 != igsm || ie1 != igeo || id != iday-*deg_day) 
                   {
                     printf("loop: n %d m %d m %d iday %d\n",n,mgsm,mgeo,id);
                     printf("read: n %d m %d m %d iday %d\n",in0,im0,im1,id);
                     printf("gsm2geo: something is wrong with input file! %d <> %d or %d <> %d or %d <> %d\n",
                            ie0,igsm,ie1,igeo,id,iday-*deg_day);
                     exit(EXIT_FAILURE);
                   }
                 for(iann = 0;iann<ncoeff_ann;iann++)
                   {
                     fscanf(cfile,"%lg",&coeff[i]);
                     /*                 printf("%12.5f\n",coeff[i]); */
                     i++;
                   }
                 fscanf(cfile,"%*[^\n]");
                 getc(cfile);
               } /* for iday */
           } /* for mgeo */
       } /* for mgsm */
   fclose(cfile);
   return ncoeff;
}


/* 
   function gsm2geo

   transforms a SH expansion from GSM to GEO for a given time.
   Use routine gsm2geo_green for the Earth-induced field. 

   input:
   double fday      time in floating point days since 00:00:00 on 1.1.2000
   int deg_want     spatial degree of the expansion to transformed
   int deg_ext      as returned by gsm2geo_init (no choice)
   int deg_day      as returned by gsm2geo_init (no choice)
   int deg_ann      as returned by gsm2geo_init (no choice)
   double coeff[]   as returned by gsm2geo_init
   double gsmcoeff[]  spherical harmonic expansion in GSM

   output:
   double geocoeff[]  corresponding spherical harmonic expansion in GEO
                      for the specified time
 */

void gsm2geo(double fday,int deg_want,                             
             int deg_ext,int deg_day,int deg_ann,double coeff[],   
             double gsmcoeff[],double geocoeff[])
{
  int ncoeff_want = deg_want*(deg_want+2);
  int ncoeff_ext = deg_ext*(deg_ext+2);
  int ncoeff_day = deg_day*2+1;
  int ncoeff_ann = deg_ann*2+1;
  int ncoeff = ncoeff_ext * ncoeff_ext * ncoeff_day * ncoeff_ann;
  int igeo,igsm,iday,md,ma,ic;
  int n,mgeo,mgsm,i1,i2;
  double harm_a,harm_d;

  if (DEBUG) 
    for (igsm=0;igsm<ncoeff_want;igsm++) 
      { 
        printf("gsmcoeff[%1d] = %12.3f ",igsm,gsmcoeff[igsm]); 
        printf("\n"); 
      } 


  for (igeo = 0; igeo < ncoeff_ext; igeo++) 
    {
      geocoeff[igeo] = 0.0;
    }


  ic = 0;
  for (n = 1;n<=deg_want;n++)
    for (i1 = 0; i1<=2*n; i1++)
      {
        mgsm = i2m(i1);
        igsm = indx(n,mgsm);
        for (i2 = 0; i2<=2*n; i2++)
          {
            mgeo = i2m(i2);
            igeo = indx(n,mgeo);
            for (md = -deg_day; md <= deg_day; md++) 
              {
                if (md == 0) harm_d = 1.0;
                if (md > 0)  harm_d = cos(md * WD * fday);
                if (md < 0)  harm_d = sin(md * WD * fday);
                for (ma = -deg_ann; ma <= deg_ann; ma++) 
                  {
                    if (ma == 0) harm_a = 1.0;
                    if (ma > 0)  harm_a = cos(ma * WA * fday);
                    if (ma < 0)  harm_a = sin(ma * WA * fday);
                    geocoeff[igeo]  += gsmcoeff[igsm] * harm_d * harm_a * coeff[ic];
                    ic++;
                    if (DEBUG) 
                      printf("%d %d %d %d %d %12.5f\n",igsm,igeo,md,ma,ic,coeff[ic]);
                  }  /*for ma*/
              }  /*for md*/
          } /* for mgeo */
      } /* for mgsm */
  return;
} /*gsm2geo*/


/* 
   function gsm2geo_green

   transforms a SH expansion from GSM to GEO for a given time.
   If the Earth is assumed to be an ideal conductor, then
   this routine is sufficient. Otherwise, use routine gsm2geo_green

   input:
   double fday      time in floating point days since 00:00:00 on 1.1.2000
   int deg_want     spatial degree of the expansion to transformed
   int deg_ext      as returned by gsm2geo_init (no choice)
   int deg_day      as returned by gsm2geo_init (no choice)
   int deg_ann      as returned by gsm2geo_init (no choice)
   double coeff[]   as returned by gsm2geo_init

   output:
   array green_day_r    real part of diurnal variation due to unit GSM field
   array green_day_i    imaginary part of diurnal variation due to unit GSM field
   array green_ann_r    real part of annual variation due to unit GSM field
   array green_ann_i    imaginary part of annual variation due to constant GSM field
   Note that these arrays have to be defined by the calling program in the same
   dimension as by 'typedef tgsm2geo_green[MAX_S_COEFF][MAX_S_COEFF][MAXPERIOD]', below

   The following lines of code show how to then afterwards get the induced field:
 
   for (igeo=0;igeo<MAX_S_COEFF;igeo++) induc[igeo] = 0;       initialize induction SH array
          
   for (n = 1; n<=deg_want; n++)                                   loop over spatial degree
     for (mgsm = -n; mgsm<=n; mgsm++)                              loop over order in GSM
        {
          igsm = indx(n,mgsm);
          for (mgeo = -n;mgeo<=n;mgeo++)                           loop over order in GEO
            {
              igeo = indx(n,mgeo);
              for (iperiod = 1; iperiod <= deg_day; iperiod++)	    
                induc[igeo] += qfac_day_r[n][iperiod] * gsm[igsm] * green_day_r[igsm][igeo][iperiod]
                             - qfac_day_i[n][iperiod] * gsm[igsm] * green_day_i[igsm][igeo][iperiod];

              for (iperiod = 1; iperiod <= deg_ann; iperiod++)     
                induc[igeo] += qfac_ann_r[n][iperiod] * gsm[igsm] * green_ann_r[igsm][igeo][iperiod]
                             - qfac_ann_i[n][iperiod] * gsm[igsm] * green_ann_i[igsm][igeo][iperiod];
            }
        } 

    where the qfactors can be taken, e.g. from the table in Maus and Luehr, GJI, 2005.
*/

#define MAX_S_COEFF 15
#define MAXPERIOD 7             /* const, 24h, 12h, 8h */

typedef double tgsm2geo_green[MAX_S_COEFF][MAX_S_COEFF][MAXPERIOD];

void gsm2geo_green(double fday,int deg_want,                           /* chosen */
                   int deg_ext,int deg_day,int deg_ann,double coeff[], /* given */
                   tgsm2geo_green green_day_r,tgsm2geo_green green_day_i,
                   tgsm2geo_green green_ann_r,tgsm2geo_green green_ann_i)
{
  int ncoeff_want = deg_want*(deg_want+2);
  int ncoeff_ext = deg_ext*(deg_ext+2);
  int ncoeff_day = deg_day*2+1;
  int ncoeff_ann = deg_ann*2+1;
  int ncoeff = ncoeff_ext * ncoeff_ext * ncoeff_day * ncoeff_ann;
  int igeo,igsm,iday,iann,md,ma,ic;
  int n,mgeo,mgsm,i1,i2;
  double harm_a_r,harm_d_r,harm_a_i,harm_d_i;

  if (DEBUG) 
    printf ("calling gsm2geo with fday: %.2f, deg_want: %1d, deg_ext: %1d, deg_day: %1d, deg_ann: %1d\n",
            fday, deg_want, deg_ext, deg_day, deg_ann);
                   
  for (igsm = 0; igsm < ncoeff_ext; igsm++) 
    for (igeo = 0; igeo < ncoeff_ext; igeo++)
      { 
        for (iday = 0;iday <= deg_day;iday++) green_day_r[igsm][igeo][iday] = 0.0;
        for (iann = 0;iann <= deg_ann;iann++) green_ann_r[igsm][igeo][iann] = 0.0;
        for (iday = 0;iday <= deg_day;iday++) green_day_i[igsm][igeo][iday] = 0.0;
        for (iann = 0;iann <= deg_ann;iann++) green_ann_i[igsm][igeo][iann] = 0.0;
      }
    
  ic = 0;
  for (n = 1;n<=deg_want;n++)
    for (i1 = 0; i1<=2*n; i1++)
      {
        mgsm = i2m(i1);
        igsm = indx(n,mgsm);
        for (i2 = 0; i2<=2*n; i2++)
          {
            mgeo = i2m(i2);
            igeo = indx(n,mgeo);
            if (DEBUG) printf("mgsm: %d, mgeo: %d\n",mgsm,mgeo);
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
                for (ma = -deg_ann; ma <= deg_ann; ma++) 
                  {
                    if (ma == 0) 
                      {
                        harm_a_r = 1.0;
                        harm_a_i = 0.0;
                      }
                    if (ma > 0)
                      {
                        harm_a_r = cos(ma * WA * fday);
                        harm_a_i = sin(ma * WA * fday);
                      }
                    if (ma < 0)
                      {
                        harm_a_r = sin(ma * WA * fday);
                        harm_a_i = cos(ma * WA * fday);
                      }
                    ic++;
                    if  (md < 0) iday = -md; else iday = md;
                    green_day_r[igsm][igeo][iday] += harm_d_r * harm_a_r * coeff[ic-1];
                    green_day_i[igsm][igeo][iday] += harm_d_i * harm_a_r * coeff[ic-1];
                    if (ma < 0) iann = -ma; else iann = ma; 
                    if (md == 0)  /* daily-varying terms average to zero, annually */
                      {
                        green_ann_r[igsm][igeo][iann] += harm_d_r * harm_a_r * coeff[ic-1];
                        green_ann_i[igsm][igeo][iann] += harm_d_r * harm_a_i * coeff[ic-1];
                        if (DEBUG)
                          if (igsm==0 && igeo ==2) printf("%d %d %d %d %d %12.5f %12.5f %12.5f %12.5f\n",igsm,igeo,md,ma,ic, 
                                                          harm_a_r,harm_a_i,
                                                          green_ann_r[igsm][igeo][iann],green_ann_i[igsm][igeo][iann] );
                      }
                  }  /*for ma*/
              }  /*for md*/
          } /* for mgeo */
      } /* for mgsm */
/*   if (DEBUG) getchar(); */
  if (DEBUG)
    {
      printf ("green_day_r[igsm=0, igeo=0, iday=0]: %.3f\n", green_day_r[0][0][0]);
      printf ("green_day_r[igsm=0, igeo=1, iday=1]: %.3f\n", green_day_r[0][1][1]);
      printf ("green_ann_r[igsm=0, igeo=0, iann=0]: %.3f\n", green_ann_r[0][0][0]);
      printf ("green_ann_r[igsm=3, igeo=3, iann=2]: %.3f\n", green_ann_r[3][3][2]);
    }
  return;
} /* gsm2geo_green */
