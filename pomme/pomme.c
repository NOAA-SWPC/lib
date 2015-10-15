#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


// Linux
#include "./POMME8_core.h"
#include "./POMME8_ext.h"
#include "./POMME8_crust.h"
#include "./POMME8_gsm2geo.h"
#include "./POMME8_sm2geo.h"
#include "./POMME8_qfac.h"
#include "./pomme_mod.h"
#include "./Gdef.h"

#define TEST 0
#define MIN_POLE_DIST   (0.001 * M_PI/180.0) /* Sperical coordinates not defined closer to poles */

extern int imin(int i1, int i2);

extern void get_ext(double g1, double g2, double g3, 
                    double phi, double theta, double rrel, 
                    double *x, double  *y, double  *z); 
extern void get_int(double g1, double g2, double g3,
                    double phi, double theta, double rrel, 
                    double *x, double *y, double *z); 
extern void plmbar_(double *p,double *dp,double *costh,int *degree,int *inorm);
extern void magfdz_(double *p,double *dp,double *theta,double *phi,double *rad,
                    int *degree,double *g,double *dx,double *dy,double *dz,
                    double *dpot,double *x,double *y,double *z,double *h,
                    double *f,double *ainc,double *d,double *pot);
extern void magfdz_ext_(double *p,double *dp,double *theta,double *phi,double *rad,
                    int *degree,double *g,double *dx,double *dy,double *dz,
                    double *dpot,double *x,double *y,double *z,double *h,
                    double *f,double *ainc,double *d,double *pot);

typedef double tgsm2geo_green[MAX_S_COEFF][MAX_S_COEFF][MAXPERIOD];

extern void sm2geo(double fday,int deg_wanted,                                          /* chosen */
                   int deg_ext,int deg_day,double coeff[SM2GEO_MAXNCOEFF],              /* given */
                   double smcoeff[MAX_EXT_NCOEFF],double geocoeff[MAX_EXT_NCOEFF]);

extern void sm2geo_green(double fday,int deg_want,                               /* chosen */
                         int deg_ext,int deg_day,double coeff[SM2GEO_MAXNCOEFF], /* given */
                         tgsm2geo_green green_day_r,tgsm2geo_green green_day_i);

extern void gsm2geo(double fday,int deg_wanted,                                               /* chosen */
                    int deg_ext,int deg_day,int deg_ann,double coeff[GSM2GEO_MAXNCOEFF],      /* given */
                    double gsmcoeff[MAX_EXT_NCOEFF],double geocoeff[MAX_EXT_NCOEFF]);

extern void gsm2geo_green(double fday,int deg_want,                                           /* chosen */
                          int deg_ext,int deg_day,int deg_ann,double coeff[SM2GEO_MAXNCOEFF], /* given */
                          tgsm2geo_green green_day_r,tgsm2geo_green green_day_i,
                          tgsm2geo_green green_ann_r,tgsm2geo_green green_ann_i);

extern void geodetic2geocentric(double phi, double h,double *latrad, double *r); /* phi is geodetic latitude in radian*/

extern void geocentric2geodetic_vec(double theta, double delta, double bx, double by, double bz,
                                    double *x, double *y, double *z);


int imin(int x,int y){ if (x < y) return x; else return y; }
int imax(int x,int y){ if (x > y) return x; else return y; }

/* psa */
double GLOBAL_Q[3]; /* save geocentric SH degree 1 representation of external field */

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


void pomme_get_internal_coeff ( double fyear, int *deg, double  *g)
{
  int i,iyear;

  double dyear;

  *deg = POMME_DEG_CRUST;

  iyear = imin(POMME_LAST_YEAR - POMME_FIRST_YEAR,(int) fyear - POMME_FIRST_YEAR);
  dyear = fyear-iyear - POMME_FIRST_YEAR;
  
  for (i=0; i<POMME_NCOEFF_MAIN; i++)
    g[i] = POMME_Core[iyear][0][i] + dyear*POMME_Core[iyear][1][i];
  for (i=POMME_NCOEFF_MAIN; i<POMME_NCOEFF_CRUST; i++)
    g[i] = POMME_Crust[i];
  
  return;
} /* pomme_get_coeff */


void pomme_internal_geocentric (int limit_deg, /* 0 = no limit */
                      double fyear, double rad, double theta, double phi,
                      double *mx, double *my, double *mz, 
                      double *cx, double *cy, double *cz)
{
  int i,iyear;
  int inorm = 1;	                   /*1=Schmid, 2=fully*/
  int deg_main = POMME_DEG_MAIN;
  int deg_crust = POMME_DEG_CRUST;
  int np, nc;
  double h,ainc,d,pot,f;
  double costh,dyear;
  double mg[POMME_NCOEFF_MAIN];
  double *p, *dp, *dx, *dy, *dz, *dpot;

  if (limit_deg > 0) deg_crust = limit_deg;

  np = (deg_crust+1) * (deg_crust+2) / 2;
  nc = deg_crust * (deg_crust+2);
  p = malloc(np*sizeof(double));
  dp = malloc(np*sizeof(double));
  dx = malloc(nc*sizeof(double));
  dy = malloc(nc*sizeof(double));
  dz = malloc(nc*sizeof(double));
  dpot = malloc(nc*sizeof(double));

  iyear = imin(POMME_LAST_YEAR - POMME_FIRST_YEAR,(int) fyear - POMME_FIRST_YEAR);
  dyear = fyear-iyear - POMME_FIRST_YEAR;
  
  for (i=0; i<POMME_NCOEFF_MAIN; i++)
    mg[i] = POMME_Core[iyear][0][i] + dyear*POMME_Core[iyear][1][i];
  
  costh = cos(theta);
  plmbar_(&p[0],&dp[0],&costh,&deg_crust,&inorm);
  
  magfdz_(&p[0],&dp[0],&theta,&phi,&rad,&deg_main,&mg[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],mx,my,mz,&h,&f,&ainc,&d,&pot);
  
  magfdz_(&p[0],&dp[0],&theta,&phi,&rad,&deg_crust,&POMME_Crust[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],cx,cy,cz,&h,&f,&ainc,&d,&pot);
  
  if (fabs(M_PI-theta) < MIN_POLE_DIST || fabs(theta) < MIN_POLE_DIST)
    {
      *mx = 0;                /* not defined at poles*/
      *my = 0;                /* not defined at poles*/
      *cx = 0;                /* not defined at poles*/
      *cy = 0;                /* not defined at poles*/
    }

  free(p);
  free(dp);
  free(dx);
  free(dy);
  free(dz);
  free(dpot);

  return;
} /* pomme_internal_geocentric */

void pomme_external_geocentric (int control[MAXCONTROL],
                                double fday, double rrel, double theta, double phi,
                                double est, double ist, 
                                double imf_by,double f107, double Em,  /* time lagged but not peaked or subtracted mean */
                                double *x, double *y, double *z)
{
  int igeo,igsm,ism,iperiod,gsm2geo_deg,sm2geo_deg;
  double tx,ty,tz,xi,yi,zi,xe,ye,ze,fac, Em_eff, f107_eff;
  double est_x,ist_x,est_y,ist_y,est_z,ist_z;
  double geo_ext_field[MAX_EXT_NCOEFF],geo_ind_field[MAX_EXT_NCOEFF];
  double induc_sm[MAX_EXT_NCOEFF],induc_gsm[MAX_EXT_NCOEFF],induc_imfby[MAX_EXT_NCOEFF];
  double est_in_geo[MAX_EXT_NCOEFF],ist_in_geo[MAX_EXT_NCOEFF];
  double sm_est_field[MAX_EXT_NCOEFF],sm_ist_field[MAX_EXT_NCOEFF],sm_f107_field[MAX_EXT_NCOEFF],gsm_imfby_field[MAX_EXT_NCOEFF],gsm_Em_field[MAX_EXT_NCOEFF];
  double sm_in_geo[MAX_EXT_NCOEFF],gsm_in_geo[MAX_EXT_NCOEFF],imfby_in_geo[MAX_EXT_NCOEFF];
  double Em_in_geo[MAX_EXT_NCOEFF],induc_Em[MAX_EXT_NCOEFF],f107_in_geo[MAX_EXT_NCOEFF],induc_f107[MAX_EXT_NCOEFF];

  tgsm2geo_green gsm_green_day_r, gsm_green_ann_r, gsm_green_day_i, gsm_green_ann_i;
  tgsm2geo_green sm_green_day_r, sm_green_day_i;
  int ext_deg = MAX_EXT_DEG;
  double costh,t;
  int inorm = 1;	/*1=schmid, 2=fully*/
  double h,ainc,d,pot,f,rad;
  double g[MAX_EXT_NCOEFF];
  double p[(MAX_EXT_DEG+1)*(MAX_EXT_DEG+2)/2],dp[(MAX_EXT_DEG+1)*(MAX_EXT_DEG+2)/2],dx[MAX_EXT_DEG*(MAX_EXT_DEG+2)];
  double dy[MAX_EXT_DEG*(MAX_EXT_DEG+2)],dz[MAX_EXT_DEG*(MAX_EXT_DEG+2)],dpot[MAX_EXT_DEG*(MAX_EXT_DEG+2)];
  int n,mgsm,msm,mgeo;
#define EST_IST_via_SM2GEO 0    /* The SM2GEO method is not as accurate. Difference is up to 1 nT for large Dst */

  Em_eff = EM_PEAK * Em / sqrt(EM_PEAK * EM_PEAK + Em * Em) - EM_MEAN;
  f107_eff = f107 - F107_MEAN;
 
  sm2geo_deg = control[POS_SM_IND];
  if (control[POS_F107] > sm2geo_deg) sm2geo_deg = control[POS_F107]; 

  sm2geo_green(fday, sm2geo_deg, SM2GEO_DEG_EXT, SM2GEO_DEG_DAY, coeff_sm2geo,
               sm_green_day_r, sm_green_day_i);

  gsm2geo_deg = control[POS_GSM_IND];
  if (control[POS_IMFBY_IND] > gsm2geo_deg) gsm2geo_deg = control[POS_IMFBY_IND];
  if (control[POS_F107_IND] > gsm2geo_deg) gsm2geo_deg = control[POS_F107_IND];

  gsm2geo_green(fday, gsm2geo_deg, GSM2GEO_DEG_EXT, GSM2GEO_DEG_DAY, GSM2GEO_DEG_ANN, coeff_gsm2geo,
                gsm_green_day_r, gsm_green_day_i, gsm_green_ann_r, gsm_green_ann_i);


  for (igeo=0;igeo<MAX_S_COEFF;igeo++) 
    {
      geo_ext_field[igeo] = 0;
      geo_ind_field[igeo] = 0;
    }

  if (TEST)
    for (igeo=0;igeo<MAX_S_COEFF;igeo++) 
      {
        sm_in_geo[igeo] = 0;
        induc_sm[igeo] = 0;
        est_in_geo[igeo] = 0;
        ist_in_geo[igeo] = 0;
        gsm_in_geo[igeo] = 0;
        induc_gsm[igeo] = 0;
        imfby_in_geo[igeo] = 0;
        induc_imfby[igeo] = 0;
        Em_in_geo[igeo] = 0;
        induc_Em[igeo] = 0;
        f107_in_geo[igeo] = 0;
        induc_f107[igeo] = 0;
      }

                                /* SM */

  if (control[POS_SM] > 0)
    {
      sm2geo(fday, control[POS_SM], SM2GEO_DEG_EXT, SM2GEO_DEG_DAY, coeff_sm2geo,
             POMME_C_SM,sm_in_geo);

      for (igeo=0;igeo<control[POS_SM]*(control[POS_SM]+2);igeo++)
        geo_ext_field[igeo] +=  sm_in_geo[igeo];
    }

                                /* Induced by SM in GEO */

  if (control[POS_SM_IND] > 0)
    {      
      for (igeo=0;igeo<MAX_S_COEFF;igeo++) induc_sm[igeo] = 0;
      
      for (n = 1;n<=control[POS_SM_IND];n++)
        for (msm = -n;msm<=n;msm++)
          {
            ism = indx(n,msm);
            for (mgeo = -n;mgeo<=n;mgeo++)
              {
                igeo = indx(n,mgeo);
                for (iperiod = 1; iperiod <= SM2GEO_DEG_DAY; iperiod++)	      /*green[,,0] is constant term*/
                  induc_sm[igeo] += qfac_day_r[n-1][iperiod - 1] * POMME_C_SM[ism] * sm_green_day_r[ism][igeo][iperiod]
                    - CONVENTION * qfac_day_i[n-1][iperiod - 1] * POMME_C_SM[ism] * sm_green_day_i[ism][igeo][iperiod];
              } /* for mgeo */
          } /* for msm */
      for (igeo=0;igeo<MAX_S_COEFF;igeo++) geo_ind_field[igeo] += induc_sm[igeo];
    } /* if sm induced */

                                /* EST/IST field */

  if (fabs(est-9999.0) > 1.0e-10)
    {
      if (control[POS_EST_FAC] > 0) fac = EST_FAC; else fac = 1.0;
#if EST_IST_via_SM2GEO

      if (control[POS_EST] > 0)
        {
          sm_est_field[0] = -fac * est;
          sm_est_field[1] = 0.0;
          sm_est_field[2] = 0.0;
          sm2geo(fday, 1, SM2GEO_DEG_EXT, SM2GEO_DEG_DAY, coeff_sm2geo,
                 sm_est_field,est_in_geo);
          for (igeo=0;igeo<3;igeo++)
            geo_ext_field[igeo] +=  est_in_geo[igeo];
        }
      if (control[POS_IST] > 0)
        {
          sm_ist_field[0] = -fac * ist;
          sm_ist_field[1] = 0.0;
          sm_ist_field[2] = 0.0;
          sm2geo(fday, 1, SM2GEO_DEG_EXT, SM2GEO_DEG_DAY, coeff_sm2geo,
                 sm_ist_field,ist_in_geo);
          for (igeo=0;igeo<3;igeo++)
            geo_ind_field[igeo] +=  ist_in_geo[igeo];
        }
#else
      if (control[POS_EST] > 0)
        get_ext(fac*est*G10/MOM, fac*est*G11/MOM, fac*est*H11/MOM, phi, theta, rrel, &est_x, &est_y, &est_z);
      else {est_x = 0; est_y = 0; est_z = 0;}
      if (control[POS_IST] > 0)
        get_int(fac*ist*G10/MOM, fac*ist*G11/MOM, fac*ist*H11/MOM, phi, theta, rrel, &ist_x, &ist_y, &ist_z);
      else {ist_x = 0; ist_y = 0; ist_z = 0;}
#endif
    } /* if not missing value for dst */

                                /* F107 */
  for (ism=0;ism<control[POS_F107]*(control[POS_F107]+2);ism++)
    sm_f107_field[ism] =  f107_eff*POMME_C_F107[ism];

  if (control[POS_F107] > 0)
    {
      sm2geo(fday, control[POS_F107], SM2GEO_DEG_EXT, SM2GEO_DEG_DAY, coeff_sm2geo,
              sm_f107_field,f107_in_geo);
      for (igeo=0;igeo<control[POS_F107]*(control[POS_F107]+2);igeo++)
        geo_ext_field[igeo] +=  f107_in_geo[igeo];
    }


                                /* F107 induced */
  if (control[POS_F107_IND] > 0)
    {
      for (igeo=0;igeo<MAX_S_COEFF;igeo++) induc_f107[igeo] = 0;

      for (n = 1;n<=control[POS_F107_IND];n++)
        for (msm = -n;msm<=n;msm++)
          {
            ism = indx(n,msm);
            for (mgeo = -n;mgeo<=n;mgeo++)
              {
                igeo = indx(n,mgeo);
                for (iperiod = 1; iperiod <= SM2GEO_DEG_DAY; iperiod++) 
                  induc_f107[igeo] += POMME_C_F107[ism] * f107_eff * (qfac_day_r[n-1][iperiod - 1] * sm_green_day_r[ism][igeo][iperiod]
                                                                        - CONVENTION * qfac_day_i[n-1][iperiod - 1] * sm_green_day_i[ism][igeo][iperiod]);
              } /* for mgeo */
          } /* for msm */
      for (igeo=0;igeo<MAX_S_COEFF;igeo++) geo_ind_field[igeo] += induc_f107[igeo];
    }


                                /* GSM */
  if (control[POS_GSM] > 0)
    {
      gsm2geo(fday, control[POS_GSM], GSM2GEO_DEG_EXT, GSM2GEO_DEG_DAY, GSM2GEO_DEG_ANN, coeff_gsm2geo,
              POMME_C_GSM,gsm_in_geo);
      for (igeo=0;igeo<control[POS_GSM]*(control[POS_GSM]+2);igeo++)
        geo_ext_field[igeo] +=  gsm_in_geo[igeo];
    }

                                /* Induced by GSM in GEO */

  if (control[POS_GSM_IND] > 0)
    {
      for (igeo=0;igeo<MAX_S_COEFF;igeo++) induc_gsm[igeo] = 0;
          
      for (n = 1;n<=control[POS_GSM_IND];n++)
        for (mgsm = -n;mgsm<=n;mgsm++)
          {
            igsm = indx(n,mgsm);
            for (mgeo = -n;mgeo<=n;mgeo++)
              {
                igeo = indx(n,mgeo);
                for (iperiod = 1; iperiod <= GSM2GEO_DEG_DAY; iperiod++)	  /*green[,,0] is constant term*/ 
                  induc_gsm[igeo] += qfac_day_r[n-1][iperiod - 1] * POMME_C_GSM[igsm] * gsm_green_day_r[igsm][igeo][iperiod]
                    - CONVENTION * qfac_day_i[n-1][iperiod - 1] * POMME_C_GSM[igsm] * gsm_green_day_i[igsm][igeo][iperiod];
                for (iperiod = 1; iperiod <= GSM2GEO_DEG_ANN; iperiod++)    /*green[,,0] is constant term*/ 
                  induc_gsm[igeo] +=   qfac_ann_r[n-1][iperiod - 1] * POMME_C_GSM[igsm] * gsm_green_ann_r[igsm][igeo][iperiod]
                    - CONVENTION * qfac_ann_i[n-1][iperiod - 1] * POMME_C_GSM[igsm] * gsm_green_ann_i[igsm][igeo][iperiod];
              } /* for mgeo */
          } /* for mgsm */
      for (igeo=0;igeo<MAX_S_COEFF;igeo++) geo_ind_field[igeo] += induc_gsm[igeo];
    } /* if gsm induced */


                                /* IMF-By */
  for (igsm=0;igsm<control[POS_IMFBY]*(control[POS_IMFBY]+2);igsm++) /* deg*(deg+2) */
    gsm_imfby_field[igsm] =  imf_by*POMME_C_IMFBY[igsm];

  if (control[POS_IMFBY] > 0)
    {
      gsm2geo(fday, control[POS_IMFBY], GSM2GEO_DEG_EXT, GSM2GEO_DEG_DAY, GSM2GEO_DEG_ANN, coeff_gsm2geo,
              gsm_imfby_field,imfby_in_geo);
      for (igeo=0;igeo<control[POS_IMFBY]*(control[POS_IMFBY]+2);igeo++)
        geo_ext_field[igeo] +=  imfby_in_geo[igeo];
    }


                                /* IMF-By induced */
  if (control[POS_IMFBY_IND] > 0)
    {
      for (igeo=0;igeo<MAX_S_COEFF;igeo++) induc_imfby[igeo] = 0;

      for (n = 1;n<=control[POS_IMFBY_IND];n++)
        for (mgsm = -n;mgsm<=n;mgsm++)
          {
            igsm = indx(n,mgsm);
            for (mgeo = -n;mgeo<=n;mgeo++)
              {
                igeo = indx(n,mgeo);
                for (iperiod = 1; iperiod <= GSM2GEO_DEG_DAY; iperiod++) 
                  induc_imfby[igeo] += POMME_C_IMFBY[igsm] * imf_by * (qfac_day_r[n-1][iperiod - 1] * gsm_green_day_r[igsm][igeo][iperiod]
                                                                        - CONVENTION * qfac_day_i[n-1][iperiod - 1] * gsm_green_day_i[igsm][igeo][iperiod]);
                for (iperiod = 1; iperiod <= GSM2GEO_DEG_ANN; iperiod++) 
                  induc_imfby[igeo] += POMME_C_IMFBY[igsm] * imf_by * (qfac_ann_r[n-1][iperiod - 1] * gsm_green_ann_r[igsm][igeo][iperiod]
                                                                        - CONVENTION * qfac_ann_i[n-1][iperiod - 1] * gsm_green_ann_i[igsm][igeo][iperiod]);
              } /* for mgeo */
          } /* for mgsm */
      for (igeo=0;igeo<MAX_S_COEFF;igeo++) geo_ind_field[igeo] += induc_imfby[igeo];
    } /* IMF-By*/


                                /* Em */
  for (igsm=0;igsm<control[POS_EM]*(control[POS_EM]+2);igsm++)
    gsm_Em_field[igsm] =  Em_eff*POMME_C_EM[igsm];

  if (control[POS_EM] > 0)
    {
      gsm2geo(fday, control[POS_EM], GSM2GEO_DEG_EXT, GSM2GEO_DEG_DAY, GSM2GEO_DEG_ANN, coeff_gsm2geo,
              gsm_Em_field,Em_in_geo);
      for (igeo=0;igeo<control[POS_EM]*(control[POS_EM]+2);igeo++)
        geo_ext_field[igeo] +=  Em_in_geo[igeo];
    }


                                /* Em induced */
  if (control[POS_EM_IND] > 0)
    {
      for (igeo=0;igeo<MAX_S_COEFF;igeo++) induc_Em[igeo] = 0;

      for (n = 1;n<=control[POS_EM_IND];n++)
        for (mgsm = -n;mgsm<=n;mgsm++)
          {
            igsm = indx(n,mgsm);
            for (mgeo = -n;mgeo<=n;mgeo++)
              {
                igeo = indx(n,mgeo);
                for (iperiod = 1; iperiod <= GSM2GEO_DEG_DAY; iperiod++) 
                  induc_Em[igeo] += POMME_C_EM[igsm] * Em_eff * (qfac_day_r[n-1][iperiod - 1] * gsm_green_day_r[igsm][igeo][iperiod]
                                                                  - CONVENTION * qfac_day_i[n-1][iperiod - 1] * gsm_green_day_i[igsm][igeo][iperiod]);
                for (iperiod = 1; iperiod <= GSM2GEO_DEG_ANN; iperiod++) 
                  induc_Em[igeo] += POMME_C_EM[igsm] * Em_eff * (qfac_ann_r[n-1][iperiod - 1] * gsm_green_ann_r[igsm][igeo][iperiod]
                                                                  - CONVENTION * qfac_ann_i[n-1][iperiod - 1] * gsm_green_ann_i[igsm][igeo][iperiod]);

              } /* for mgeo */
          } /* for mgsm */
      for (igeo=0;igeo<MAX_S_COEFF;igeo++) geo_ind_field[igeo] += induc_Em[igeo];
    }

  /* psa - save geo_ext_field in GLOBAL_Q */
  GLOBAL_Q[0] = geo_ext_field[0] + fac*est*G10/MOM;
  GLOBAL_Q[1] = geo_ext_field[1] + fac*est*G11/MOM;
  GLOBAL_Q[2] = geo_ext_field[2] + fac*est*H11/MOM;

  costh = cos(theta);
  rad = rrel*6371.2;
  plmbar_(&p[0],&dp[0],&costh,&ext_deg,&inorm);
  magfdz_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&geo_ind_field[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&xi,&yi,&zi,&h,&f,&ainc,&d,&pot);
  magfdz_ext_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&geo_ext_field[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&xe,&ye,&ze,&h,&f,&ainc,&d,&pot);
  *x = xe + xi;
  *y = ye + yi;
  *z = ze + zi;
#if !EST_IST_via_SM2GEO
  *x += est_x + ist_x;
  *y += est_y + ist_y;
  *z += est_z + ist_z;
#endif

  if (TEST)
    {
      printf("fday,rrel,theta,phi,est,ist,imf_by,f107,Em:%f %f %f %f %f %f %f %f %f\n",fday,rrel,theta,phi,est,ist,imf_by,f107,Em);
      magfdz_ext_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&sm_in_geo[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&tx,&ty,&tz,&h,&f,&ainc,&d,&pot);
      printf("SM extern:  %10.4f %10.4f %10.4f\n",tx,ty,tz);
      magfdz_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&induc_sm[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&tx,&ty,&tz,&h,&f,&ainc,&d,&pot);
      printf("SM induced: %10.4f %10.4f %10.4f\n",tx,ty,tz);
#if EST_IST_via_SM2GEO
      magfdz_ext_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&est_in_geo[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&xe,&ye,&ze,&h,&f,&ainc,&d,&pot);
      printf("SM Est:     %10.4f %10.4f %10.4f\n",xe,ye,ze);
      magfdz_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&ist_in_geo[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&tx,&ty,&tz,&h,&f,&ainc,&d,&pot);
      printf("SM Ist:     %10.4f %10.4f %10.4f\n",tx,ty,tz);
      printf("SM Est+Ist: %10.4f %10.4f %10.4f\n",tx+xe,ty+ye,tz+ze);
#else
      printf("SM Est:     %10.4f %10.4f %10.4f\n",est_x,est_y,est_z);
      printf("SM Ist:     %10.4f %10.4f %10.4f\n",ist_x,ist_y,ist_z);
      printf("SM Est+Ist: %10.4f %10.4f %10.4f\n",est_x+ist_x,est_y+ist_y,est_z+ist_z);
#endif
      magfdz_ext_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&gsm_in_geo[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&tx,&ty,&tz,&h,&f,&ainc,&d,&pot);
      printf("GSM extern: %10.4f %10.4f %10.4f\n",tx,ty,tz);
      magfdz_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&induc_gsm[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&tx,&ty,&tz,&h,&f,&ainc,&d,&pot);
      printf("GSM ind:    %10.4f %10.4f %10.4f\n",tx,ty,tz);

      magfdz_ext_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&imfby_in_geo[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&tx,&ty,&tz,&h,&f,&ainc,&d,&pot);
      printf("GSM IMFBy:  %10.4f %10.4f %10.4f\n",tx,ty,tz);
      magfdz_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&induc_imfby[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&tx,&ty,&tz,&h,&f,&ainc,&d,&pot);
      printf("IMFBy ind:  %10.4f %10.4f %10.4f\n",tx,ty,tz);

      magfdz_ext_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&Em_in_geo[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&tx,&ty,&tz,&h,&f,&ainc,&d,&pot);
      printf("GSM Em:  %10.4f %10.4f %10.4f\n",tx,ty,tz);
      magfdz_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&induc_Em[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&tx,&ty,&tz,&h,&f,&ainc,&d,&pot);
      printf("Em ind:  %10.4f %10.4f %10.4f\n",tx,ty,tz);

      magfdz_ext_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&f107_in_geo[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&tx,&ty,&tz,&h,&f,&ainc,&d,&pot);
      printf("SM F107:  %10.4f %10.4f %10.4f\n",tx,ty,tz);
      magfdz_(&p[0],&dp[0],&theta,&phi,&rad,&ext_deg,&induc_f107[0],
          &dx[0],&dy[0],&dz[0],&dpot[0],&tx,&ty,&tz,&h,&f,&ainc,&d,&pot);
      printf("F107 ind:  %10.4f %10.4f %10.4f\n",tx,ty,tz);
      
      printf("\nComplete result: %10.4f %10.4f %10.4f\n",*x,*y,*z);

      printf("lon = %12.5f  <RETURN> ",phi/M_PI*180.0);
      getchar();
    } /* if TEST and right longitude */

  if (fabs(M_PI-theta) < MIN_POLE_DIST || fabs(theta) < MIN_POLE_DIST)
    {
      *x = 0;                /* not defined at poles*/
      *y = 0;                /* not defined at poles*/
    }
  return;
} /* pomme_external_geocentric */


void pomme_geodetic (int control[MAXCONTROL],
                     double fday, double geod_alt, double geod_lat, double lon,
                     double est, double ist, double imf_by, double f107, double Em,  /* time lagged but not peaked or subtracted mean */
                     double *geod_x, double *geod_y, double *geod_z, double *h, double *f, double *incl, double *decl)
{
  double theta, phi, lat_rad, delta, x,y,z, r, rrel;
  double mx, my, mz, cx, cy, cz, ex, ey, ez, fyear;

  geodetic2geocentric(geod_lat*M_PI/180.0, geod_alt, &lat_rad, &r); /* input and output latitudes in radian*/

  theta = M_PI/2.0 - lat_rad;
  phi = lon * M_PI / 180.0;
  fyear = 2000.0 + fday/365.25;
  pomme_internal_geocentric (control[POS_0], fyear, r, theta, phi, &mx, &my, &mz, &cx, &cy, &cz);

  rrel = r/6371.2;
  pomme_external_geocentric (control, fday, rrel, theta, phi, est, ist, imf_by, f107, Em, &ex, &ey, &ez);
  x = mx + cx + ex;
  y = my + cy + ey;
  z = mz + cz + ez;
  *f = sqrt(x*x + y*y + z*z);

  delta = (90.0-geod_lat) * M_PI / 180.0; /* delta is the geodetic co-latitude */
  geocentric2geodetic_vec(theta, delta, x, y, z, geod_x, geod_y, geod_z);

  *h     = sqrt(*geod_x * *geod_x + *geod_y * *geod_y);
  *incl  = 0.0;
  if (*f > 0.0) *incl = asin(*geod_z / *f);
  *incl *= 180.0 / M_PI;
  *decl  = atan2(*geod_y, *geod_x);
  *decl *= 180.0 / M_PI;
  return;
} /* pomme_geodetic */

void pomme_geocentric (int control[MAXCONTROL],
                     double fday, double r, double lat, double lon,
                     double est, double ist, double imf_by, double f107, double Em,  /* time lagged but not peaked or subtracted mean */
                     double *x, double *y, double *z)
{
  double theta, phi, lat_rad, rrel;
  double mx, my, mz, cx, cy, cz, ex, ey, ez, fyear;

  lat_rad = lat * M_PI / 180.0;
  theta = M_PI/2.0 - lat_rad;
  phi = lon * M_PI / 180.0;
  fyear = 2000.0 + fday/365.25;
  pomme_internal_geocentric (control[POS_0], fyear, r, theta, phi, &mx, &my, &mz, &cx, &cy, &cz);

  rrel = r/6371.2;
  pomme_external_geocentric (control, fday, rrel, theta, phi, est, ist, imf_by, f107, Em, &ex, &ey, &ez);
  *x = mx + cx + ex;
  *y = my + cy + ey;
  *z = mz + cz + ez;
  return;
} /* pomme_geocentric */
