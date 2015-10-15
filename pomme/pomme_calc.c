/* -------------------------------------------------------------- */
/* Simple main program for POMME7 model computation               */
/*                                                                */
/* The program asks for time, location and magnetic indices.      */
/* In case the magnetic indices are unavailable, use their        */
/* default value.                                                 */
/* Contributions to the model can be turned on and off using      */
/* control parameters. These can also be used to set the          */
/* desired degree of the internal and external field.             */
/*                                                                */
/*                                          Stefan Maus, Feb-2012 */
/* -------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>

#include "./pomme_mod.h"
#include "./POMME8_ext.h"


static tcontrol control;        /* array of user-control parameters */

extern void pomme_geodetic (int control[MAXCONTROL],
                            double fday, double geod_alt, double geod_lat, double lon,
                            double est, double ist, double imf_by, double f107, double Em,
                            double *geod_x, double *geod_y, double *geod_z, 
                            double *h, double *f, double *incl, double *decl);

static  double lon, geod_lat, geod_alt, fday;                      /* Position and time */
static  double E_st=0, I_st=0, IMF_By=0, f107=120.0, Em=0.5;       /* Magnetic indices for time varying external field */
static  double x, y, z,  h, f, incl, decl;                         /* Resulting field */

int main(void)
{
  int i, model_nbytes;

 /* mainfield_open_file_mp initializes the control parameters. You can now modify these: */

  control[POS_0] =         133;   /* default: model.deg_int, 0=off */
  control[POS_1] =         15;    /* default: model.deg_1st, 0=off */

  control[POS_SM] =        DEG_SM;       /* degree of SM field, default = (*model).deg_sm, 0=off */
  control[POS_SM_IND] =    DEG_SM;       /* degree of field induced by SM field, default = (*model).deg_sm, 0=off */
  control[POS_F107] =      DEG_F107;     /* degree of SM field, default = (*model).deg_sm, 0=off */
  control[POS_F107_IND] =  DEG_F107;     /* degree of field induced by SM field, default = (*model).deg_sm, 0=off */
  control[POS_EST] =       1;                  /* Est dependent field: 1=on, 0=off */
  control[POS_IST] =       1;                  /* Ist dependent field: 1=on, 0=off */
  control[POS_EST_FAC] =   1;                  /* 1 = use model Est factor, 0 = use factor 1.0 */

  control[POS_GSM] =       DEG_GSM;     /* degree of GSM field, default = (*model).deg_gsm, 0=off */
  control[POS_GSM_IND] =   DEG_GSM;     /* degree of field induced by GSM field, default = (*model).deg_gsm, 0=off */
  control[POS_IMFBY] =     1;                 /* IMF-By dependent field: 1=on, 0=off */
  control[POS_IMFBY_IND] = 1;                 /* Induction by IMF-By dependent field: 1=on, 0=off */
  control[POS_EM] =        DEG_EM;      /* degree of GSM field, default = (*model).deg_gsm, 0=off */
  control[POS_EM_IND] =    DEG_EM;      /* degree of field induced by GSM field, default = (*model).deg_gsm, 0=off */


  printf("Universal time in floating point days since 00:00:00 on 1/1/2000? ");
  scanf("%lg%*[^\n]",&fday);
  printf("Longitude (degrees, positive to East)? ");
  scanf("%lg%*[^\n]",&lon);
  printf("Geodetic Latitude (degrees, positive to North)? ");
  scanf("%lg%*[^\n]",&geod_lat);
  printf("Altitude above mean sea level in km (actually WGS84 ellipsoid)? ");
  scanf("%lg%*[^\n]",&geod_alt);

  printf("IMF-By 35 minutes earlier (enter '0' if unknown)? ");
  scanf("%lg%*[^\n]",&IMF_By);
  printf("Merging electric field, 60 minutes earlier (enter '0.5' if unknown)? ");
  scanf("%lg%*[^\n]",&Em);
  printf("E_st (enter '0' if unknown)? ");
  scanf("%lg%*[^\n]",&E_st);
  printf("I_st (enter '0' if unknown)? ");
  scanf("%lg%*[^\n]",&I_st);
  printf("F10.7, 81-day average, 20 months earlier (enter '120' if unknown)? ");
  scanf("%lg%*[^\n]",&f107);

  
  pomme_geodetic(control, fday, geod_alt, geod_lat, lon, E_st, I_st,
                 IMF_By, f107, Em, &x, &y, &z, &h, &f, &incl, &decl);                                         
  
  printf("\n\n");
  printf("    fday     lon     lat     alt    E_st    I_st   IMF_By     Em   F10.7-av\n");
  printf("----------------------------------------------------------------------------\n");
  printf("%8.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %9.1f\n\n",fday,lon,geod_lat,geod_alt,E_st,I_st,IMF_By,Em,f107);

  printf("        X          Y          Z          H          F       Incl       Decl\n");
  printf("----------------------------------------------------------------------------\n");
  printf("%10.2f %10.2f %10.2f %10.2f %10.2f %10.5f %10.5f\n\n",x,y,z,h,f,incl,decl);
  
} /* End. */
