#define R_E               6371.2  /* Magnetic reference radius  */
#define MAX_0_DEG         720     /* Max degree of static field */
#define MAX_0_NCOEFF      (MAX_0_DEG * (MAX_0_DEG + 2))
#define MAX_V_DEG         20      /* Max degree of time derivatives (all same) */
#define MAX_V_NCOEFF      (MAX_V_DEG * (MAX_V_DEG + 2))

#define MAX_EXT_DEG       3       /* Max degree of external field, limited to 3 */
#define MAX_EXT_NCOEFF    (MAX_EXT_DEG * (MAX_EXT_DEG + 2))  

#define POS_0             0       /* static internal */
#define POS_1             1       /* SV */
#define POS_SM            4       /* static ring current */
#define POS_SM_IND        5       /* if asymmetric then could induce field, not implemented */
#define POS_EST           6       /* external RC field */
#define POS_IST           7       /* induced RC field */
#define POS_GSM           8       /* external GSM static field */
#define POS_GSM_IND       9       /* induced internal*/
#define POS_IMFBY         10      /* external IMFBY correlated field */
#define POS_IMFBY_IND     11      /* induced internal*/
#define POS_EST_FAC       12      /* factor for Est/Ist field */
#define POS_F107          13      /* external F107 correlated field in SM */
#define POS_F107_IND      14      /* induced internal*/
#define POS_EM            15      /* external Em correlated field in GSM */
#define POS_EM_IND        16      /* induced internal*/
#define MAXCONTROL        17      /* number of control parameters */

#define MAX_S_COEFF       15       /* spatial coefficients, must change also in gsm2geo and sm2geo  */
#define MAX_T_COEFF       13       /* temporal coefficients */
#define MAXPERIOD         7        /* const, 24h, 12h, 8h ..., must change also in gsm2geo and sm2geo*/

typedef int tcontrol[MAXCONTROL];
