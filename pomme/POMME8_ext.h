#define DEG_IMFBY  1
#define DEG_EM     1
#define DEG_F107   2
#define DEG_SM     2
#define DEG_GSM    2
#define CONVENTION 1.0     /* 1.0 for exp(iwt), -1.0 for exp(-iwt) */
#define EM_MEAN    0.5     /* is subtracted to get variation around mean*/
#define EM_PEAK    8       /* Em approaches this peak assymptotically*/
#define F107_MEAN  120.0   /* is subtracted to get variation around mean*/
#define EST_FAC    0.8655
static double POMME_C_F107[8] = {
   0.0833, 0.0000, 0.0000,
   0.0000, -0.0018, 0.0010, 0.0000, 0.0000};
static double POMME_C_EM[3] = {
   1.4453, 0.0000, 0.0000};
static double POMME_C_IMFBY[3] = {
   0.0000, 0.0688, -0.2657};
static double POMME_C_SM[8] = {
   5.7304, 0.0000, 0.0000,
   0.0000, -0.9266, 1.1032, 0.0000, 0.0000};
static double POMME_C_GSM[8] = {
   8.2474, 0.1169, 0.5475, 
   0.1253, -0.6716, -1.4288, -0.0140, 0.5469};
