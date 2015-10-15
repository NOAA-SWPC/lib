/*
 * estist2.h
 * Patrick Alken
 */

#ifndef INCLUDED_estist2_h
#define INCLUDED_estist2_h

#include <time.h>

/* Est/Ist parameters */

#define ESTIST2_MAX_DATA      20000
#define ESTIST2_MAX_BUFFER    8192
#define ESTIST2_MAX_YEAR      500

/*
 * first year to store data (data is stored sequentially for
 * [ESTIST2_FIRST_YEAR, ESTIST2_FIRST_YEAR + ESTIST2_MAX_YEAR]
 */
#define ESTIST2_FIRST_YEAR    1900

/* missing data */
#define ESTIST2_BAD_VALUE     (9999.99)

/*
 * this struct keeps track of which year is stored where in the workspace
 * in case of missing data files etc
 */
typedef struct
{
  /*
   * data for a whole year
   * first index: month
   * second index: day of month
   * third index: ut hour (0-23)
   */
  double E_st[12][31][24];
  double I_st[12][31][24];
} estist2_data;

typedef struct
{
  /* each index represents a year of data starting at year 1900 */
  estist2_data estist[ESTIST2_MAX_YEAR];
} estist2_workspace;

#define ESTIST2_IDX_FILE       "/nfs/satmag_work/palken/data/Est_Ist_index.pli"

/*
 * Prototypes
 */

estist2_workspace *estist2_alloc(const char *datadir);
void estist2_free(estist2_workspace *w);
int estist2_get(time_t t, double *E_st, double *I_st, estist2_workspace *w);

#endif /* INCLUDED_estist2_h */
