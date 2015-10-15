/*
 * ut1.h
 */

#ifndef INCLUDED_ut1_h
#define INCLUDED_ut1_h

#define UT1_DATA_FILE    "/nfs/satmag_work/palken/data/DUT1/dut1.dat"

#define UT1_BAD_VALUE    (-9999.9)

#define UT1_MAX_DATA     50000
#define UT1_MAX_BUFFER   8192
#define UT1_MAX_YEAR     500

/*
 * first year to store data (data is stored sequentially for
 * [UT1_FIRST_YEAR, UT1_FIRST_YEAR + UT1_MAX_YEAR]
 */
#define UT1_FIRST_YEAR   2000

typedef struct
{
  double data[12][31]; /* UT1-UTC in s */
} ut1_data;

typedef struct
{
  ut1_data ut1[UT1_MAX_YEAR];
} ut1_workspace;

/*
 * Prototypes
 */

ut1_workspace *ut1_alloc(const char *filename);
void ut1_free(ut1_workspace *w);
int ut1_get(const time_t t, double *result, ut1_workspace *w);

#endif /* INCLUDED_ut1_h */
