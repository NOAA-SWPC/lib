/*
 * estist_calc.h
 */

#ifndef INCLUDED_estist_calc_h
#define INCLUDED_estist_calc_h

#include <indices/indices.h>

#define ESTIST_NUM_DAYS      365

/* 1-hr samples */
#define ESTIST_NDST          (ESTIST_NUM_DAYS * 24)

typedef struct
{
  double dst[ESTIST_NDST];
  double est[ESTIST_NDST];
  double ist[ESTIST_NDST];
  int model;

  dst_workspace *dst_workspace_p;
} estist_calc_workspace;

/*
 * Prototypes
 */

void weidelt_dst_(int *model, int *ndst, double *dstvec, double *est,
                  double *ist);

estist_calc_workspace *estist_calc_alloc(void);
void estist_calc_free(estist_calc_workspace *w);
int estist_calc(int ndst, double dst[], double est[], double ist[],
                estist_calc_workspace *w);
int estist_calc_get(const time_t t, double *est, double *ist,
                    estist_calc_workspace *w);

#endif /* INCLUDED_estist_calc_h */
