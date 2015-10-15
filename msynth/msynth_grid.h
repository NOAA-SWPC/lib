/*
 * msynth_grid.h
 */

#ifndef INCLUDED_msynth_grid_h
#define INCLUDED_msynth_grid_h

#ifndef INCLUDED_msynth_h
#include "msynth.h"
#define INCLUDED_msynth_h
#endif

typedef struct
{
  msynth_workspace *msynth_workspace_p;
} msynth_grid_workspace;

/*
 * Prototypes
 */

msynth_grid_workspace *msynth_grid_alloc(const size_t nphi, const msynth_workspace *msynth_p);
void msynth_grid_free(msynth_grid_workspace *w);

#endif /* INCLUDED_msynth_grid_h */
