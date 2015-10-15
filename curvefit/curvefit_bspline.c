/*
 * curvefit_bspline.c
 *
 * Routines for B-spline fitting of data
 */

#include <gsl/gsl_bspline.h>

typedef struct
  {
    size_t n;
    size_t p;
    size_t nbreak;
    gsl_bspline_workspace *bspline_p;
  }
bspline_state_t;

static int
bspline_alloc(void *vstate, const size_t n, const size_t p)
{
  bspline_state_t *state = (bspline_state_t *) vstate;
  const size_t k = 4; /* cubic B-splines */

  state->p = p;
  state->n = n;
  state->nbreak = p + 2 - k;

  state->bspline_p = gsl_bspline_alloc(k, state->nbreak);

  return GSL_SUCCESS;
} /* bspline_alloc() */

static void
bspline_free(void *vstate)
{
  bspline_state_t *state = (bspline_state_t *) vstate;
  
  if (state->bspline_p)
    gsl_bspline_free(state->bspline_p);
} /* bspline_free() */

static int
bspline_init(void *vstate, const double *x, const double *y)
{
  bspline_state_t *state = (bspline_state_t *) vstate;
  int s;
  const size_t n = state->n;
  double xmin, xmax;
  gsl_vector_const_view v = gsl_vector_const_view_array(x, n);

  gsl_vector_minmax(&v.vector, &xmin, &xmax);
  s = gsl_bspline_knots_uniform(xmin, xmax, state->bspline_p);

  return s;
} /* bspline_init() */

/*
bspline_row()
  Construct a row of the least squares design matrix
for bspline model

Inputs: vstate - state
        x      - data value
        v      - (output) where to store matrix row

Return: success or error
*/

static int
bspline_row(void *vstate, const double x, gsl_vector *v)
{
  bspline_state_t *state = (bspline_state_t *) vstate;

  /* fill in row i of A: A(i,j) = B_j(xi) */
  gsl_bspline_eval(x, v, state->bspline_p);

  return GSL_SUCCESS;
} /* bspline_row() */

static const curvefit_type bspline_type =
{
  "bspline",
  sizeof(bspline_state_t),
  &bspline_alloc,
  &bspline_init,
  &bspline_row,
  &bspline_free
};

const curvefit_type *curvefit_bspline = &bspline_type;
