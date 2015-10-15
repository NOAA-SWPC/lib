typedef struct
  {
    size_t n;
    size_t p;
  }
poly_state_t;

static int
poly_alloc (void *vstate, const size_t n, const size_t p)
{
  poly_state_t *state = (poly_state_t *) vstate;

  state->n = n;
  state->p = p;

  return GSL_SUCCESS;
} /* poly_alloc() */

static void
poly_free (void *vstate)
{
} /* poly_free() */

static int
poly_init (void *vstate, const double *x, const double *y)
{
  return GSL_SUCCESS;
}

/*
poly_row()
  Construct a row of the least squares design matrix
for polynomial model

Inputs: x - data value
        v - (output) where to store matrix row
            v(j) = x^j

Return: success or error
*/

static int
poly_row(void *vstate, const double x, gsl_vector *v)
{
  size_t p = v->size;
  size_t j;
  double term = 1.0;

  for (j = 0; j < p; ++j)
    {
      gsl_vector_set(v, j, term);
      term *= x;
    }

  return GSL_SUCCESS;
} /* poly_row() */

static const curvefit_type poly_type =
{
  "poly",
  sizeof(poly_state_t),
  &poly_alloc,
  &poly_init,
  &poly_row,
  &poly_free
};

const curvefit_type *curvefit_poly = &poly_type;
