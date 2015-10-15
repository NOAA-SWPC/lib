const size_t powell_N = 2;
const size_t powell_P = 2;

double powell_X[2][2] = {
};

double powell_F[2] = {
};

double powell_x0[4] = { 3.0, 1.0 };

int
powell_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);

  gsl_vector_set(f, 0, x0);
  gsl_vector_set(f, 1, 10.0*x0/(x0 + 0.1) + 2.0*x1*x1);

  return GSL_SUCCESS;
}

int
powell_df (const gsl_vector * x, void *params, gsl_matrix * df)
{
  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);
  double term = x0 + 0.1;

  gsl_matrix_set(df, 0, 0, 1.0);
  gsl_matrix_set(df, 0, 1, 0.0);
  gsl_matrix_set(df, 1, 0, 1.0 / (term * term));
  gsl_matrix_set(df, 1, 1, 4.0 * x1);

  return GSL_SUCCESS;
}

int
powell_fdf (const gsl_vector * x, void *params,
           gsl_vector * f, gsl_matrix * df)
{
  powell_f (x, params, f);
  powell_df (x, params, df);

  return GSL_SUCCESS;
}
