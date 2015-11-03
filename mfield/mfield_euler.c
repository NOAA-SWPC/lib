
static size_t mfield_euler_bin(const size_t sat_idx, const double t,
                               const mfield_workspace *w);

/*
mfield_euler_idx()
  Return index of Euler angles for a given satellite and time period

Inputs: sat_idx - satellite index in [0,nsat-1]
        t       - time (CDF_EPOCH)
        w       - workspace

Return: index so that:
alpha = w->c[index]
beta = w->c[index + 1]
gamma = w->c[index + 2]

Notes:
1) The Euler angles are organized inside the coefficient vector as:

c = [ MF | SV | SA | Euler | Ext ]

and

Euler = [ Euler_0 | Euler_1 | ... | Euler_n ]

where Euler_k are the Euler angles for satellite k. Now

Euler_k = [ Euler_{k0} | Euler_{k1} | ... | Euler_{km} ]

and m = nbins - 1. The appropriate bin is selected according to
the time, as Euler angles are computed for regular time intervals,
specified by the euler_period parameter. Finally, for satellite k and
bin l, the Euler angles are:

Euler_{kl} = [ alpha beta gamma ]
*/

size_t
mfield_euler_idx(const size_t sat_idx, const double t, const mfield_workspace *w)
{
  /*size_t idx = w->euler_offset + 3 * sat_idx;*/
  size_t bin = mfield_euler_bin(sat_idx, t, w);
  size_t idx = w->euler_offset + w->offset_euler[sat_idx] + 3 * bin;

  assert(idx >= w->euler_offset && idx < w->euler_offset + w->neuler);

  return idx;
} /* mfield_euler_idx() */

int
mfield_euler_print(const char *filename, const size_t sat_idx,
                   const mfield_workspace *w)
{
  FILE *fp;
  size_t i;
  double t0 = w->data_workspace_p->t0[sat_idx];
  double t1 = w->data_workspace_p->t1[sat_idx];
  double t, dt;
  double alpha_mean = 0.0,
         beta_mean = 0.0,
         gamma_mean = 0.0;
  size_t nmean = 0;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "mfield_euler_print: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Euler angles for satellite %zu\n", sat_idx);
  fprintf(fp, "# Field %zu: timestamp (CDF_EPOCH)\n", i++);
  fprintf(fp, "# Field %zu: time (decimal years)\n", i++);
  fprintf(fp, "# Field %zu: alpha (degrees)\n", i++);
  fprintf(fp, "# Field %zu: beta (degrees)\n", i++);
  fprintf(fp, "# Field %zu: gamma (degrees)\n", i++);
  fprintf(fp, "# Field %zu: alpha - mean (arcseconds)\n", i++);
  fprintf(fp, "# Field %zu: beta - mean (arcseconds)\n", i++);
  fprintf(fp, "# Field %zu: gamma - mean (arcseconds)\n", i++);

  if (w->params.euler_period > 0.0)
    dt = w->params.euler_period * 86400000; /* convert to ms */
  else
    dt = t1 - t0;

  /* compute means of Euler angles */
  for (t = t0; t < t1; t += dt)
    {
      size_t euler_idx = mfield_euler_idx(sat_idx, t, w);
      double alpha = gsl_vector_get(w->c, euler_idx);
      double beta = gsl_vector_get(w->c, euler_idx + 1);
      double gamma = gsl_vector_get(w->c, euler_idx + 2);

      if (alpha == 0.0 || beta == 0.0 || gamma == 0.0)
        continue;

      alpha_mean += alpha;
      beta_mean += beta;
      gamma_mean += gamma;
      ++nmean;
    }

  if (nmean > 0)
    {
      alpha_mean /= (double) nmean;
      beta_mean /= (double) nmean;
      gamma_mean /= (double) nmean;
    }

  for (t = t0; t < t1; t += dt)
    {
      size_t euler_idx = mfield_euler_idx(sat_idx, t, w);
      double alpha = gsl_vector_get(w->c, euler_idx);
      double beta = gsl_vector_get(w->c, euler_idx + 1);
      double gamma = gsl_vector_get(w->c, euler_idx + 2);

      fprintf(fp, "%f %f %.10f %.10f %.10f %f %f %f\n",
              t,
              satdata_epoch2year(t),
              wrap180(alpha * 180.0 / M_PI),
              wrap180(beta * 180.0 / M_PI),
              wrap180(gamma * 180.0 / M_PI),
              (alpha - alpha_mean) * 180.0 / M_PI * 3600.0,
              (beta - beta_mean) * 180.0 / M_PI * 3600.0,
              (gamma - gamma_mean) * 180.0 / M_PI * 3600.0);
    }

  fclose(fp);

  return 0;
}

/*
mfield_euler_bin()
  Return bin number of Euler angles for a given satellite and time

Inputs: sat_idx - satellite index in [0,nsat-1]
        t       - time (CDF_EPOCH)
        w       - workspace

Return: bin number corresponding to Euler angle triplet
*/

static size_t
mfield_euler_bin(const size_t sat_idx, const double t, const mfield_workspace *w)
{
  double t0 = w->data_workspace_p->t0[sat_idx];
  double t1 = w->data_workspace_p->t1[sat_idx];
  double ratio = (t - t0) / (t1 - t0);
  size_t bin = (size_t) (w->nbins_euler[sat_idx] * ratio);

  if (bin >= w->nbins_euler[sat_idx])
    bin = w->nbins_euler[sat_idx] - 1;

  return bin;
}
