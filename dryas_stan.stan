data {
  int<lower=1> N; //obs
  int<lower=1> J; //plots
  int<lower=1> Y;//yrs category
  vector[N] yrs;

  array[N] int<lower=1, upper=J> plot_id;
  array[N] int<lower=1, upper=Y> year_id;

  array[N] int<lower=0> tot_flwr;
  array[N] int<lower=0> tot_trials;

  vector[N] DOYs;
  vector[N] DOYsqs;

  vector[N] raw_DOYs;

}

parameters {
  real beta0; //intercept
  real beta_yr; //year effect

  real beta_DOYs_global;
  vector[Y] z_beta_DOYs; //year-specific effect

  real beta_DOYsqs_global;
  vector[Y] z_beta_DOYsqs; //year-specific effect

  real<lower=0> sigma_plot;
  vector[J] z_plot; //plot-specific effect
}

transformed parameters {
  vector[J] plot_intercept;

  plot_intercept = sigma_plot * z_plot;

}

model {
  vector[N] eta;

  // Priors
  beta0 ~ normal(0, 5);

  beta_DOYs_global ~ normal(0, 5);
  z_beta_DOYs ~ normal(0, 1);

  beta_DOYsqs_global ~ normal(0, 1); // encourage negative curvature
  z_beta_DOYsqs ~ normal(0, 1);

  z_plot ~ normal(0, 1);
  sigma_plot ~ exponential(1);

  // Linear predictor
  for (n in 1:N) {
    int y = year_id[n];
    eta[n] = beta0 + beta_DOYs_global+ beta_DOYsqs_global+
             beta_yr+
             z_beta_DOYs[year_id[n]] * yrs[n] +
             z_beta_DOYsqs[year_id[n]] * yrs[n] +
             plot_intercept[plot_id[n]];
  }

  // Likelihood
  tot_flwr ~ binomial_logit(tot_trials, eta);

  // Soft constraint: peak_doy ∈ [120, 260]
  // for (y in 1:Y) {
  //   real peak_std = -beta_DOYs[y] / (2 * beta_DOYsqs[y]);
  //   real peak_doy_y = peak_std * doy_sd + doy_mean;
  //
  //   // Encourage peak within 120–260:
  //   target += normal_lpdf(peak_doy_y | 190, 30);  // center at midpoint 190, wider SD
  // }
}

generated quantities {
  array[N] int y_pred;
  vector[N] eta;
  vector[Y] peak_std_y;
  vector[Y] peak_doy;

  for (n in 1:N) {
    int y = year_id[n];
   eta[n] = beta0 + beta_DOYs_global+ beta_DOYsqs_global+
             z_beta_DOYs[year_id[n]] * yrs[n] +
             z_beta_DOYsqs[year_id[n]] * yrs[n] +
             plot_intercept[plot_id[n]];
    y_pred[n] = binomial_rng(tot_trials[n], inv_logit(eta[n]));
  }

  for (y in 1:Y) {
    peak_std_y[y] = -z_beta_DOYs[y] / (2 * z_beta_DOYsqs[y]);
    peak_doy[y] = peak_std_y[y] * sd(raw_DOYs) + mean(raw_DOYs);

  //   // Clamp peak_doy to [120, 260] if desired in output:
  //   peak_doy[y] = fmin(260, fmax(120, peak_doy[y]));
   }
}
