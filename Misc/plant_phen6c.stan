data {
  int<lower=1> N;          // observations
  int<lower=1> Nplots;     // number of plots
  int<lower=1> Nyr;        // number of years

  array[N] int<lower=1, upper=Nplots> plot_id;
  array[N] int<lower=1, upper=Nyr> year_id;

  array[N] int<lower=0> tot_F;
  array[N] int<lower=0> tot_NF;

  vector[N] DOYs;
  vector[Nyr] year_vec;

  real DOY_mean;
  real DOY_sd;
}

parameters {
  real alpha;

  vector[Nyr] alpha_year_raw;
  real<lower=0> sigma_year;

  vector[Nplots] u_plot_raw;
  real<lower=0> sigma_plot;

  real mu_bar;                // mean of mu/peak time
  real beta_mu;               // slope on year
  vector[Nyr] mu_raw;         // non-centered year deviations
  real<lower=0> sigma_mu;     // SD among years

  vector[Nplots] u_plot_mu_raw;   // plot deviation in mu
  real<lower=0> sigma_mu_plot;

  real width_bar;
  vector[Nyr] width_raw;
  real<lower=0> sigma_width;
}

transformed parameters {
  vector[Nyr] alpha_year = alpha + sigma_year * alpha_year_raw;

  vector[Nplots] u_plot = sigma_plot * u_plot_raw;

  vector[Nyr] mu_year= mu_bar + beta_mu * year_vec+ sigma_mu * mu_raw;

  vector[Nplots] u_plot_mu = sigma_mu_plot * u_plot_mu_raw;

  vector[Nyr] width_year= exp(width_bar + sigma_width * width_raw);//exp() to constrain to +

  vector[N] eta;
  for (n in 1:N) {
    real mu_full = mu_year[year_id[n]] + u_plot_mu[plot_id[n]]; //clean

    eta[n] = alpha_year[year_id[n]] +u_plot[plot_id[n]]-square(DOYs[n] - mu_full)/square(width_year[year_id[n]]);
  }
}

model {
  alpha ~ normal(0, 2);

  alpha_year_raw ~ normal(0, 1);
  sigma_year ~ normal(0, 1);

  u_plot_raw ~ normal(0, 1);
  sigma_plot ~ normal(0, 1);

  mu_bar ~ normal(0, 1);
  beta_mu ~ normal(0, 1);
  mu_raw ~ normal(0, 1);
  sigma_mu ~ normal(0, 1);

  u_plot_mu_raw ~ normal(0, 1);
  sigma_mu_plot ~ normal(0, 1);

  width_bar ~ normal(0, 1);
  width_raw ~ normal(0, 1);
  sigma_width ~ normal(0, 1);

  for (n in 1:N)
    tot_F[n] ~ binomial_logit(tot_F[n] + tot_NF[n], eta[n]);
}

generated quantities {
  array[N] int y_pred;

  for (n in 1:N)
    y_pred[n] = binomial_rng(tot_F[n] + tot_NF[n],inv_logit(eta[n]));
}
