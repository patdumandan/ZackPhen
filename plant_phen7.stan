data {
  int<lower=1> N;
  int<lower=1> Nplots;
  int<lower=1> Nyr;

  array[N] int<lower=1, upper=Nplots> plot_id;
  array[N] int<lower=1, upper=Nyr> year_id;

  array[N] int<lower=0> tot_F;
  array[N] int<lower=0> tot_NF;

  vector[N] DOYs;
  vector[Nyr] year_vec;
}

parameters {
  // Global intercept
  real alpha;

  // Year-level intercepts
  vector[Nyr] alpha_year_raw;
  real<lower=0> sigma_year;

  // Plot-level random intercepts
  vector[Nplots] u_plot_raw;
  real<lower=0> sigma_plot;

  // Peak timing
  real mu_bar;
  real beta_mu;
  real<lower=0> sigma_mu;
  vector[Nyr] mu_raw;

  vector[Nplots] u_plot_mu_raw;
  real<lower=0> sigma_mu_plot;

  //pre-and post peak widths
  real log_width_L_bar;
  real log_width_R_bar;
  real<lower=0> sigma_width_L;
  real<lower=0> sigma_width_R;
  vector[Nyr] width_L_raw;
  vector[Nyr] width_R_raw;

  // Plot-specific width sd
 vector[Nplots] width_L_plot_raw;
 vector[Nplots] width_R_plot_raw;
 real<lower=0> sigma_width_L_plot;
 real<lower=0> sigma_width_R_plot;

}

transformed parameters {

  vector[Nplots] u_plot= sigma_plot * (u_plot_raw - mean(u_plot_raw));

  vector[Nyr] alpha_year= alpha + sigma_year * alpha_year_raw;

  vector[Nyr] mu= mu_bar + beta_mu * year_vec + sigma_mu * mu_raw;

  vector[Nplots] u_plot_mu= sigma_mu_plot * u_plot_mu_raw;

  vector[Nyr] width_L_year= exp(log_width_L_bar + sigma_width_L * width_L_raw);

  vector[Nyr] width_R_year= exp(log_width_R_bar + sigma_width_R * width_R_raw);

  vector[Nplots] width_L_plot= exp(sigma_width_L_plot * width_L_plot_raw);

  vector[Nplots] width_R_plot = exp(sigma_width_R_plot * width_R_plot_raw);

  array[N] real eta;

  for (n in 1:N) {

  real d = DOYs[n] - (mu[year_id[n]] + u_plot_mu[plot_id[n]]);

  real width_L_np = width_L_year[year_id[n]] * width_L_plot[plot_id[n]];
  real width_R_np = width_R_year[year_id[n]] * width_R_plot[plot_id[n]];

  if (d < 0) {
    eta[n] = alpha_year[year_id[n]] + u_plot[plot_id[n]] - square(d) / (2 * square(width_L_np));
  } else {
    eta[n] = alpha_year[year_id[n]] + u_plot[plot_id[n]] - square(d) / (2 * square(width_R_np));
  }
}

model {
  // Intercepts
  alpha ~ normal(0, 2);
  alpha_year_raw ~ normal(0, 1);
  sigma_year ~ normal(0, 0.5);

  // Plot effects
  u_plot_raw ~ normal(0, 1);
  sigma_plot ~ normal(0, 0.5);

  // Peak timing
  mu_bar ~ normal(0, 1);
  beta_mu ~ normal(0, 0.5);
  mu_raw ~ normal(0, 1);
  sigma_mu ~ normal(0, 0.3);

  u_plot_mu_raw ~ normal(0, 1);
  sigma_mu_plot ~ normal(0, 0.3);

  // Width
  log_width_L_bar ~ normal(log(6), 0.25);
  log_width_R_bar ~ normal(log(6), 0.25);

  sigma_width_L ~ normal(0, 0.15);
  sigma_width_R ~ normal(0, 0.15);

  width_L_raw ~ normal(0, 1);
  width_R_raw ~ normal(0, 1);

  width_L_plot_raw ~ normal(0, 1);
  width_R_plot_raw ~ normal(0, 1);

  sigma_width_L_plot ~ normal(0, 0.1);
  sigma_width_R_plot ~ normal(0, 0.1);


  // Likelihood
  for (n in 1:N)
    tot_F[n] ~ binomial_logit(tot_F[n] + tot_NF[n], eta[n]);
}

generated quantities {
  array[N] int y_pred;
  for (n in 1:N)
    y_pred[n] =
      binomial_rng(tot_F[n] + tot_NF[n], inv_logit(eta[n]));
}
