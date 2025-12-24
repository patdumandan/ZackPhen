//Plant only phenology model
data {
  int<lower=1> ND; //DOY pred grid, can be removed
  int<lower=1> N; //nrow/nobs
  int<lower=1> Nplots; //no. opf plots
  int<lower=1> Nyr; //no. of obs years
  real DOY_sd;
  real DOY_mean;

  array[N] int<lower=1, upper=Nplots> plot_id;
  array[N] int<lower=1, upper=Nyr> year_id;

  array[N] int<lower=0> tot_F;
  array[N] int<lower=0> tot_NF;

  vector[N] DOYs;
  vector[N] DOYsqs;
}

parameters {
  real alpha; // Global intercept

  vector[Nyr] alpha_year_raw; // Year-level intercepts
  real<lower=0> sigma_year;

  vector[Nplots] u_plot_raw; // Plot-level random intercepts, also for NCP
  real<lower=0> sigma_plot;

  real mu_bar;    // mean peak time across years
  real<lower=0> sigma_mu; // sd of yearly peak times across years
  vector[Nyr] mu_raw;

  vector[Nyr] beta_DOYsqs;  // year specific curvature (second order effect)
}

transformed parameters {
  vector[Nplots] u_plot = u_plot_raw*sigma_plot;

  vector[Nyr] alpha_year = alpha+sigma_year*alpha_year_raw;

  vector[Nyr] mu = mu_bar + sigma_mu*mu_raw; //year-specific peak times

  vector[Nyr] beta_DOYs = -2*(mu.*beta_DOYsqs);

  array[N] real eta;

  for (n in 1:N) {
    eta[n] =alpha_year[year_id[n]]+ //is this still correct with the way beta_DOYs is estimated? why?how?
            beta_DOYs[year_id[n]]*DOYs[n]+ //identifiability issues?
            beta_DOYsqs[year_id[n]]*DOYsqs[n] +
            u_plot[plot_id[n]];}
}

model {
  // Intercepts
  alpha ~ normal(0, 5);
  alpha_year_raw ~ normal(0, 1);
  sigma_year ~ normal(0, 2);

  // Plot effects
  u_plot_raw ~ normal(0, 1);
  sigma_plot ~ normal(0, 2);

  // Squared term
  beta_DOYsqs ~ normal(-1, 1);      // favors concave-down shape

  // peak time -- Hierarchical priors ----
  mu_bar   ~ normal(0, 2); //or 0,2
  //sigma_mu ~ normal(0,1);
  sigma_mu ~ student_t(4, 0, 0.2);
  mu_raw ~ normal(0, 1);

  // Likelihood
  for (i in 1:N)
    tot_F[i] ~ binomial_logit(tot_F[i] + tot_NF[i], eta[i]);
}

generated quantities {
  array[N] int y_pred;

  for (n in 1:N)
    y_pred[n] = binomial_rng(tot_F[n] + tot_NF[n], inv_logit(eta[n]));
}

// Warning: 4747 of 8000 (59.0%) transitions ended with a divergence.
// See https://mc-stan.org/misc/warnings for details.
//
// Warning: 3027 of 8000 (38.0%) transitions hit the maximum treedepth limit of 10.
// See https://mc-stan.org/misc/warnings for details.
