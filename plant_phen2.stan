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

  // ---- Hierarchical slopeS, IS THIS RIGHT?? do we just want independent slopes for DOYsqs but no partial pooling?
  real mu_DOY; //mean slope for DOY
  real mu_DOYsq; //mean slope for DOYsqs

  real<lower=0> sigma_DOY;//sd of DOY
  real<lower=0> sigma_DOYsq;

  vector[Nyr] beta_DOYs_raw; //year-specific slope values
  vector[Nyr] beta_DOYsqs_raw;
}

transformed parameters {
  vector[Nplots] u_plot = u_plot_raw*sigma_plot;

  vector[Nyr] alpha_year = alpha+sigma_year*alpha_year_raw;

  vector[Nyr] beta_DOYs    = mu_DOY+sigma_DOY*beta_DOYs_raw;

  vector[Nyr] beta_DOYsqs_uncon = mu_DOYsq + sigma_DOYsq * beta_DOYsqs_raw;
  vector<upper=0>[Nyr] beta_DOYsqs; //constrain to bwe negative

  for (y in 1:Nyr)
    beta_DOYsqs[y] = -abs(beta_DOYsqs_uncon[y]); // or beta_DOYsqs[y] = fmin(beta_DOYsqs_uncon[y], 0);

  array[N] real eta;

  for (n in 1:N) {
    eta[n] =alpha_year[year_id[n]]+
            beta_DOYs[year_id[n]]*DOYs[n]+
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

  // ---- Hierarchical slope priors ----
  mu_DOY   ~ normal(0, 5); //or 0,2
  mu_DOYsq ~ normal(-1, 1);      // favors concave-down shape

  sigma_DOY   ~ normal(0, 1); //or 0,2
  sigma_DOYsq ~ normal(0, 1); //or 0,2

  beta_DOYs_raw ~ normal(0, 1);
  beta_DOYsqs_raw ~ normal(0, 1);

  // Likelihood
  for (i in 1:N)
    tot_F[i] ~ binomial_logit(tot_F[i] + tot_NF[i], eta[i]);
}

generated quantities {
  array[N] int y_pred;

  for (n in 1:N)
    y_pred[n] = binomial_rng(tot_F[n] + tot_NF[n], inv_logit(eta[n]));
}

// //warnings for this code:
// Warning: 4 of 8000 (0.0%) transitions ended with a divergence.
// See https://mc-stan.org/misc/warnings for details.
//
// Warning: 3604 of 8000 (45.0%) transitions hit the maximum treedepth limit of 10.
// See https://mc-stan.org/misc/warnings for details.
