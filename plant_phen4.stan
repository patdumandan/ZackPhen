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

  vector[Nyr] alpha_year_z; // Year-level intercepts (NC)
  real<lower=0> sigma_year;

  vector[Nplots] u_plot_z; // Plot-level random intercepts, also for NCP
  real<lower=0> sigma_plot;

  real mu_bar;    // mean predicted peak
  real<lower=0> sigma_mu; // sd of yearly peak times across years
  vector[Nyr] mu_z; //year-specific slope values (NC)

//to fix divergence issues, maybe reparameterize the quadratic as a gaussian curve with mean and sd params,
//and add hierarchical structure

  //vector[Nyr] beta_DOYsqs;  // year specific curvature (second order effect)
  real mu_beta_DOYsq;
  real <lower=0> sigma_beta_DOYsqs;
  vector [Nyr] sigma_beta_DOYsqs_z;
}

transformed parameters {
  vector[Nplots] u_plot = u_plot_z*sigma_plot;

  vector[Nyr] alpha_year = alpha+sigma_year*alpha_year_z;

  vector[Nyr] mu_year = mu_bar+sigma_mu*mu_z;

  vector[Nyr] beta_DOYsqs_bar=exp(mu_beta_DOYsq+sigma_beta_DOYsqs*sigma_beta_DOYsqs_z); //curvature

 // vector[Nyr] beta_DOYs    = -2*(mu.*beta_DOYsqs);

//  vector[Nyr] beta_DOYsqs_uncon = mu_DOYsq + sigma_DOYsq * beta_DOYsqs_raw;
//  vector<upper=0>[Nyr] beta_DOYsqs; //constrain to bwe negative

  //for (y in 1:Nyr)
    //beta_DOYsqs[y] = -abs(beta_DOYsqs_uncon[y]); // or beta_DOYsqs[y] = fmin(beta_DOYsqs_uncon[y], 0);

  array[N] real eta;

  for (n in 1:N) {
    int y = year_id[n];

    eta[n] =alpha_year[y]- square(DOYs[n] - mu_year[y]) / (2 * square(beta_DOYsqs_bar[y]))+ u_plot[plot_id[n]];
}
}

model {
  // Intercepts
  alpha ~ normal(0, 5);
  alpha_year_z ~ normal(0, 1);
  sigma_year ~ normal(0, 2);

  // Plot effects
  u_plot_z ~ normal(0, 1);
  sigma_plot ~ normal(0, 2);

  //beta_DOYsqs ~ normal(-1, 1);      // negative to force concave down

  // peak time -- Hierarchical priors ----
  mu_bar   ~ normal(0, 5);
  sigma_mu ~ normal(0,1);
  mu_z ~ normal(mu_bar, sigma_mu);

  mu_beta_DOYsq ~ normal(0, 1);
  sigma_beta_DOYsqs ~ normal (0,1);
  sigma_beta_DOYsqs_z ~ normal(0, 1);

  // Likelihood
  for (i in 1:N)
    tot_F[i] ~ binomial_logit(tot_F[i] + tot_NF[i], eta[i]);
}

generated quantities {
  array[N] int y_pred;

  for (n in 1:N)
    y_pred[n] = binomial_rng(tot_F[n] + tot_NF[n], inv_logit(eta[n]));
}

// Warning: 16 of 8000 (0.0%) transitions ended with a divergence.
// See https://mc-stan.org/misc/warnings for details.
//
// Warning: 1238 of 8000 (15.0%) transitions hit the maximum treedepth limit of 10.
// See https://mc-stan.org/misc/warnings for details.
