data {
  int<lower=1> N;                  // number of observations
  int<lower=0> tot_flwr[N];       // successes
  int<lower=0> tot_NF[N];         // failures
  int<lower=1> J;                 // number of plots
  int<lower=1, upper=J> plot[N];  // plot index

  vector[N] DOYs;
  vector[N] DOYsqs;
  vector[N] yearc;
}

transformed data {
  //needed to do cbind(tot_NF, tot_flwr)
  
  int<lower=0> trials[N]; 
  
  for (i in 1:N)
    trials[i] = tot_flwr[i] + tot_NF[i];  // total trials per obs
}

parameters {
  
  real beta_DOYs;
  real beta_DOYsqs; //need to constrain so effect doesn't balloon)
  real beta_yearc;
  real beta_DOYs_yearc;
  real beta_DOYsqs_yearc;

  real alpha; // intercept
  real<lower=0> sigma_plot;         // std dev of random intercepts
  vector[J] z_plot;                 // non-centered random effects
  
}

transformed parameters {
  
  vector[J] alpha_plot;
  alpha_plot = alpha + sigma_plot * z_plot;
}

model {
  vector[N] eta;

  // Priors
  beta_DOYs ~ normal(0, 5);
  beta_DOYsqs ~ normal(0, 5);
  beta_yearc ~ normal(0, 5);
  beta_DOYs_yearc ~ normal(0, 5);
  beta_DOYsqs_yearc ~ normal(0, 5);
  alpha ~ normal(0, 5);
  z_plot ~ normal(0, 1);
  sigma_plot ~ exponential(1);

  // Linear predictor
  for (i in 1:N) {
    eta[i] = alpha_plot[plot[i]]
           + beta_DOYs * DOYs[i]
           + beta_DOYsqs * DOYsqs[i]
           + beta_yearc * yearc[i]
           + beta_DOYs_yearc * DOYs[i] * yearc[i]
           + beta_DOYsqs_yearc * DOYsqs[i] * yearc[i];
  }

  // Likelihood
  tot_flwr ~ binomial_logit(trials, eta);
}
