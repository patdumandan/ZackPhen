data {
  int<lower=1> N;                    // number of observations
  int<lower=1> Nplots;              // number of plots
  int<lower=1> Nyr;                 // number of years
  real DOY_sd;                      // sd of raw DOY
  real DOY_mean;                    // mean of raw DOY

  array[N] int<lower=1, upper=Nplots> plot_id;
  array[N] int<lower=1, upper=Nyr> year_id;

  array[N] int<lower=0> y;               // total abundance
  array[N] real<lower=0> obs_days;       // trapping effort (offset)

  vector[N] DOYs;                        // standardized DOY
  vector[N] DOYsqs;                      // DOY^2 (quadratic term)
}

parameters {
  real alpha;                            // global intercept
  real<lower=0> phi;                     // dispersion parameter

  vector[Nyr] alpha_year_raw; // Year-level intercepts
  real<lower=0> sigma_year;

  vector[Nplots] u_plot_raw;             // NCP for plot (small samples for some plots)
  real <lower=0> sigma_plot;             //

// ---- Hierarchical slopeS, IS THIS RIGHT?? do we just want independent slopes for DOYsqs but no partial pooling?
  real mu_DOY; //mean slope for DOY
  real mu_DOYsq; //mean slope for DOYsqs

  real<lower=0> sigma_DOY;//sd of DOY
  real<lower=0> sigma_DOYsq;

  vector[Nyr] beta_DOYs_raw; //year-specific slope values
  vector[Nyr] beta_DOYsqs_raw;
}

transformed parameters {
  vector[Nplots] u_plot=u_plot_raw * sigma_plot;      //plot-level random effect

  vector[Nyr] alpha_year = alpha+sigma_year*alpha_year_raw;

  vector[Nyr] beta_DOYs    = mu_DOY+sigma_DOY*beta_DOYs_raw;

  vector[Nyr] beta_DOYsqs_uncon = mu_DOYsq + sigma_DOYsq * beta_DOYsqs_raw;
  vector<upper=0>[Nyr] beta_DOYsqs; //constrain to bwe negative

  for (i in 1:Nyr)
    beta_DOYsqs[i] = -abs(beta_DOYsqs_uncon[i]); // or beta_DOYsqs[y] = fmin(beta_DOYsqs_uncon[y], 0);

  array [N] real eta;

  for (n in 1:N) {
       eta[n]= alpha +
               beta_DOYs[year_id[n]] * DOYs[n] +
               beta_DOYsqs[year_id[n]] * DOYsqs[n] +
               u_plot[plot_id[n]] +
               log(obs_days[n]);           // offset on log scale

}
}

model {
  alpha ~ normal(0, 5);
  alpha_year_raw ~ normal(0, 1);
  sigma_year ~ normal(0, 2);

  u_plot_raw ~ normal(0, 1);
  sigma_plot ~ normal(0, 2);

  mu_DOY   ~ normal(0, 5); //or 0,2
  mu_DOYsq ~ normal(-1, 1);      // favors concave-down shape

  sigma_DOY   ~ normal(0, 1); //or 0,2
  sigma_DOYsq ~ normal(0, 1); //or 0,2

  beta_DOYs_raw ~ normal(0, 1);
  beta_DOYsqs_raw ~ normal(0, 1);
  phi ~ exponential(1);         //apparently more stable than cauchy(0,2.5)

for (i in 1:N) {
    y[i] ~ neg_binomial_2(exp(eta[i]), phi);  //likelihood of abundance
    }
}

generated quantities {

  array[N] int y_pred;   //PPD of each obs

  for (n in 1:N) {

    y_pred[n] = neg_binomial_2_rng(exp(eta[n]), phi);
  }
}

