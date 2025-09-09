data {
  int<lower=1> N;                      // observations
  int<lower=1> Nplots;                      // number of plots
  int<lower=1> Nyr;                      // number of years
  real DOY_sd;                           //sd of raw DOY
  real DOY_mean;                        //mean of raw DOY

  array[N] int<lower=1, upper=Nplots> plot_id;
  array[N] int<lower=1, upper=Nyr> year_id;

  array[N] int<lower=0> tot_F;      // flowers
  array[N] int<lower=0> tot_NF;        // buds+senescent

  vector[N] DOYs;                      //DOY standardized
  vector[N] DOYsqs;                    // DOY squared, standardized (quadratic term)
}

parameters {
  real alpha;                          // global intercept

  vector[Nplots] u_plot_raw;          // for non.centered paramterization of plot-level random effects
  real<lower=0> sigma_plot;

  vector[Nyr] beta_DOYs;         // DOYs slope per year
  vector <upper=0> [Nyr] beta_DOYsqs;       // DOYsqs slope per year, constrained to be +
}

transformed parameters {
  vector[Nplots] u_plot;
  u_plot = u_plot_raw * sigma_plot;
}


model {
    alpha ~ normal(0, 5);
    beta_DOYs ~ normal(0, 2);          // adjust SD as needed
    beta_DOYsqs ~ normal(-1, 1);       // negative to enforce concave-down
    sigma_plot ~ normal(0, 2);
    u_plot_raw ~ normal(0, 1);

  // Likelihood
  for (n in 1:N) {
    int y_total = tot_F[n] + tot_NF[n]; // total individuals

    real eta = alpha + // eta=probability of being flower
               beta_DOYs[year_id[n]] * DOYs[n] +
               beta_DOYsqs[year_id[n]] * DOYsqs[n] +
               u_plot[plot_id[n]];

    tot_F[n] ~ binomial_logit(y_total, eta);
  }
}

generated quantities {
  array[N] int y_pred;
  vector[N] eta_out;
  vector[Nyr] DOY_peak_std;
  vector[Nyr] DOY_peak_unscaled;

  for (n in 1:N) {
    int y = year_id[n];
    eta_out[n] = alpha +
                 beta_DOYs[y] * DOYs[n] +
                 beta_DOYsqs[y] * DOYsqs[n] +
                 u_plot[plot_id[n]];
    y_pred[n] = binomial_rng(tot_F[n] + tot_NF[n], inv_logit(eta_out[n]));
  }

  for (y in 1:Nyr) {
    if (beta_DOYsqs[y] != 0) {
      DOY_peak_std[y] = -beta_DOYs[y] / (2 * beta_DOYsqs[y]);
      DOY_peak_unscaled[y] = DOY_peak_std[y] * DOY_sd + DOY_mean;
      DOY_peak_unscaled[y]=fmin(fmax(DOY_peak_unscaled[y], 150), 270);//to constrain
    //  predictions to observation period
    }
    else {
      DOY_peak_std[y] = negative_infinity();
      DOY_peak_unscaled[y] = negative_infinity();
    }
  }
}
