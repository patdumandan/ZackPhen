data {
  int<lower=1> N;                      // observations
  int<lower=1> Nplots;                      // number of plots
  int<lower=1> Nyr;                      // number of years
  real DOY_sd;                           //sd of raw DOY
  real DOY_mean;                        //mean of raw DOY

  array[N] int<lower=1, upper=Nplots> plot_id;
  array[N] int<lower=1, upper=Nyr> year_id;

  array[N] int<lower=0> tot_flwr;      // successes
  array[N] int<lower=0> tot_NF;        // failures

  vector[N] DOYs;
  vector[N] DOYsqs;
}

parameters {
  real alpha;                          // global intercept

  vector[Nplots] u_plot;                    // random intercept per plot
  real<lower=0> sigma_plot;            // plot SD

  vector[Nyr] beta_DOYs;         // DOYs slope per year
  vector[Nyr] beta_DOYsqs;       // DOYsqs slope per year

}

// transformed parameters {
//   vector[J] plot_intercept = sigma_plot * z_plot;
// }

model {
    alpha ~ normal(0, 5);
    beta_DOYs ~ normal(0, 2);          // adjust SD as needed
    beta_DOYsqs ~ normal(-1, 0.3);       // negative to enforce concave-down
    u_plot ~ normal(0, sigma_plot);
    sigma_plot ~ normal(0, 2);

  // Likelihood
  for (n in 1:N) {
    int y_total = tot_flwr[n] + tot_NF[n];

    real eta = alpha +
               beta_DOYs[year_id[n]] * DOYs[n] +
               beta_DOYsqs[year_id[n]] * DOYsqs[n] +
               u_plot[plot_id[n]];

    tot_flwr[n] ~ binomial_logit(y_total, eta);
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
    y_pred[n] = binomial_rng(tot_flwr[n] + tot_NF[n], inv_logit(eta_out[n]));
  }

  for (y in 1:Nyr) {
    if (beta_DOYsqs[y] != 0) {
      DOY_peak_std[y] = -beta_DOYs[y] / (2 * beta_DOYsqs[y]);
      DOY_peak_unscaled[y] = DOY_peak_std[y] * DOY_sd + DOY_mean;
      DOY_peak_unscaled[y]=fmin(fmax(DOY_peak_unscaled[y], 150), 270);//to constrain
    //  predictions to 1 to 365
    }
    else {
      DOY_peak_std[y] = negative_infinity();
      DOY_peak_unscaled[y] = negative_infinity();
    }
  }
}
