// Combined Plant–Arthropod Phenology Model

data {
  int<lower=1> ND;                   //  days in prediction grid=121

  int<lower=1> Nyr_pl;               //  plant years
  int<lower=1> Nyr_ar;               //  arthropod years

  // int<lower=1> N_shared;
  // array[N_shared] int<lower=1,upper=Nyr_pl> pl_shared_id;
  // array[N_shared] int<lower=1,upper=Nyr_ar> ar_shared_id;

  real DOY_mean_pl;
  real DOY_sd_pl;

  real DOY_mean_ar;
  real DOY_sd_ar;

  int<lower=1> Nplant;
  int<lower=1> Nplant_plots;
  array[Nplant] int<lower=1, upper=Nplant_plots> plant_plot_id;
  array[Nplant] int<lower=1, upper=Nyr_pl> plant_year_id;
  array[Nplant] int<lower=0> tot_F;
  array[Nplant] int<lower=0> tot_NF;
  vector[Nplant] plant_DOY;
  vector[Nplant] plant_DOYsqs;

  int<lower=1> Narth;
  int<lower=1> Narth_plots;
  array[Narth] int<lower=1, upper=Narth_plots> arth_plot_id;
  array[Narth] int<lower=1, upper=Nyr_ar> arth_year_id;
  array[Narth] int<lower=0> arth_y;
  array[Narth] real<lower=0> arth_effort;
  vector[Narth] arth_DOY;
  vector[Narth] arth_DOYsqs;
}

parameters {

  real alpha_pl;
  vector[Nyr_pl] beta_DOYs_pl;
  vector<upper=0>[Nyr_pl] beta_DOYsqs_pl;

  real<lower=0> sigma_pl_plot;
  vector[Nplant_plots] u_pl_plot_raw;

  real alpha_ar;
  vector[Nyr_ar] beta_DOYs_ar;
  vector<upper=0>[Nyr_ar] beta_DOYsqs_ar;

  real<lower=0> sigma_ar_plot;
  vector[Narth_plots] u_ar_plot_raw;

  real<lower=0> phi_ar;
}

transformed parameters {

  vector[Nplant_plots] u_pl_plot = sigma_pl_plot * u_pl_plot_raw;
  vector[Narth_plots] u_ar_plot = sigma_ar_plot * u_ar_plot_raw;

  array[Nplant] real eta_pl;
  array[Narth] real eta_ar;

  for (n in 1:Nplant) {
    eta_pl[n] =alpha_pl +
              beta_DOYs_pl[ plant_year_id[n] ]   * plant_DOY[n] +
              beta_DOYsqs_pl[ plant_year_id[n] ] * plant_DOYsqs[n] +
              u_pl_plot[ plant_plot_id[n] ];
  }

  for (n in 1:Narth) {
    eta_ar[n] = alpha_ar +
                beta_DOYs_ar[ arth_year_id[n] ]   * arth_DOY[n] +
                beta_DOYsqs_ar[ arth_year_id[n] ] * arth_DOYsqs[n] +
                u_ar_plot[ arth_plot_id[n] ] + log(arth_effort[n]);
  }
}

model {

  alpha_pl ~ normal(0,5);
  beta_DOYs_pl ~ normal(0,2);
  beta_DOYsqs_pl ~ normal(-1,1);
  sigma_pl_plot ~ normal(0,1);
  u_pl_plot_raw ~ normal(0,1);

  alpha_ar ~ normal(0,5);
  beta_DOYs_ar ~ normal(0,2);
  beta_DOYsqs_ar ~ normal(-1,1);
  sigma_ar_plot ~ normal(0,1);
  u_ar_plot_raw ~ normal(0,1);
  phi_ar ~ exponential(1);

  // Likelihood: Plants
  for (i in 1:Nplant) {
    int Ntot = tot_F[i] + tot_NF[i];
    tot_F[i] ~ binomial_logit(Ntot, eta_pl[i]);
  }

  // Likelihood: Arthropods
  for (i in 1:Narth) {
    arth_y[i] ~ neg_binomial_2(exp(eta_ar[i]), phi_ar);
  }
}

generated quantities {

  array[Nplant] int y_plant;
  vector[Nyr_pl] plant_DOY_peak_std;
  vector[Nyr_pl] plant_DOY_peak_unscaled;

  array[Narth] int y_arth;
  vector[Nyr_ar] arth_DOY_peak_std;
  vector[Nyr_ar] arth_DOY_peak_unscaled;

//plant peak DOY
  for (j in 1:Nplant) {
    y_plant[j] = binomial_rng(tot_F[j] + tot_NF[j], inv_logit(eta_pl[j]));
  }

  for (y in 1:Nyr_pl) {
    if (beta_DOYsqs_pl[y] != 0) {
      plant_DOY_peak_std[y] = -beta_DOYs_pl[y] / (2 * beta_DOYsqs_pl[y]);
      plant_DOY_peak_unscaled[y] = plant_DOY_peak_std[y] * DOY_sd_pl + DOY_mean_pl;

      //constrain to observation period
      plant_DOY_peak_unscaled[y]=fmin(fmax(plant_DOY_peak_unscaled[y], 150), 270);
    }
    else {
      plant_DOY_peak_std[y] = negative_infinity();
      plant_DOY_peak_unscaled[y] = negative_infinity();
    }
  }

  //arth peak DOY
   for (n in 1:Narth) {

    y_arth[n] = neg_binomial_2_rng(exp(eta_ar[n]), phi_ar);
  }

  for (w in 1:Nyr_ar) {
    if (beta_DOYsqs_ar[w] != 0) {
      arth_DOY_peak_std[w] = -beta_DOYs_ar[w] / (2 * beta_DOYsqs_ar[w]);
      arth_DOY_peak_unscaled[w] = arth_DOY_peak_std[w] * DOY_sd_ar + DOY_mean_ar;

      // Constrain to observed period
      arth_DOY_peak_unscaled[w] = fmin(fmax(arth_DOY_peak_unscaled[w], 150), 270);
    } else {
      arth_DOY_peak_std[w] = negative_infinity(); //-Inf;forbid estimation of improbable values
      arth_DOY_peak_unscaled[w] = negative_infinity();
    }
  }
}
