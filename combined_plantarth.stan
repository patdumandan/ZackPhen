// Combined Plant–Arthropod Phenology Model

data {
  int<lower=1> ND;                   // # days in prediction grid

  int<lower=1> Nyr_pl;               // # plant years
  int<lower=1> Nyr_ar;               // # arthropod years

  int<lower=1> N_shared;
  array[N_shared] int<lower=1,upper=Nyr_pl> pl_shared_id;
  array[N_shared] int<lower=1,upper=Nyr_ar> ar_shared_id;

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
  sigma_pl_plot ~ normal(0,2);
  u_pl_plot_raw ~ normal(0,1);

  alpha_ar ~ normal(0,5);
  beta_DOYs_ar ~ normal(0,2);
  beta_DOYsqs_ar ~ normal(-1,1);
  sigma_ar_plot ~ normal(0,2);
  u_ar_plot_raw ~ normal(0,1);
  phi_ar ~ cauchy(0,2.5);

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

  // Overlap only for shared years
  vector[N_shared] overlap_joint;

  // Peak DOY (raw calendar day) for each year in each taxon
  vector[Nyr_pl] peak_DOY_pl;
  vector[Nyr_ar] peak_DOY_ar;

  // Temporary curves
  vector[ND] f_pl;
  vector[ND] f_ar;

  for (s in 1:N_shared) {

    int yp = pl_shared_id[s];
    int ya = ar_shared_id[s];

    real best_pl = -1;
    real best_ar = -1;
    int best_d_pl = 1;
    int best_d_ar = 1;

    for (d in 1:ND) {

      real DOY_raw = 150 + (d - 1);

      // Plant curve
      real zpl = (DOY_raw - DOY_mean_pl) / DOY_sd_pl;
      real zpl2 = zpl * zpl;

      real lp_pl =alpha_pl +
        beta_DOYs_pl[yp]   * zpl +
        beta_DOYsqs_pl[yp] * zpl2;

      f_pl[d] = inv_logit(lp_pl);

      // Arthropod curve
      real zar = (DOY_raw - DOY_mean_ar) / DOY_sd_ar;
      real zar2 = zar * zar;

      real lp_ar = alpha_ar +
        beta_DOYs_ar[ya]   * zar +
        beta_DOYsqs_ar[ya] * zar2;

      f_ar[d] = exp(lp_ar);

      // Track peak DOYs
      if (f_pl[d] > best_pl) {
        best_pl = f_pl[d];
        best_d_pl = d;
      }

      if (f_ar[d] > best_ar) {
        best_ar = f_ar[d];
        best_d_ar = d;
      }
    }

    // Save peaks for these shared-year indices
    peak_DOY_pl[yp] = 150 + (best_d_pl - 1);
    peak_DOY_ar[ya] = 150 + (best_d_ar - 1);

    // Compute overlap
    real o = 0;
    for (d in 1:ND)
      o += fmin(f_pl[d], f_ar[d]);

    overlap_joint[s] = o;
  }

  // For plant years with no shared arthropod year → assign peak DOY via the same loop
  for (y in 1:Nyr_pl) {
    if (peak_DOY_pl[y] == 0) {   // not yet filled
      real best = -1;
      int best_d = 1;

      for (d in 1:ND) {
        real DOY_raw = 150 + (d - 1);
        real z = (DOY_raw - DOY_mean_pl) / DOY_sd_pl;
        real z2 = z * z;
        real lp =alpha_pl +
          beta_DOYs_pl[y] * z +
          beta_DOYsqs_pl[y] * z2;
        real fx = inv_logit(lp);
        if (fx > best) {
          best = fx;
          best_d = d;
        }
      }
      peak_DOY_pl[y] = 150 + (best_d - 1);
    }
  }

  // Same for arthropod years
  for (y in 1:Nyr_ar) {
    if (peak_DOY_ar[y] == 0) {
      real best = -1;
      int best_d = 1;

      for (d in 1:ND) {
        real DOY_raw = 150 + (d - 1);
        real z = (DOY_raw - DOY_mean_ar) / DOY_sd_ar;
        real z2 = z * z;
        real lp =
          alpha_ar +
          beta_DOYs_ar[y] * z +
          beta_DOYsqs_ar[y] * z2;
        real fx = exp(lp);
        if (fx > best) {
          best = fx;
          best_d = d;
        }
      }
      peak_DOY_ar[y] = 150 + (best_d - 1);
    }
  }
}
