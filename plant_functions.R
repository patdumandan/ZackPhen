rolling_mod=function(split) {

  analysis_set=analysis(split)
  Plot=analysis_set[, "Plot"]

  fit_model= lme4::glmer(analysis_set[, "DOY"]~analysis_set[, "Year"]+
                           (1|Plot),data=analysis_set)
  return(summary(fit_model))

}
rolling_plot=function(split, slope, intercept) {

  analysis_set=analysis(split)
  slope=slope
  intercept=intercept

  plot_model=ggplot2::ggplot(data=analysis_set,
                             aes(x=as.integer(analysis_set[, "Year"]), y=analysis_set[, "DOY"]))+
    geom_point()+
    geom_abline(slope=slope, intercept=intercept)+theme_classic()+
    ylab("DOY")+ xlab("Year")+
    stat_smooth(method="lm")

  print(plot_model)
}

get_year=function(split) {

  analysis_set=analysis(split)

  start_yr=as.numeric(analysis_set[,"Year"][1])

}

get_dat=function(split) {

  analysis_set=analysis(split)

  yr=as.vector(analysis_set[,"Year"])
  doys=as.vector(analysis_set[,"DOY"])

  return(cbind(yr, doys))
}

generate_fitted_curves <- function(n_years, alpha_mean, beta_doy, beta_doy_sq,
                                   doy_std, doy_sq_std, doy_seq, year_levels) {
  lapply(1:n_years, function(i) {
    eta <- alpha_mean + beta_doy[i] * doy_std + beta_doy_sq[i] * doy_sq_std
    p <- plogis(eta)

    data.frame(
      DOY = doy_seq,
      prob = p,
      year = rep(year_levels[i], length(doy_seq))
    )
  })
}
