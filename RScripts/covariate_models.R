#covar data

library(brms)
library(posterior)
library(dplyr)
library(tibble)
library(ggpubr)

#functions
fit_mod_bayes=function(dat_full, covar) {

  form=reformulate("year", response = covar)

  fit <- brm(formula = form,data = dat_full,family = gaussian(),

             prior = c(prior(normal(0, 1), class = "b")),
             chains = 4,iter = 2000,warmup = 1000,cores = 4,seed = 123,refresh = 0)

  return(fit)
}

# 2. Function to fit ONE moving window

get_slope_per_window_bayes <- function(dat,covar,template_fit,
                                       probs = c(0.025, 0.975)) {

  fit=update(template_fit,newdata = dat,recompile = FALSE,refresh = 0)

  draws <- as_draws_df(fit)

  slope_draws <- draws$b_year

  slope_mean <- mean(slope_draws)

  slope_median <- median(slope_draws)

  ci <- quantile(slope_draws,probs = probs)

  prob_positive <- mean(slope_draws > 0)

  prob_negative <- mean(slope_draws < 0)

  direction_certainty <- max(prob_positive,prob_negative)

  ci_width <- ci[2] - ci[1]

  precision <- 1 / ci_width

  tibble(covar = covar,
         start_year = min(dat$year),
         end_year = max(dat$year),
         n_years = length(unique(dat$year)),
         n_obs = nrow(dat),
         slope = slope_mean,
         slope_median = slope_median,
         lower = unname(ci[1]),
         upper = unname(ci[2]),
         prob_positive = prob_positive,
         prob_negative = prob_negative,
         direction_certainty = direction_certainty,
         ci_width = unname(ci_width),
         precision = unname(precision))
}

#spring AT####
step_size <- 1

years_all <- sort(unique(airtemp_dat_apr_jul$year))

springtemp_slide_draws <- list()

for (nyr in seq(5, 27, by = step_size)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs= years_all[start:(start + nyr - 1)]

    springtemp_slide_draws[[length(springtemp_slide_draws) + 1]] <-
      airtemp_dat_apr_jul %>%
      filter(year %in% yrs)
  }
}

#summer AT####
#summer temp(June-Aug)
summertemp_slide_draws <- list()

for (nyr in seq(5, 27, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    summertemp_slide_draws[[length(summertemp_slide_draws) + 1]] <-
      airtemp_dat_apr_jul %>% filter(year %in% yrs)
  }}

#snow cover####
step_size=1
nyr=length(unique(snow_ave$year))

years_all <- sort(unique(snow_ave$year))

snow_slide_draws <- list()

for (nyr in seq(5, 29, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    snow_slide_draws[[length(snow_slide_draws) + 1]] <-
      snow_ave %>% filter(year %in% yrs)
  }}

##spring AT####
#sanity check
fit1_springAT_full=fit_mod_bayes(dat_full = airtemp_dat_apr_jul,
                                          covar = "spring_mean_temp")

fit_springAT_win1=get_slope_per_window_bayes(dat=springtemp_slide_draws[[1]],
                                             covar="spring_mean_temp",
                                             template_fit=fit1_springAT_full)

spring_draws_df=purrr::map_dfr(
  springtemp_slide_draws,
  ~ get_slope_per_window_bayes(
    dat = .x,
    covar = "spring_mean_temp",
    template_fit = fit1_springAT_full))

spring_direction <- ggplot(
  spring_draws_df,
  aes(x = n_years, y = start_year, fill = slope)
) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "darkblue",
    mid = "grey95",
    high = "darkred",
    midpoint = 0,
    name = "Slope"
  ) +ggtitle("direction")+
  theme_classic()+
  ylab("Start year")+ xlab("Time-series length (years)")

spring_direction

spring_certainty <- ggplot(
  spring_draws_df,
  aes(x = n_years, y = start_year, fill = direction_certainty)
) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    limits = c(0.5, 1),
    name = "Direction\ncertainty"
  ) +ggtitle("certainty")+
  theme_classic()+
  ylab("Start year")+ xlab("Time-series length (years)")

spring_certainty

spring_precision <- ggplot(
  spring_draws_df,
  aes(x = n_years, y = start_year, fill = precision)
) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    name = expression(1 / "CI width")
  ) +
  theme_classic()+ggtitle("precision")+
  ylab("Start year")+ xlab(" Time-series length (years)")

spring_precision

spring_plots=ggarrange(spring_direction,
                       spring_certainty,
                       spring_precision,
                       ncol = 3,nrow = 1)

spring_plots=annotate_figure(spring_plots,
  top = text_grob("Spring air temperature",
    face = "bold",
    size = 14))

spring_plots


##summer AT####
fit1_summerAT_full=fit_mod_bayes(dat_full = airtemp_dat_apr_jul,
                                 covar = "summer_mean_temp")
#sanity check
fit_summerAT_win1=get_slope_per_window_bayes(dat=summertemp_slide_draws[[1]],
                                             covar="summer_mean_temp",
                                             template_fit=fit1_summerAT_full)
summary(fit_summerAT_win1)
summary(fit1_summerAT_full)

summer_draws_df=purrr::map_dfr(
  summertemp_slide_draws,
  ~ get_slope_per_window_bayes(
    dat = .x,
    covar = "summer_mean_temp",
    template_fit = fit1_summerAT_full))

summer_direction <- ggplot(summer_draws_df,
                           aes(x = n_years, y = start_year, fill = slope)) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "darkblue",
    mid = "grey95",
    high = "darkred",
    midpoint = 0,
    name = "Slope"
  ) +ggtitle("direction")+
  theme_classic()+
  ylab("Start year")+ xlab("Time-series length (years)")

summer_direction

summer_certainty <- ggplot(
  summer_draws_df,
  aes(x = n_years, y = start_year, fill = direction_certainty)
) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    limits = c(0.5, 1),
    name = "Direction\ncertainty"
  ) +ggtitle("certainty")+
  theme_classic()+
  ylab("Start year")+ xlab("Time-series length (years)")

summer_certainty

summer_precision <- ggplot(
  summer_draws_df,
  aes(x = n_years, y = start_year, fill = precision)
) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    name = expression(1 / "CI width")
  ) +
  theme_classic()+ggtitle("precision")+
  ylab("Start year")+ xlab(" Time-series length (years)")

summer_precision

summer_plots=ggarrange(summer_direction,
                       summer_certainty,
                       summer_precision,
                       ncol = 3,nrow = 1)

summer_plots=annotate_figure(summer_plots,
                             top = text_grob("Summer air temperature",
                                             face = "bold",
                                             size = 14))

summer_plots


##snow_cover####
fit1_snow_cover_full=fit_mod_bayes(dat_full = snow_ave,
                                   covar = "Jun10_cover")
#sanity check
fit_snow_cover_win1=get_slope_per_window_bayes(dat=snow_slide_draws[[1]],
                                               covar="snow_ave",
                                               template_fit=fit1_snow_cover_full)
summary(fit_snow_cover_win1)
summary(fit1_snow_cover_full)

snow_draws_df=purrr::map_dfr(
  snow_slide_draws,
  ~ get_slope_per_window_bayes(
    dat = .x,
    covar = "Jun10_cover",
    template_fit = fit1_snow_cover_full))

snow_direction <- ggplot(snow_draws_df,
                         aes(x = n_years, y = start_year, fill = slope)) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "darkblue",
    mid = "grey95",
    high = "darkred",
    midpoint = 0,
    name = "Slope"
  ) +ggtitle("direction")+
  theme_classic()+
  ylab("Start year")+ xlab("Time-series length (years)")

snow_direction

snow_certainty <- ggplot(
  snow_draws_df,
  aes(x = n_years, y = start_year, fill = direction_certainty)
) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    limits = c(0.5, 1),
    name = "Direction\ncertainty"
  ) +ggtitle("certainty")+
  theme_classic()+
  ylab("Start year")+ xlab("Time-series length (years)")

snow_certainty

snow_precision <- ggplot(
  snow_draws_df,
  aes(x = n_years, y = start_year, fill = precision)
) +
  geom_tile(color = "grey90", linewidth = 0.2) +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    name = expression(1 / "CI width")
  ) +
  theme_classic()+ggtitle("precision")+
  ylab("Start year")+ xlab(" Time-series length (years)")

snow_precision

snow_plots=ggarrange(snow_direction,
                     snow_certainty,
                     snow_precision,
                     ncol = 3,nrow = 1)

snow_plots=annotate_figure(snow_plots,
                           top = text_grob("Snow cover",
                                           face = "bold",
                                           size = 14))

snow_plots
