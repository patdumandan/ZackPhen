library(cmdstanr)
library(dplyr)

# restructure data
file_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen"
dat_name=paste(file_path, '\\plant_datA','.csv', sep = '')

plant_datA=read.csv(dat_name, header=T, sep=',',  stringsAsFactors = F)

sal_datA <- plant_datA %>%filter(species=="Salix")

sal_data <- list(
  N = nrow(sal_datA),
  tot_F = sal_datA$tot_F,
  tot_NF = sal_datA$tot_NF,
  Nyr = length(unique(sal_datA$year)),
  year_id = as.integer(factor(sal_datA$year)),
  DOYs=sal_datA$DOYs,
  DOYsqs=sal_datA$DOYsqs,
  yearc=sal_datA$yearc,
  Nplots = length(unique(sal_datA$plot_id)),
  plot_id = sal_datA$plot_id,
  DOY_sd= sd(plant_datA$DOY),
  DOY_mean=mean(plant_datA$DOY)
)

#compile model
plant_mod=cmdstan_model("plant_phen.stan")

#fit model
sal_mod <- plant_mod$sample(
  data = sal_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500)

#predictions

sal_draws <- sal_mod$draws(format="df", variables="y_pred")
sal_pred_matrix <- as.data.frame(sal_draws)
sal_pred_matrix=sal_pred_matrix[,-939:-941]

sal_pred_mean <- apply(sal_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
sal_pred_lower <- apply(sal_pred_matrix, MARGIN=2, quantile, probs = 0.025)
sal_pred_upper <- apply(sal_pred_matrix, MARGIN=2, quantile, probs = 0.975)


#plot predictions of flower totals

sal_observed <- sal_datA$tot_F
salN <- length(sal_observed)

saldf_plot <- data.frame(
  pred_mean = sal_pred_mean,
  lower = sal_pred_lower,
  upper = sal_pred_upper
)%>%cbind(sal_datA)

#to check weird years
highlight_years <- c( "1998", "2018", "1997", "2014", "2015", "2020", "2021")

ggplot(saldf_plot, aes(x = DOY, y = pred_mean, group = year, col = as.factor(year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_line(data = subset(saldf_plot, year %in% highlight_years),
            aes(x = DOY, y = pred_mean, group = year, color = as.factor(year)),
            linewidth = 1.2) +  # bold lines for highlighted years
  geom_point(aes(y = tot_F), size = 1.5) +  # points for observed data
  facet_wrap(~Plot) +
  scale_color_manual(
    values = c( "1998" = "orange", "2018"="red", "1997"= "blue", "2014"="green", "2015"="pink", "2020"="black", "2021"="violet"),
    breaks = highlight_years,
    guide = guide_legend(title = "Odd Years")) +
  theme_classic() +labs(
    title = "Predicted Flowering Curve per Year",
    y = "Predicted Flower Totals",
    color = "Year")

#extract peak DOY

saldraws_doy=sal_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(sal_datA$year))

require(posterior)

salpeak_df <- as_draws_df(saldraws_doy) |>
  mutate(.draw = row_number()) |>
  tidyr::pivot_longer(
    cols = starts_with("DOY_peak_unscaled"),
    names_pattern = "DOY_peak_unscaled\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak"
  ) |>
  mutate(
    year_index = as.integer(year_index),
    year = years[year_index]
  )

# peak timing summary table
salsummary_peak <- salpeak_df |>
  group_by(year) |>
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%
  mutate(
    highlight_group = ifelse(year %in% highlight_years, as.character(year), "Other"))

print(salsummary_peak)

# #HPDI
# require(bayestestR)
#
# summary_peak_hpdi <- peak_df |>
#   group_by(year) |>
#   summarise(
#     median = median(DOY_peak),
#     hpdi = list(hdi(DOY_peak, ci = 0.90)))%>% # hdi() returns a data frame
#    mutate(
#     highlight_group = ifelse(year %in% highlight_years, as.character(year), "Other"))%>%
#   unnest_wider(hpdi)  # expands the "hdi" list column into lower/upper

# visualise
sal_peak=ggplot(salsummary_peak, aes(x = year, y = mean, col=as.factor(year))) +
  geom_line(data = subset(salsummary_peak, year %in% highlight_years),
            aes(x = year, y = median, group = year, color = as.factor(year)),
            linewidth = 1.2) +  # bold lines for highlighted years
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
  scale_color_manual(
    values = c( "1998" = "orange", "2018"="red", "1997"= "blue", "2014"="green", "2015"="pink", "2020"="black", "2021"="violet"),
    breaks = highlight_years,
    guide = guide_legend(title = "Odd Years")) +labs(
      title = "peak timing",
      x = "Year",
      y = "Peak Day of Year (DOY)",
      caption = "Mean and 90% CI")+
  theme_classic() #1998 uncertainty (coldest yearand longest flowering duration): https://www.nature.com/articles/nclimate1909

salslope_samples <- salpeak_df |>
  group_by(.draw) |>
  summarise(
    slope = coef(lm(DOY_peak ~ year))[2]  # Extract the slope
  )

salslope_summary <- salslope_samples |>
  summarise(
    mean    = mean(slope),
    median  = median(slope),
    lower90 = quantile(slope, 0.05),
    upper90 = quantile(slope, 0.95),
    lower95 = quantile(slope, 0.025),
    upper95 = quantile(slope, 0.975)
  )

print(salslope_summary)

library(ggplot2)

sal_slope=ggplot(salslope_samples, aes(x = slope)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "slope distribution",
    x = "Change in Peak DOY per Year",
    y = "Density"
  ) +
  theme_classic()

#curvatures####
saldraws_ct=sal_mod$draws(variable="beta_DOYsqs", format="df")

salct_df <- as_draws_df(saldraws_ct) |>
  mutate(.draw = row_number()) |>
  tidyr::pivot_longer(
    cols = starts_with("beta_DOYsqs"),
    names_pattern = "beta_DOYsqs\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "ct"
  ) |>
  mutate(
    year_index = as.integer(year_index),
    year = years[year_index]
  )

salsummary_ct <- salct_df |>
  group_by(year) |>
  summarise(
    mean  = mean(ct),
    median= median(ct),
    lower= quantile(ct, 0.05),
    upper= quantile(ct, 0.95)
  )%>%
  mutate(
    highlight_group = ifelse(year %in% highlight_years, as.character(year), "Other"))

sal_curve=ggplot(salsummary_ct, aes(x = year, y = mean, col=as.factor(year))) +
  geom_line(data = subset(salsummary_ct, year %in% highlight_years),
            aes(x = year, y = mean, group = year, color = as.factor(year)),
            linewidth = 1.2) +  # bold lines for highlighted years
  geom_point(size = 2) +
  # stat_smooth(aes(x = year, y = mean, group = 1), method = "lm",
  #             color = "black", se = TRUE, linewidth = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
  scale_color_manual(
    values = c( "1998" = "orange", "2018"="red", "1997"= "blue", "2014"="green", "2015"="pink", "2020"="black", "2021"="violet"),
    breaks = highlight_years,
    guide = guide_legend(title = "Odd Years")) +labs(
      title = "curvature",
      x = "Year",
      y = "curvature (beta_DOYsqs)",
      caption = "Mean and 90% CI")+
  theme_classic()

sal_res=ggarrange(sal_peak, sal_slope, sal_curve, nrow=1, ncol=3)
annotate_figure(sal_res, top = text_grob("Salix", face = "bold", size = 20))
