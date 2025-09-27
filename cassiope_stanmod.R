library(cmdstanr)
library(dplyr)

# restructure data
file_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\raw"
dat_name=paste(file_path, '\\plant_datA','.csv', sep = '')

plant_datA=read.csv(dat_name, header=T, sep=',',  stringsAsFactors = F)

cas_datA <- plant_datA %>%filter(species=="Cassiope")

cass_data <- list(
  N = nrow(cas_datA),
  tot_F = cas_datA$tot_F,
  tot_NF = cas_datA$tot_NF,
  Nyr = length(unique(cas_datA$year)),
  year_id = as.integer(factor(cas_datA$year)),
  DOYs=cas_datA$DOYs,
  DOYsqs=cas_datA$DOYsqs,
  yearc=cas_datA$yearc,
  Nplots = length(unique(cas_datA$plot_id)),
  plot_id = cas_datA$plot_id,
  DOY_sd= sd(plant_datA$DOY),
  DOY_mean=mean(plant_datA$DOY)
)

#compile model
plant_mod=cmdstan_model("plant_phen.stan")

#fit model
cas_mod <- plant_mod$sample(
  data = cass_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500
)

#predictions

cas_draws <- cas_mod$draws(format="df", variables="y_pred")
cas_pred_matrix <- as.data.frame(cas_draws)
cas_pred_matrix=cas_pred_matrix[,-865:-867]

cas_pred_mean <- apply(cas_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
cas_pred_lower <- apply(cas_pred_matrix, MARGIN=2, quantile, probs = 0.025)
cas_pred_upper <- apply(cas_pred_matrix, MARGIN=2, quantile, probs = 0.975)


#plot predictions of flower totals

cas_observed <- cas_datA$tot_F
casN <- length(cas_observed)

casdf_plot <- data.frame(
  pred_mean = cas_pred_mean,
  lower = cas_pred_lower,
  upper = cas_pred_upper
)%>%cbind(cas_datA)

#to check weird years
highlight_years <- c( "1998", "2018", "1997", "2014", "2015", "2020", "2021")

ggplot(casdf_plot, aes(x = DOY, y = pred_mean, group = year, col = as.factor(year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_line(data = subset(casdf_plot, year %in% highlight_years),
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

casdraws_doy=cas_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(cas_datA$year))

require(posterior)

caspeak_df <- as_draws_df(casdraws_doy) |>
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
cassummary_peak <- caspeak_df |>
  group_by(year) |>
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%
  mutate(
    highlight_group = ifelse(year %in% highlight_years, as.character(year), "Other"))

print(cassummary_peak)

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
cas_peak=ggplot(cassummary_peak, aes(x = year, y = mean, col=as.factor(year))) +
  geom_line(data = subset(cassummary_peak, year %in% highlight_years),
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

casslope_samples <- caspeak_df |>
  group_by(.draw) |>
  summarise(
    slope = coef(lm(DOY_peak ~ year))[2]  # Extract the slope
  )

casslope_summary <- casslope_samples |>
  summarise(
    mean    = mean(slope),
    median  = median(slope),
    lower90 = quantile(slope, 0.05),
    upper90 = quantile(slope, 0.95),
    lower95 = quantile(slope, 0.025),
    upper95 = quantile(slope, 0.975)
  )

print(casslope_summary)

library(ggplot2)

cas_slope=ggplot(casslope_samples, aes(x = slope)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "slope distribution",
    x = "Change in Peak DOY per Year",
    y = "Density"
  ) +
  theme_classic()

#curvatures####
casdraws_ct=cas_mod$draws(variable="beta_DOYsqs", format="df")

casct_df <- as_draws_df(casdraws_ct) |>
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

cassummary_ct <- casct_df |>
  group_by(year) |>
  summarise(
    mean  = mean(ct),
    median= median(ct),
    lower= quantile(ct, 0.05),
    upper= quantile(ct, 0.95)
  )%>%
  mutate(
    highlight_group = ifelse(year %in% highlight_years, as.character(year), "Other"))

cas_curve=ggplot(cassummary_ct, aes(x = year, y = mean, col=as.factor(year))) +
  geom_line(data = subset(cassummary_ct, year %in% highlight_years),
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

cas_res=ggarrange(cas_peak, cas_slope, cas_curve, nrow=1, ncol=3)
annotate_figure(cas_res, top = text_grob("Cassiope", face = "bold", size = 20))
