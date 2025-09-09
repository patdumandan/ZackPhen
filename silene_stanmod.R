library(cmdstanr)
library(dplyr)

# restructure data
file_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen"
dat_name=paste(file_path, '\\plant_datA','.csv', sep = '')

plant_datA=read.csv(dat_name, header=T, sep=',',  stringsAsFactors = F)

sil_datA <- plant_datA %>%filter(species=="Silene")

sil_data <- list(
  N = nrow(sil_datA),
  tot_F = sil_datA$tot_F,
  tot_NF = sil_datA$tot_NF,
  Nyr = length(unique(sil_datA$year)),
  year_id = as.integer(factor(sil_datA$year)),
  DOYs=sil_datA$DOYs,
  DOYsqs=sil_datA$DOYsqs,
  yearc=sil_datA$yearc,
  Nplots = length(unique(sil_datA$plot_id)),
  plot_id = sil_datA$plot_id,
  DOY_sd= sd(plant_datA$DOY),
  DOY_mean=mean(plant_datA$DOY)
)

#compile model
plant_mod=cmdstan_model("pheno_quad.stan")

#fit model
sil_mod <- plant_mod$sample(
  data = sil_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500
)

#predictions

sil_draws <- sil_mod$draws(format="df", variables="y_pred")
sil_pred_matrix <- as.data.frame(sil_draws)
sil_pred_matrix=sil_pred_matrix[,-1107:-1109]

sil_pred_mean <- apply(sil_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
sil_pred_lower <- apply(sil_pred_matrix, MARGIN=2, quantile, probs = 0.025)
sil_pred_upper <- apply(sil_pred_matrix, MARGIN=2, quantile, probs = 0.975)


#plot predictions of flower totals

sil_observed <- sil_datA$tot_F
silN <- length(sil_observed)

sildf_plot <- data.frame(
  pred_mean = sil_pred_mean,
  lower = sil_pred_lower,
  upper = sil_pred_upper
)%>%cbind(sil_datA)

#to check weird years
highlight_years <- c( "1998", "2018", "1997", "2014", "2015", "2020", "2021")

ggplot(sildf_plot, aes(x = DOY, y = pred_mean, group = year, col = as.factor(year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_line(data = subset(sildf_plot, year %in% highlight_years),
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

sildraws_doy=sil_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(sil_datA$year))

require(posterior)

silpeak_df <- as_draws_df(sildraws_doy) |>
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
silsummary_peak <- silpeak_df |>
  group_by(year) |>
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%
  mutate(
    highlight_group = ifelse(year %in% highlight_years, as.character(year), "Other"))

print(silsummary_peak)

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
sil_peak=ggplot(silsummary_peak, aes(x = year, y = mean, col=as.factor(year))) +
  geom_line(data = subset(silsummary_peak, year %in% highlight_years),
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

silslope_samples <- silpeak_df |>
  group_by(.draw) |>
  summarise(
    slope = coef(lm(DOY_peak ~ year))[2]  # Extract the slope
  )

silslope_summary <- silslope_samples |>
  summarise(
    mean    = mean(slope),
    median  = median(slope),
    lower90 = quantile(slope, 0.05),
    upper90 = quantile(slope, 0.95),
    lower95 = quantile(slope, 0.025),
    upper95 = quantile(slope, 0.975)
  )

print(silslope_summary)

library(ggplot2)

sil_slope=ggplot(silslope_samples, aes(x = slope)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "slope distribution",
    x = "Change in Peak DOY per Year",
    y = "Density"
  ) +
  theme_classic()

#curvatures####
sildraws_ct=sil_mod$draws(variable="beta_DOYsqs", format="df")

silct_df <- as_draws_df(sildraws_ct) |>
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

silsummary_ct <- silct_df |>
  group_by(year) |>
  summarise(
    mean  = mean(ct),
    median= median(ct),
    lower= quantile(ct, 0.05),
    upper= quantile(ct, 0.95)
  )%>%
  mutate(
    highlight_group = ifelse(year %in% highlight_years, as.character(year), "Other"))

sil_curve=ggplot(silsummary_ct, aes(x = year, y = mean, col=as.factor(year))) +
  geom_line(data = subset(silsummary_ct, year %in% highlight_years),
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

sil_res=ggarrange(sil_peak, sil_slope, sil_curve, nrow=1, ncol=3)
annotate_figure(sil_res, top = text_grob("Silene", face = "bold", size = 20))
