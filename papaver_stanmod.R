library(cmdstanr)
library(dplyr)

# restructure data
file_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\data"
dat_name=paste(file_path, '\\plant_datA','.csv', sep = '')

plant_datA=read.csv(dat_name, header=T, sep=',',  stringsAsFactors = F)

pap_datA <- plant_datA %>%filter(species=="Papaver")

pap_data <- list(
  N = nrow(pap_datA),
  tot_F = pap_datA$tot_F,
  tot_NF = pap_datA$tot_NF,
  Nyr = length(unique(pap_datA$year)),
  year_id = as.integer(factor(pap_datA$year)),
  DOYs=pap_datA$DOYs,
  DOYsqs=pap_datA$DOYsqs,
  yearc=pap_datA$yearc,
  Nplots = length(unique(pap_datA$plot_id)),
  plot_id = pap_datA$plot_id,
  DOY_sd= sd(plant_datA$DOY),
  DOY_mean=mean(plant_datA$DOY)
)

#compile model
plant_mod=cmdstan_model("plant_phen.stan")

#fit model
pap_mod <- plant_mod$sample(
  data = pap_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500
)

#predictions

pap_draws <- pap_mod$draws(format="df", variables="y_pred")
pap_pred_matrix <- as.data.frame(pap_draws)
pap_pred_matrix=pap_pred_matrix[,-834:-836]

pap_pred_mean <- apply(pap_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
pap_pred_lower <- apply(pap_pred_matrix, MARGIN=2, quantile, probs = 0.025)
pap_pred_upper <- apply(pap_pred_matrix, MARGIN=2, quantile, probs = 0.975)


#plot predictions of flower totals

pap_observed <- pap_datA$tot_F
papN <- length(pap_observed)

papdf_plot <- data.frame(
  pred_mean = pap_pred_mean,
  lower = pap_pred_lower,
  upper = pap_pred_upper
)%>%cbind(pap_datA)

#to check weird years
highlight_years <- c( "1998", "2018", "1997", "2014", "2015", "2020", "2021", "2016")

ggplot(papdf_plot, aes(x = DOY, y = pred_mean, group = year, col = as.factor(year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_line(data = subset(papdf_plot, year %in% highlight_years),
            aes(x = DOY, y = pred_mean, group = year, color = as.factor(year)),
            linewidth = 1.2) +  # bold lines for highlighted years
  geom_point(aes(y = tot_F), size = 1.5) +  # points for observed data
  facet_wrap(~Plot) +
  scale_color_manual(
    values = c( "1998" = "orange", "2018"="red", "1997"= "blue", "2014"="green", "2015"="pink", "2016"="black", "2021"="violet"),
    breaks = highlight_years,
    guide = guide_legend(title = "Odd Years")) +
  theme_classic() +labs(
    title = "Predicted Flowering Curve per Year",
    y = "Predicted Flower Totals",
    color = "Year")

#extract peak DOY

papdraws_doy=pap_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(pap_datA$year))

require(posterior)

pappeak_df <- as_draws_df(papdraws_doy) |>
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
papsummary_peak <- pappeak_df |>
  group_by(year) |>
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%filter(!(mean==150 | lower==150| upper==150))

print(papsummary_peak)

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
pap_peak=ggplot(papsummary_peak, aes(x = year, y = mean, col=as.factor(year))) +
  geom_line(data = subset(papsummary_peak, year %in% highlight_years),
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

papslope_samples <- pappeak_df |>
  group_by(.draw) |>
  summarise(
    slope = coef(lm(DOY_peak ~ year))[2]  # Extract the slope
  )

papslope_summary <- papslope_samples |>
  summarise(
    mean    = mean(slope),
    median  = median(slope),
    lower90 = quantile(slope, 0.05),
    upper90 = quantile(slope, 0.95),
    lower95 = quantile(slope, 0.025),
    upper95 = quantile(slope, 0.975)
  )

print(papslope_summary)

library(ggplot2)

pap_slope=ggplot(papslope_samples, aes(x = slope)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "slope distribution",
    x = "Change in Peak DOY per Year",
    y = "Density"
  ) +
  theme_classic()

#curvatures####
papdraws_ct=pap_mod$draws(variable="beta_DOYsqs", format="df")

papct_df <- as_draws_df(papdraws_ct) |>
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

papsummary_ct <- papct_df |>
  group_by(year) |>
  summarise(
    mean  = mean(ct),
    median= median(ct),
    lower= quantile(ct, 0.05),
    upper= quantile(ct, 0.95)
  )%>%
  mutate(
    highlight_group = ifelse(year %in% highlight_years, as.character(year), "Other"))

pap_curve=ggplot(papsummary_ct, aes(x = year, y = mean, col=as.factor(year))) +
  geom_line(data = subset(papsummary_ct, year %in% highlight_years),
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

pap_res=ggarrange(pap_peak, pap_slope, pap_curve, nrow=1, ncol=3)
annotate_figure(pap_res, top = text_grob("Papaver", face = "bold", size = 20))

#posterior preds of phenological curves

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

pap_res=as_draws_df(pap_mod$draws())
papyr_lvls=sort(unique(pap_datA$year))
names(papyr_lvls)=1:length(papyr_lvls)

papalpha_mean=pap_res%>%summarise(alpha=mean(alpha))%>%as.numeric(alpha)
papbeta_DOYs_mean=pap_res%>%summarise(across(starts_with("beta_DOYs["), mean))%>%unlist()
papbeta_DOYsqs_mean=pap_res%>%summarise(across(starts_with("beta_DOYsqs["), mean))%>%unlist()

DOY_mean=mean(plant_datA$DOY)
DOY_sd=sd(plant_datA$DOY)
DOY_seq=seq(150, 270, by=1)
DOY_std=(DOY_seq-DOY_mean)/DOY_sd
DOY_sq_std=DOY_std^2
n_DOY=length(DOY_std)
papnyr=length(papbeta_DOYs_mean)

papfitted_curves=generate_fitted_curves(papnyr, papalpha_mean, papbeta_DOYs_mean,
                                        papbeta_DOYsqs_mean,DOY_std, DOY_sq_std, DOY_seq, papyr_lvls)

incyears <- sort(unique(papsummary_peak$year))

papfitted_df=do.call(rbind, papfitted_curves)%>%filter(year %in% incyears)


papp=ggplot(papfitted_df, aes(x=DOY, y=prob, col=as.factor(year)))+
  geom_line(linewidth=0.6, alpha=2)+theme_classic()+
  labs(x="DOY", y="P(flower)", title="Papaver")+
  scale_color_viridis_d()+
  xlim(150,270)

