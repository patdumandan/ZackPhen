shared_years <- intersect(unique(dry_datA$year),
                          unique(pho_datA$Year))

drypho$N_shared = length(shared_years)

drypho$pl_shared_id = match(shared_years,
                            sort(unique(dry_datA$year)))

drypho$ar_shared_id = match(shared_years,
                            sort(unique(pho_datA$Year)))
drypho=list(
DOY_mean_pl=mean(plant_datA$DOY),
DOY_sd_pl=sd(plant_datA$DOY),
Nyr_pl=length(unique(dry_datA$year)),
Nplant=nrow(dry_datA),
Nplant_plots=length(unique(dry_datA$plot_id)),
plant_plot_id=dry_datA$plot_id,
plant_year_id=as.integer(factor(dry_datA$year)),
tot_F=dry_datA$tot_F,
tot_NF=dry_datA$tot_NF,
plant_DOY=dry_datA$DOYs,
plant_DOYsqs=dry_datA$DOYsqs,
DOY_mean_ar=mean(arth_df$DOY),
DOY_sd_ar=sd(arth_df$DOY),
Nyr_ar=length(unique(pho_datA$Year)),
Narth=nrow(pho_datA),
Narth_plots=length(unique(pho_datA$plot_id)),
arth_plot_id=pho_datA$plot_id,
arth_year_id=as.integer(factor(pho_datA$Year)),
arth_y=pho_datA$TotalCatch1,
arth_effort=pho_datA$TrapDays,
arth_DOY=pho_datA$DOYs,
arth_DOYsqs=pho_datA$DOYsqs,
ND=121)


drypho_mod=cmdstan_model("combined_plantarth.stan")

#fit model
drypho_mods <- drypho_mod$sample(
  data = drypho,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 200,
  iter_warmup = 50)

library(ggplot2)
library(dplyr)
library(tidyr)
library(posterior)  # if using cmdstanr

# Example: shared_years vector
# shared_years <- c(2010, 2011, 2012, ...)  # length = N_shared

# ----- Extract posterior draws -----
overlap <- drypho_mods$draws("overlap_joint") |>
  as_draws_matrix()

# Convert to long/tidy format
overlap_df <- as.data.frame(overlap) |>
  pivot_longer(cols = everything(),
               names_to = "shared_index",
               values_to = "overlap") |>
  mutate(shared_index = as.integer(gsub("[^0-9]", "", shared_index)),
         year = shared_years[shared_index])

# Summarize for plotting (mean and 95% CI)
overlap_sum <- overlap_df |>
  group_by(year) |>
  summarise(
    mean  = mean(overlap),
    lower = quantile(overlap, 0.025),
    upper = quantile(overlap, 0.975),
    .groups = "drop"
  )

# ----- Plot -----
ggplot(overlap_sum, aes(x = year, y = mean)) +
  geom_point(size = 3) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(
    x = "Year",
    y = "Phenology Overlap",
    title = "Overlap Between Dryas & Phoridae Phenology",
    # subtitle = "Points = posterior mean, shaded = 95% credible interval"
  ) +
  theme_classic() +
  theme(text = element_text(size = 14))

