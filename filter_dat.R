library(dplyr)
library(tidyr)
library(mgcv)
library(purrr)

# Step 1: Aggregate and calculate proportions
dry_dat_sum <- dry_dat_raw %>%
  group_by(Plot, year, DOY) %>%
  summarise(
    tot_bud = sum(Buds, na.rm = TRUE),
    tot_flwr = sum(Flowers, na.rm = TRUE),
    tot_sen = sum(Senescent, na.rm = TRUE),
    tot_all = sum(tot_bud, tot_flwr, tot_sen),
    tot_NF=sum(tot_bud, tot_sen, na.rm=T),
    .groups = 'drop'
  ) %>%
  group_by(Plot, year) %>%
  mutate(
    tot_yr = sum(tot_all, na.rm = TRUE),
    tot_buds = sum(tot_bud, na.rm = TRUE),
    tot_flwrs = sum(tot_flwr, na.rm = TRUE),
    tot_sens = sum(tot_sen, na.rm = TRUE)
  ) %>%
  filter(tot_buds > 0 & tot_flwrs > 0 & tot_sens > 0) %>%
  mutate(
    prop_bud = tot_bud / tot_all,
    prop_flwr = tot_flwr / tot_all,
    prop_sen = tot_sen / tot_all
  ) %>%
  ungroup()

# Step 2: Nest by Plot-Year
nested_data <- dry_dat_sum %>%
  group_by(Plot, year) %>%
  filter(n() > 2) %>%
  nest()

# Step 3: Function to check whether GAM for a phenophase is estimable
check_single_gam <- function(df, response_var) {
  tryCatch({
    formula <- as.formula(paste0(response_var, " ~ s(DOY, k=4)"))
    model <- gam(formula, data = df, family = poisson(link="log"))
    edf_val <- summary(model)$s.table[1, "edf"]
    pred <- predict(model, newdata = df)
    pred_var <- var(pred, na.rm = TRUE)

    # Accept the model if it’s either flexible or shows any variation
    return(edf_val > 1 || pred_var > 1e-4)
  }, error = function(e) {
    return(FALSE)
  })
}

# Step 4: Apply the GAM check to all three phenophases
nested_data <- nested_data %>%
  mutate(
    bud_ok  = map_lgl(data, ~check_single_gam(.x, "tot_bud")),
    flwr_ok = map_lgl(data, ~check_single_gam(.x, "tot_flwr")),
    sen_ok  = map_lgl(data, ~check_single_gam(.x, "tot_sen")),
    all_ok  = bud_ok & flwr_ok & sen_ok
  ) %>%
  filter(flwr_ok)

# Step 5: Unnest valid data and continue with calculations
dry_dat_filtered <- nested_data %>%
  select(Plot, year, data) %>%
  unnest(cols = data) %>%
  group_by(Plot, year, DOY) %>%
  mutate(
    cumprop_10 = tot_bud / tot_yr,
    cumprop_50 = tot_flwr / tot_yr,
    cumprop_90 = tot_sen / tot_yr,
    props_10 = tot_bud / tot_buds,
    props_50 = tot_flwr / tot_flwrs,
    props_90 = tot_sen / tot_sens
  )

dry_datA=dry_dat_filtered%>%
  #select(,-c(15:20))%>%
  replace_na(list(value=0))
