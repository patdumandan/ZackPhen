#arthropod analysis

#load packages
require(dplyr)
require(tidyr)
require(lme4)
require(ggplot2)
require(ggpubr)

#load cleana arthropod data
file_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\data"
art_name=paste(file_path, '\\arth_dat_clean','.csv', sep = '')

arth_dat=read.csv(art_name, header=T, sep=',',  stringsAsFactors = F)
arth_dat=arth_dat%>%select(-c(6:23))

#inspect data crudely
#plots per taxon

arth_plots=arth_dat%>%group_by(HoyeTaxon, Year)%>%
  summarise(plot_ids=list(unique(Plot.ID)), n_plots=length(unique(Plot.ID)))

#remove Art1, Art4 and Art6

#Art1: window trap
#Art4 and 6:discontinuous opening
#Arth 7 missing from 1996-1998?

arth_data=arth_dat%>%filter(!Plot.ID%in%c("Art1","Art4", "Art6"), !Year<1996)%>%
  group_by(HoyeTaxon) %>%select(-1)%>%
  mutate(plot_id = dense_rank(Plot.ID)) %>%
  ungroup()

#check plot and species totals, remove those with <50, and those with <5 years of data

arth_tot=arth_data%>%
  group_by(Plot.ID, Year, HoyeTaxon)%>%
  summarise(tot=sum(TotalCatch1), .groups = "drop")%>%
  filter(tot>=50)%>%
  group_by(Plot.ID, HoyeTaxon)%>%
  mutate(nyr = n_distinct(Year), .groups = "drop")%>%
  filter(nyr>=5)

arth_data50 <- arth_data %>%
  semi_join(arth_tot, by = c("Plot.ID", "Year", "HoyeTaxon"))

#check for unimodality
#include only plots and years that pass Hartigan's dip test for unimodality
library(diptest)

unimod_results <- arth_data50 %>%
  group_by(HoyeTaxon, Plot.ID, Year) %>%
  summarise(
    D_statistic = tryCatch({
      if (length(TotalCatch1) > 5) dip.test(TotalCatch1)$statistic else NA
    }, error = function(e) NA),
    p_value = tryCatch({
      if (length(TotalCatch1) > 5) dip.test(TotalCatch1)$p.value else NA
    }, error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(
    unimodal = ifelse(!is.na(p_value), p_value > 0.05, NA),
    Include = ifelse(unimodal == TRUE, 1, 0))

arth_df=arth_data50%>%
  semi_join(unimod_results, by=c("HoyeTaxon", "Year", "Plot.ID"))

#rescale data
arth_df$yearc=as.integer(factor(arth_df$Year))
#arth_data$years = (arth_data$Year - mean(arth_data$Year))/sd(arth_data$Year)
arth_df$DOYs = scale(arth_df$DOY, center = TRUE, scale = TRUE)[,1]
arth_df$DOYsqs = arth_df$DOYs^2
arth_df$TrapDays[arth_df$TrapDays <= 0] <- 0.001

#unique(arth_df$Plot.ID)
#unique(arth_df$HoyeTaxon)

#models####

#muscid####

mus_datA <- arth_df %>%filter(HoyeTaxon=="Muscidae")

mus_data <- list(
  N = nrow(mus_datA),
  y = mus_datA$TotalCatch1,
  obs_days=mus_datA$TrapDays,
  Nyr = length(unique(mus_datA$Year)),
  year_id = as.integer(factor(mus_datA$Year)),
  DOYs=mus_datA$DOYs,
  DOYsqs=mus_datA$DOYsqs,
  yearc=mus_datA$yearc,
  Nplots = length(unique(mus_datA$plot_id)),
  plot_id = mus_datA$plot_id,
  DOY_sd= sd(arth_df$DOY),
  DOY_mean=mean(arth_df$DOY))

#compile model
arth_mod=cmdstan_model("arthropod_phen.stan")

#fit model
mus_mod <- arth_mod$sample(
  data = mus_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500)

mus_draws <- mus_mod$draws(format="df", variables="y_pred")
mus_pred_matrix <- as.data.frame(mus_draws)
mus_pred_matrix=mus_pred_matrix[,-1436:-1438]

mus_pred_mean <- apply(mus_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
mus_pred_lower <- apply(mus_pred_matrix, MARGIN=2, quantile, probs = 0.025)
mus_pred_upper <- apply(mus_pred_matrix, MARGIN=2, quantile, probs = 0.975)


mus_observed <- mus_datA$TotalCatch1
musN <- length(mus_observed)

musdf_plot <- data.frame(
  pred_mean = mus_pred_mean,
  lower = mus_pred_lower,
  upper = mus_pred_upper
)%>%cbind(mus_datA)

ggplot(musdf_plot, aes(x = DOY, y = pred_mean, group = Year, col = as.factor(Year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_point(aes(y = TotalCatch1), size = 1.5) +  # points for observed data
  facet_wrap(~Plot.ID) +
  theme_classic() +labs(
    title = "Muscidae",
    y = "Predicted Abudance",
    color = "Year")

#extract peak DOY

musdraws_doy=mus_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(mus_datA$Year))

require(posterior)

muspeak_df <- as_draws_df(musdraws_doy)  %>%
  mutate(.draw = row_number())  %>%
  tidyr::pivot_longer(
    cols = starts_with("DOY_peak_unscaled"),
    names_pattern = "DOY_peak_unscaled\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak") %>%
  mutate(
    year_index = as.integer(year_index),
    year = years[year_index])

# peak timing summary table
mussummary_peak <- muspeak_df  %>%
  group_by(year)  %>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%filter(!(mean==150 | lower==150| upper==150))

print(mussummary_peak)

#posterior preds of phenological curves

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

mus_res=as_draws_df(mus_mod$draws())
musyr_lvls=sort(unique(mus_datA$Year))
names(musyr_lvls)=1:length(musyr_lvls)

musalpha_mean=mus_res%>%summarise(alpha=mean(alpha))%>%as.numeric(alpha)
musbeta_DOYs_mean=mus_res%>%summarise(across(starts_with("beta_DOYs["), mean))%>%unlist()
musbeta_DOYsqs_mean=mus_res%>%summarise(across(starts_with("beta_DOYsqs["), mean))%>%unlist()

DOY_mean=mean(arth_df$DOY)
DOY_sd=sd(arth_df$DOY)
DOY_seq=seq(150, 270, by=1)
DOY_std=(DOY_seq-DOY_mean)/DOY_sd
DOY_sq_std=DOY_std^2
n_DOY=length(DOY_std)
musnyr=length(musbeta_DOYs_mean)

musfitted_curves=generate_fitted_curves(musnyr, musalpha_mean, musbeta_DOYs_mean,
                                        musbeta_DOYsqs_mean,DOY_std, DOY_sq_std, DOY_seq, musyr_lvls)

incyears <- sort(unique(mussummary_peak$year))

musfitted_df=do.call(rbind, musfitted_curves)%>%filter(year %in% incyears)%>%
  mutate(year=as.factor(year))

musp=ggplot(musfitted_df, aes(x=DOY, y=prob, col=year))+
  geom_line(linewidth=0.6, alpha=2)+theme_classic()+
  labs(x="DOY", y="P(abundance)", title="Muscidae")+
  scale_color_viridis_d()+
  xlim(150,270)

#collembola####


col_datA <- arth_df %>%filter(HoyeTaxon=="Collembola")

col_data <- list(
  N = nrow(col_datA),
  y = col_datA$TotalCatch1,
  obs_days=col_datA$TrapDays,
  Nyr = length(unique(col_datA$Year)),
  year_id = as.integer(factor(col_datA$Year)),
  DOYs=col_datA$DOYs,
  DOYsqs=col_datA$DOYsqs,
  yearc=col_datA$yearc,
  Nplots = length(unique(col_datA$plot_id)),
  plot_id = col_datA$plot_id,
  DOY_sd= sd(arth_df$DOY),
  DOY_mean=mean(arth_df$DOY))

#compile model
arth_mod=cmdstan_model("arthropod_phen.stan")

#fit model
col_mod <- arth_mod$sample(
  data = col_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500)

col_draws <- col_mod$draws(format="df", variables="y_pred")
col_pred_matrix <- as.data.frame(col_draws)
col_pred_matrix=col_pred_matrix[,-1405:-1407]

col_pred_mean <- apply(col_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
col_pred_lower <- apply(col_pred_matrix, MARGIN=2, quantile, probs = 0.025)
col_pred_upper <- apply(col_pred_matrix, MARGIN=2, quantile, probs = 0.975)


col_observed <- col_datA$TotalCatch1
colN <- length(col_observed)

coldf_plot <- data.frame(
  pred_mean = col_pred_mean,
  lower = col_pred_lower,
  upper = col_pred_upper
)%>%cbind(col_datA)

ggplot(coldf_plot, aes(x = DOY, y = pred_mean, group = Year, col = as.factor(Year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_point(aes(y = TotalCatch1), size = 1.5) +  # points for observed data
  facet_wrap(~Plot.ID) +
  theme_classic() +labs(
    title = "Collembola",
    y = "Predicted Abudance",
    color = "Year")

#extract peak DOY

coldraws_doy=col_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(col_datA$Year))

require(posterior)

colpeak_df <- as_draws_df(coldraws_doy)  %>%
  mutate(.draw = row_number())  %>%
  tidyr::pivot_longer(
    cols = starts_with("DOY_peak_unscaled"),
    names_pattern = "DOY_peak_unscaled\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak") %>%
  mutate(
    year_index = as.integer(year_index),
    year = years[year_index])

# peak timing summary table
colsummary_peak <- colpeak_df  %>%
  group_by(year)  %>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%filter(!(mean==150 | lower==150| upper==150))

print(colsummary_peak)

#posterior preds of phenological curves

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

col_res=as_draws_df(col_mod$draws())
colyr_lvls=sort(unique(col_datA$Year))
names(colyr_lvls)=1:length(colyr_lvls)

colalpha_mean=col_res%>%summarise(alpha=mean(alpha))%>%as.numeric(alpha)
colbeta_DOYs_mean=col_res%>%summarise(across(starts_with("beta_DOYs["), mean))%>%unlist()
colbeta_DOYsqs_mean=col_res%>%summarise(across(starts_with("beta_DOYsqs["), mean))%>%unlist()

DOY_mean=mean(arth_df$DOY)
DOY_sd=sd(arth_df$DOY)
DOY_seq=seq(150, 270, by=1)
DOY_std=(DOY_seq-DOY_mean)/DOY_sd
DOY_sq_std=DOY_std^2
n_DOY=length(DOY_std)
colnyr=length(colbeta_DOYs_mean)

colfitted_curves=generate_fitted_curves(colnyr, colalpha_mean, colbeta_DOYs_mean,
                                        colbeta_DOYsqs_mean,DOY_std, DOY_sq_std, DOY_seq, colyr_lvls)

incyears <- sort(unique(colsummary_peak$year))

colfitted_df=do.call(rbind, colfitted_curves)%>%filter(year %in% incyears)%>%
  mutate(year=as.factor(year))

colp=ggplot(colfitted_df, aes(x=DOY, y=prob, col=year))+
  geom_line(linewidth=0.6, alpha=2)+theme_classic()+
  labs(x="DOY", y="P(abundance)", title="Collembola")+
  scale_color_viridis_d()+
  xlim(150,270)

#chironomid####

chi_datA <- arth_df %>%filter(HoyeTaxon=="Chironomidae")

chi_data <- list(
  N = nrow(chi_datA),
  y = chi_datA$TotalCatch1,
  obs_days=chi_datA$TrapDays,
  Nyr = length(unique(chi_datA$Year)),
  year_id = as.integer(factor(chi_datA$Year)),
  DOYs=chi_datA$DOYs,
  DOYsqs=chi_datA$DOYsqs,
  yearc=chi_datA$yearc,
  Nplots = length(unique(chi_datA$plot_id)),
  plot_id = chi_datA$plot_id,
  DOY_sd= sd(arth_df$DOY),
  DOY_mean=mean(arth_df$DOY))

#compile model
arth_mod=cmdstan_model("arthropod_phen.stan")

#fit model
chi_mod <- arth_mod$sample(
  data = chi_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500
)

chi_draws <- chi_mod$draws(format="df", variables="y_pred")
chi_pred_matrix <- as.data.frame(chi_draws)
chi_pred_matrix=chi_pred_matrix[,-1228:-1230]

chi_pred_mean <- apply(chi_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
chi_pred_lower <- apply(chi_pred_matrix, MARGIN=2, quantile, probs = 0.025)
chi_pred_upper <- apply(chi_pred_matrix, MARGIN=2, quantile, probs = 0.975)

chi_observed <- chi_datA$TotalCatch1
chiN <- length(chi_observed)

chidf_plot <- data.frame(
  pred_mean = chi_pred_mean,
  lower = chi_pred_lower,
  upper = chi_pred_upper
)%>%cbind(chi_datA)

ggplot(chidf_plot, aes(x = DOY, y = pred_mean, group = Year, col = as.factor(Year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_point(aes(y = TotalCatch1), size = 1.5) +  # points for observed data
  facet_wrap(~Plot.ID) +
  theme_classic() +labs(
    title = "Chironomidae",
    y = "Predicted Abudance",
    color = "Year")

#extract peak DOY

chidraws_doy=chi_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(chi_datA$Year))

require(posterior)

chipeak_df <- as_draws_df(chidraws_doy)  %>%
  mutate(.draw = row_number())  %>%
  tidyr::pivot_longer(
    cols = starts_with("DOY_peak_unscaled"),
    names_pattern = "DOY_peak_unscaled\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak") %>%
  mutate(
    year_index = as.integer(year_index),
    year = years[year_index])

# peak timing summary table
chisummary_peak <- chipeak_df  %>%
  group_by(year)  %>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%filter(!(mean==150 | lower==150| upper==150))

print(chisummary_peak)

#posterior preds of phenological curves

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

chi_res=as_draws_df(chi_mod$draws())
chiyr_lvls=sort(unique(chi_datA$Year))
names(chiyr_lvls)=1:length(chiyr_lvls)

chialpha_mean=chi_res%>%summarise(alpha=mean(alpha))%>%as.numeric(alpha)
chibeta_DOYs_mean=chi_res%>%summarise(across(starts_with("beta_DOYs["), mean))%>%unlist()
chibeta_DOYsqs_mean=chi_res%>%summarise(across(starts_with("beta_DOYsqs["), mean))%>%unlist()

DOY_mean=mean(arth_df$DOY)
DOY_sd=sd(arth_df$DOY)
DOY_seq=seq(150, 270, by=1)
DOY_std=(DOY_seq-DOY_mean)/DOY_sd
DOY_sq_std=DOY_std^2
n_DOY=length(DOY_std)
chinyr=length(chibeta_DOYs_mean)

chifitted_curves=generate_fitted_curves(chinyr, chialpha_mean, chibeta_DOYs_mean,
                                        chibeta_DOYsqs_mean,DOY_std, DOY_sq_std, DOY_seq, chiyr_lvls)

incyears <- sort(unique(chisummary_peak$year))

chifitted_df=do.call(rbind, chifitted_curves)%>%filter(year %in% incyears)%>%
  mutate(year=as.factor(year))

chip=ggplot(chifitted_df, aes(x=DOY, y=prob, col=year))+
  geom_line(linewidth=0.6, alpha=2)+theme_classic()+
  labs(x="DOY", y="P(abundance)", title="Chironomidae")+
  scale_color_viridis_d()+
  xlim(150,270)

#ichnumonidae####
unique(arth_df$HoyeTaxon)

ich_datA <- arth_df %>%filter(HoyeTaxon=="Ichneumonidae")

ich_data <- list(
  N = nrow(ich_datA),
  y = ich_datA$TotalCatch1,
  obs_days=ich_datA$TrapDays,
  Nyr = length(unique(ich_datA$Year)),
  year_id = as.integer(factor(ich_datA$Year)),
  DOYs=ich_datA$DOYs,
  DOYsqs=ich_datA$DOYsqs,
  yearc=ich_datA$yearc,
  Nplots = length(unique(ich_datA$plot_id)),
  plot_id = ich_datA$plot_id,
  DOY_sd= sd(arth_df$DOY),
  DOY_mean=mean(arth_df$DOY))

#compile model
arth_mod=cmdstan_model("arthropod_phen.stan")

#fit model
ich_mod <- arth_mod$sample(
  data = ich_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500)

ich_draws <- ich_mod$draws(format="df", variables="y_pred")
ich_pred_matrix <- as.data.frame(ich_draws)
ich_pred_matrix=ich_pred_matrix[,-797:-799]

ich_pred_mean <- apply(ich_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
ich_pred_lower <- apply(ich_pred_matrix, MARGIN=2, quantile, probs = 0.025)
ich_pred_upper <- apply(ich_pred_matrix, MARGIN=2, quantile, probs = 0.975)

ich_observed <- ich_datA$TotalCatch1
ichN <- length(ich_observed)

ichdf_plot <- data.frame(
  pred_mean = ich_pred_mean,
  lower = ich_pred_lower,
  upper = ich_pred_upper
)%>%cbind(ich_datA)

ggplot(ichdf_plot, aes(x = DOY, y = pred_mean, group = Year, col = as.factor(Year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_point(aes(y = TotalCatch1), size = 1.5) +  # points for observed data
  facet_wrap(~Plot.ID) +
  theme_classic() +labs(
    title = "ichronomidae",
    y = "Predicted Abudance",
    color = "Year")

#extract peak DOY

ichdraws_doy=ich_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(ich_datA$Year))

require(posterior)

ichpeak_df <- as_draws_df(ichdraws_doy)  %>%
  mutate(.draw = row_number())  %>%
  tidyr::pivot_longer(
    cols = starts_with("DOY_peak_unscaled"),
    names_pattern = "DOY_peak_unscaled\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak") %>%
  mutate(
    year_index = as.integer(year_index),
    year = years[year_index])

# peak timing summary table
ichsummary_peak <- ichpeak_df  %>%
  group_by(year)  %>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%filter(!(mean==150 | lower==150| upper==150))

print(ichsummary_peak)

#posterior preds of phenological curves

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

ich_res=as_draws_df(ich_mod$draws())
ichyr_lvls=sort(unique(ich_datA$Year))
names(ichyr_lvls)=1:length(ichyr_lvls)

ichalpha_mean=ich_res%>%summarise(alpha=mean(alpha))%>%as.numeric(alpha)
ichbeta_DOYs_mean=ich_res%>%summarise(across(starts_with("beta_DOYs["), mean))%>%unlist()
ichbeta_DOYsqs_mean=ich_res%>%summarise(across(starts_with("beta_DOYsqs["), mean))%>%unlist()

DOY_mean=mean(arth_df$DOY)
DOY_sd=sd(arth_df$DOY)
DOY_seq=seq(150, 270, by=1)
DOY_std=(DOY_seq-DOY_mean)/DOY_sd
DOY_sq_std=DOY_std^2
n_DOY=length(DOY_std)
ichnyr=length(ichbeta_DOYs_mean)

ichfitted_curves=generate_fitted_curves(ichnyr, ichalpha_mean, ichbeta_DOYs_mean,
                                        ichbeta_DOYsqs_mean,DOY_std, DOY_sq_std, DOY_seq, ichyr_lvls)

incyears <- sort(unique(ichsummary_peak$year))

ichfitted_df=do.call(rbind, ichfitted_curves)%>%filter(year %in% incyears)%>%
  mutate(year=as.factor(year))

ichp=ggplot(ichfitted_df, aes(x=DOY, y=prob, col=year))+
  geom_line(linewidth=0.6, alpha=2)+theme_classic()+
  labs(x="DOY", y="P(abundance)", title="Ichneumonidae")+
  scale_color_viridis_d()+
  xlim(150,270)

#linyphiidae

unique(arth_df$HoyeTaxon)

lin_datA <- arth_df %>%filter(HoyeTaxon=="Linyphiidae")

lin_data <- list(
  N = nrow(lin_datA),
  y = lin_datA$TotalCatch1,
  obs_days=lin_datA$TrapDays,
  Nyr = length(unique(lin_datA$Year)),
  year_id = as.integer(factor(lin_datA$Year)),
  DOYs=lin_datA$DOYs,
  DOYsqs=lin_datA$DOYsqs,
  yearc=lin_datA$yearc,
  Nplots = length(unique(lin_datA$plot_id)),
  plot_id = lin_datA$plot_id,
  DOY_sd= sd(arth_df$DOY),
  DOY_mean=mean(arth_df$DOY))

#compile model
arth_mod=cmdstan_model("arthropod_phen.stan")

#fit model
lin_mod <- arth_mod$sample(
  data = lin_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500)

lin_draws <- lin_mod$draws(format="df", variables="y_pred")
lin_pred_matrix <- as.data.frame(lin_draws)
lin_pred_matrix=lin_pred_matrix[,-943:-945]

lin_pred_mean <- apply(lin_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
lin_pred_lower <- apply(lin_pred_matrix, MARGIN=2, quantile, probs = 0.025)
lin_pred_upper <- apply(lin_pred_matrix, MARGIN=2, quantile, probs = 0.975)

lin_observed <- lin_datA$TotalCatch1
linN <- length(lin_observed)

lindf_plot <- data.frame(
  pred_mean = lin_pred_mean,
  lower = lin_pred_lower,
  upper = lin_pred_upper
)%>%cbind(lin_datA)

ggplot(lindf_plot, aes(x = DOY, y = pred_mean, group = Year, col = as.factor(Year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_point(aes(y = TotalCatch1), size = 1.5) +  # points for observed data
  facet_wrap(~Plot.ID) +
  theme_classic() +labs(
    title = "Linyphiidae",
    y = "Predicted Abudance",
    color = "Year")

#extract peak DOY

lindraws_doy=lin_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(lin_datA$Year))

require(posterior)

linpeak_df <- as_draws_df(lindraws_doy)  %>%
  mutate(.draw = row_number())  %>%
  tidyr::pivot_longer(
    cols = starts_with("DOY_peak_unscaled"),
    names_pattern = "DOY_peak_unscaled\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak") %>%
  mutate(
    year_index = as.integer(year_index),
    year = years[year_index])

# peak timing summary table
linsummary_peak <- linpeak_df  %>%
  group_by(year)  %>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%filter(!(mean==150 | lower==150| upper==150))

print(linsummary_peak)

#posterior preds of phenological curves

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

lin_res=as_draws_df(lin_mod$draws())
linyr_lvls=sort(unique(lin_datA$Year))
names(linyr_lvls)=1:length(linyr_lvls)

linalpha_mean=lin_res%>%summarise(alpha=mean(alpha))%>%as.numeric(alpha)
linbeta_DOYs_mean=lin_res%>%summarise(across(starts_with("beta_DOYs["), mean))%>%unlist()
linbeta_DOYsqs_mean=lin_res%>%summarise(across(starts_with("beta_DOYsqs["), mean))%>%unlist()

DOY_mean=mean(arth_df$DOY)
DOY_sd=sd(arth_df$DOY)
DOY_seq=seq(150, 270, by=1)
DOY_std=(DOY_seq-DOY_mean)/DOY_sd
DOY_sq_std=DOY_std^2
n_DOY=length(DOY_std)
linnyr=length(linbeta_DOYs_mean)

linfitted_curves=generate_fitted_curves(linnyr, linalpha_mean, linbeta_DOYs_mean,
                                        linbeta_DOYsqs_mean,DOY_std, DOY_sq_std, DOY_seq, linyr_lvls)

incyears <- sort(unique(linsummary_peak$year))

linfitted_df=do.call(rbind, linfitted_curves)%>%filter(year %in% incyears)%>%
  mutate(year=as.factor(year))

linp=ggplot(linfitted_df, aes(x=DOY, y=prob, col=year))+
  geom_line(linewidth=0.6, alpha=2)+theme_classic()+
  labs(x="DOY", y="P(abundance)", title="Linyphiidae")+
  scale_color_viridis_d()+
  xlim(150,270)

#sciaridae####
unique(arth_df$HoyeTaxon)

sci_datA <- arth_df %>%filter(HoyeTaxon=="Sciaridae")

sci_data <- list(
  N = nrow(sci_datA),
  y = sci_datA$TotalCatch1,
  obs_days=sci_datA$TrapDays,
  Nyr = length(unique(sci_datA$Year)),
  year_id = as.integer(factor(sci_datA$Year)),
  DOYs=sci_datA$DOYs,
  DOYsqs=sci_datA$DOYsqs,
  yearc=sci_datA$yearc,
  Nplots = length(unique(sci_datA$plot_id)),
  plot_id = sci_datA$plot_id,
  DOY_sd= sd(arth_df$DOY),
  DOY_mean=mean(arth_df$DOY))

#compile model
arth_mod=cmdstan_model("arthropod_phen.stan")

#fit model
sci_mod <- arth_mod$sample(
  data = sci_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500)

sci_draws <- sci_mod$draws(format="df", variables="y_pred")
sci_pred_matrix <- as.data.frame(sci_draws)
sci_pred_matrix=sci_pred_matrix[,-1152:-1154]

sci_pred_mean <- apply(sci_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
sci_pred_lower <- apply(sci_pred_matrix, MARGIN=2, quantile, probs = 0.025)
sci_pred_upper <- apply(sci_pred_matrix, MARGIN=2, quantile, probs = 0.975)

sci_observed <- sci_datA$TotalCatch1
sciN <- length(sci_observed)

scidf_plot <- data.frame(
  pred_mean = sci_pred_mean,
  lower = sci_pred_lower,
  upper = sci_pred_upper
)%>%cbind(sci_datA)

ggplot(scidf_plot, aes(x = DOY, y = pred_mean, group = Year, col = as.factor(Year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_point(aes(y = TotalCatch1), size = 1.5) +  # points for observed data
  facet_wrap(~Plot.ID) +
  theme_classic() +labs(
    title = "Sciaridae",
    y = "Predicted Abudance",
    color = "Year")

#extract peak DOY

scidraws_doy=sci_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(sci_datA$Year))

require(posterior)

scipeak_df <- as_draws_df(scidraws_doy)  %>%
  mutate(.draw = row_number())  %>%
  tidyr::pivot_longer(
    cols = starts_with("DOY_peak_unscaled"),
    names_pattern = "DOY_peak_unscaled\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak") %>%
  mutate(
    year_index = as.integer(year_index),
    year = years[year_index])

# peak timing summary table
scisummary_peak <- scipeak_df  %>%
  group_by(year)  %>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%filter(!(mean==150 | lower==150| upper==150))

print(scisummary_peak)

#posterior preds of phenological curves

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

sci_res=as_draws_df(sci_mod$draws())
sciyr_lvls=sort(unique(sci_datA$Year))
names(sciyr_lvls)=1:length(sciyr_lvls)

scialpha_mean=sci_res%>%summarise(alpha=mean(alpha))%>%as.numeric(alpha)
scibeta_DOYs_mean=sci_res%>%summarise(across(starts_with("beta_DOYs["), mean))%>%unlist()
scibeta_DOYsqs_mean=sci_res%>%summarise(across(starts_with("beta_DOYsqs["), mean))%>%unlist()

DOY_mean=mean(arth_df$DOY)
DOY_sd=sd(arth_df$DOY)
DOY_seq=seq(150, 270, by=1)
DOY_std=(DOY_seq-DOY_mean)/DOY_sd
DOY_sq_std=DOY_std^2
n_DOY=length(DOY_std)
scinyr=length(scibeta_DOYs_mean)

scifitted_curves=generate_fitted_curves(scinyr, scialpha_mean, scibeta_DOYs_mean,
                                        scibeta_DOYsqs_mean,DOY_std, DOY_sq_std, DOY_seq, sciyr_lvls)

incyears <- sort(unique(scisummary_peak$year))

scifitted_df=do.call(rbind, scifitted_curves)%>%filter(year %in% incyears)%>%
  mutate(year=as.factor(year))

scip=ggplot(scifitted_df, aes(x=DOY, y=prob, col=year))+
  geom_line(linewidth=0.6, alpha=2)+theme_classic()+
  labs(x="DOY", y="P(abundance)", title="Sciaridae")+
  scale_color_viridis_d()+
  xlim(150,270)

#coccoidea####

unique(arth_df$HoyeTaxon)

coc_datA <- arth_df %>%filter(HoyeTaxon=="Coccoidea")%>%
  mutate(plot_id = dense_rank(Plot.ID)) %>%
  ungroup()

coc_data <- list(
  N = nrow(coc_datA),
  y = coc_datA$TotalCatch1,
  obs_days=coc_datA$TrapDays,
  Nyr = length(unique(coc_datA$Year)),
  year_id = as.integer(factor(coc_datA$Year)),
  DOYs=coc_datA$DOYs,
  DOYsqs=coc_datA$DOYsqs,
  yearc=coc_datA$yearc,
  Nplots = length(unique(coc_datA$plot_id)),
  plot_id = coc_datA$plot_id,
  DOY_sd= sd(arth_df$DOY),
  DOY_mean=mean(arth_df$DOY))

#compile model
arth_mod=cmdstan_model("arthropod_phen.stan")

#fit model
coc_mod <- arth_mod$sample(
  data = coc_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500)

coc_draws <- coc_mod$draws(format="df", variables="y_pred")
coc_pred_matrix <- as.data.frame(coc_draws)
coc_pred_matrix=coc_pred_matrix[,-389:-391]

coc_pred_mean <- apply(coc_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
coc_pred_lower <- apply(coc_pred_matrix, MARGIN=2, quantile, probs = 0.025)
coc_pred_upper <- apply(coc_pred_matrix, MARGIN=2, quantile, probs = 0.975)

coc_observed <- coc_datA$TotalCatch1
cocN <- length(coc_observed)

cocdf_plot <- data.frame(
  pred_mean = coc_pred_mean,
  lower = coc_pred_lower,
  upper = coc_pred_upper
)%>%cbind(coc_datA)

ggplot(cocdf_plot, aes(x = DOY, y = pred_mean, group = Year, col = as.factor(Year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_point(aes(y = TotalCatch1), size = 1.5) +  # points for observed data
  facet_wrap(~Plot.ID) +
  theme_classic() +labs(
    title = "Coccoidea",
    y = "Predicted Abudance",
    color = "Year")

#extract peak DOY

cocdraws_doy=coc_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(coc_datA$Year))

require(posterior)

cocpeak_df <- as_draws_df(cocdraws_doy)  %>%
  mutate(.draw = row_number())  %>%
  tidyr::pivot_longer(
    cols = starts_with("DOY_peak_unscaled"),
    names_pattern = "DOY_peak_unscaled\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak") %>%
  mutate(
    year_index = as.integer(year_index),
    year = years[year_index])

# peak timing summary table
cocsummary_peak <- cocpeak_df  %>%
  group_by(year)  %>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%filter(!(mean==150 | lower==150| upper==150))

print(cocsummary_peak)

#posterior preds of phenological curves

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

coc_res=as_draws_df(coc_mod$draws())
cocyr_lvls=sort(unique(coc_datA$Year))
names(cocyr_lvls)=1:length(cocyr_lvls)

cocalpha_mean=coc_res%>%summarise(alpha=mean(alpha))%>%as.numeric(alpha)
cocbeta_DOYs_mean=coc_res%>%summarise(across(starts_with("beta_DOYs["), mean))%>%unlist()
cocbeta_DOYsqs_mean=coc_res%>%summarise(across(starts_with("beta_DOYsqs["), mean))%>%unlist()

DOY_mean=mean(arth_df$DOY)
DOY_sd=sd(arth_df$DOY)
DOY_seq=seq(150, 270, by=1)
DOY_std=(DOY_seq-DOY_mean)/DOY_sd
DOY_sq_std=DOY_std^2
n_DOY=length(DOY_std)
cocnyr=length(cocbeta_DOYs_mean)

cocfitted_curves=generate_fitted_curves(cocnyr, cocalpha_mean, cocbeta_DOYs_mean,
                                        cocbeta_DOYsqs_mean,DOY_std, DOY_sq_std, DOY_seq, cocyr_lvls)

incyears <- sort(unique(cocsummary_peak$year))

cocfitted_df=do.call(rbind, cocfitted_curves)%>%filter(year %in% incyears)%>%
  mutate(year=as.factor(year))

cocp=ggplot(cocfitted_df, aes(x=DOY, y=prob, col=year))+
  geom_line(linewidth=0.6, alpha=2)+theme_classic()+
  labs(x="DOY", y="P(abundance)", title="Coccoidea")+
  scale_color_viridis_d()+
  xlim(150,270)

#nymphalidae####
unique(arth_df$HoyeTaxon)

nym_datA <- arth_df %>%filter(HoyeTaxon=="Nymphalidae")%>%
  mutate(plot_id = dense_rank(Plot.ID)) %>%
  ungroup()

nym_data <- list(
  N = nrow(nym_datA),
  y = nym_datA$TotalCatch1,
  obs_days=nym_datA$TrapDays,
  Nyr = length(unique(nym_datA$Year)),
  year_id = as.integer(factor(nym_datA$Year)),
  DOYs=nym_datA$DOYs,
  DOYsqs=nym_datA$DOYsqs,
  yearc=nym_datA$yearc,
  Nplots = length(unique(nym_datA$plot_id)),
  plot_id = nym_datA$plot_id,
  DOY_sd= sd(arth_df$DOY),
  DOY_mean=mean(arth_df$DOY))

#compile model
arth_mod=cmdstan_model("arthropod_phen.stan")

#fit model
nym_mod <- arth_mod$sample(
  data = nym_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500)

nym_draws <- nym_mod$draws(format="df", variables="y_pred")
nym_pred_matrix <- as.data.frame(nym_draws)
nym_pred_matrix=nym_pred_matrix[,-162:-164]

nym_pred_mean <- apply(nym_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
nym_pred_lower <- apply(nym_pred_matrix, MARGIN=2, quantile, probs = 0.025)
nym_pred_upper <- apply(nym_pred_matrix, MARGIN=2, quantile, probs = 0.975)

nym_observed <- nym_datA$TotalCatch1
nymN <- length(nym_observed)

nymdf_plot <- data.frame(
  pred_mean = nym_pred_mean,
  lower = nym_pred_lower,
  upper = nym_pred_upper
)%>%cbind(nym_datA)

ggplot(nymdf_plot, aes(x = DOY, y = pred_mean, group = Year, col = as.factor(Year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_point(aes(y = TotalCatch1), size = 1.5) +  # points for observed data
  facet_wrap(~Plot.ID) +
  theme_classic() +labs(
    title = "Nymphalidae",
    y = "Predicted Abudance",
    color = "Year")

#extract peak DOY

nymdraws_doy=nym_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(nym_datA$Year))

require(posterior)

nympeak_df <- as_draws_df(nymdraws_doy)  %>%
  mutate(.draw = row_number())  %>%
  tidyr::pivot_longer(
    cols = starts_with("DOY_peak_unscaled"),
    names_pattern = "DOY_peak_unscaled\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak") %>%
  mutate(
    year_index = as.integer(year_index),
    year = years[year_index])

# peak timing summary table
nymsummary_peak <- nympeak_df  %>%
  group_by(year)  %>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%filter(!(mean==150 | lower==150| upper==150))

print(nymsummary_peak)

#posterior preds of phenological curves

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

nym_res=as_draws_df(nym_mod$draws())
nymyr_lvls=sort(unique(nym_datA$Year))
names(nymyr_lvls)=1:length(nymyr_lvls)

nymalpha_mean=nym_res%>%summarise(alpha=mean(alpha))%>%as.numeric(alpha)
nymbeta_DOYs_mean=nym_res%>%summarise(across(starts_with("beta_DOYs["), mean))%>%unlist()
nymbeta_DOYsqs_mean=nym_res%>%summarise(across(starts_with("beta_DOYsqs["), mean))%>%unlist()

DOY_mean=mean(arth_df$DOY)
DOY_sd=sd(arth_df$DOY)
DOY_seq=seq(150, 270, by=1)
DOY_std=(DOY_seq-DOY_mean)/DOY_sd
DOY_sq_std=DOY_std^2
n_DOY=length(DOY_std)
nymnyr=length(nymbeta_DOYs_mean)

nymfitted_curves=generate_fitted_curves(nymnyr, nymalpha_mean, nymbeta_DOYs_mean,
                                        nymbeta_DOYsqs_mean,DOY_std, DOY_sq_std, DOY_seq, nymyr_lvls)

incyears <- sort(unique(nymsummary_peak$year))

nymfitted_df=do.call(rbind, nymfitted_curves)%>%filter(year %in% incyears)%>%
  mutate(year=as.factor(year))

nymp=ggplot(nymfitted_df, aes(x=DOY, y=prob, col=year))+
  geom_line(linewidth=0.6, alpha=2)+theme_classic()+
  labs(x="DOY", y="P(abundance)", title="Nymphalidae")+
  scale_color_viridis_d()+
  xlim(150,270)

#phoridae####
unique(arth_df$HoyeTaxon)

pho_datA <- arth_df %>%filter(HoyeTaxon=="Phoridae")%>%
  mutate(plot_id = dense_rank(Plot.ID)) %>%
  ungroup()

pho_data <- list(
  N = nrow(pho_datA),
  y = pho_datA$TotalCatch1,
  obs_days=pho_datA$TrapDays,
  Nyr = length(unique(pho_datA$Year)),
  year_id = as.integer(factor(pho_datA$Year)),
  DOYs=pho_datA$DOYs,
  DOYsqs=pho_datA$DOYsqs,
  yearc=pho_datA$yearc,
  Nplots = length(unique(pho_datA$plot_id)),
  plot_id = pho_datA$plot_id,
  DOY_sd= sd(arth_df$DOY),
  DOY_mean=mean(arth_df$DOY))

#compile model
arth_mod=cmdstan_model("arthropod_phen.stan")

#fit model
pho_mod <- arth_mod$sample(
  data = pho_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500)

pho_draws <- pho_mod$draws(format="df", variables="y_pred")
pho_pred_matrix <- as.data.frame(pho_draws)
pho_pred_matrix=pho_pred_matrix[,-605:-607]

pho_pred_mean <- apply(pho_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
pho_pred_lower <- apply(pho_pred_matrix, MARGIN=2, quantile, probs = 0.025)
pho_pred_upper <- apply(pho_pred_matrix, MARGIN=2, quantile, probs = 0.975)

pho_observed <- pho_datA$TotalCatch1
phoN <- length(pho_observed)

phodf_plot <- data.frame(
  pred_mean = pho_pred_mean,
  lower = pho_pred_lower,
  upper = pho_pred_upper
)%>%cbind(pho_datA)

ggplot(phodf_plot, aes(x = DOY, y = pred_mean, group = Year, col = as.factor(Year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_point(aes(y = TotalCatch1), size = 1.5) +  # points for observed data
  facet_wrap(~Plot.ID) +
  theme_classic() +labs(
    title = "Phoridae",
    y = "Predicted Abudance",
    color = "Year")

#extract peak DOY

phodraws_doy=pho_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(pho_datA$Year))

require(posterior)

phopeak_df <- as_draws_df(phodraws_doy)  %>%
  mutate(.draw = row_number())  %>%
  tidyr::pivot_longer(
    cols = starts_with("DOY_peak_unscaled"),
    names_pattern = "DOY_peak_unscaled\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak") %>%
  mutate(
    year_index = as.integer(year_index),
    year = years[year_index])

# peak timing summary table
phosummary_peak <- phopeak_df  %>%
  group_by(year)  %>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%filter(!(mean==150 | lower==150| upper==150))

print(phosummary_peak)

#posterior preds of phenological curves

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

pho_res=as_draws_df(pho_mod$draws())
phoyr_lvls=sort(unique(pho_datA$Year))
names(phoyr_lvls)=1:length(phoyr_lvls)

phoalpha_mean=pho_res%>%summarise(alpha=mean(alpha))%>%as.numeric(alpha)
phobeta_DOYs_mean=pho_res%>%summarise(across(starts_with("beta_DOYs["), mean))%>%unlist()
phobeta_DOYsqs_mean=pho_res%>%summarise(across(starts_with("beta_DOYsqs["), mean))%>%unlist()

DOY_mean=mean(arth_df$DOY)
DOY_sd=sd(arth_df$DOY)
DOY_seq=seq(150, 270, by=1)
DOY_std=(DOY_seq-DOY_mean)/DOY_sd
DOY_sq_std=DOY_std^2
n_DOY=length(DOY_std)
phonyr=length(phobeta_DOYs_mean)

phofitted_curves=generate_fitted_curves(phonyr, phoalpha_mean, phobeta_DOYs_mean,
                                        phobeta_DOYsqs_mean,DOY_std, DOY_sq_std, DOY_seq, phoyr_lvls)

incyears <- sort(unique(phosummary_peak$year))

phofitted_df=do.call(rbind, phofitted_curves)%>%filter(year %in% incyears)%>%
  mutate(year=as.factor(year))

phop=ggplot(phofitted_df, aes(x=DOY, y=prob, col=year))+
  geom_line(linewidth=0.6, alpha=2)+theme_classic()+
  labs(x="DOY", y="P(abundance)", title="Phoridae")+
  scale_color_viridis_d()+
  xlim(150,270)

#Acari####
unique(arth_df$HoyeTaxon)

aca_datA <- arth_df %>%filter(HoyeTaxon=="Acari")%>%
  mutate(plot_id = dense_rank(Plot.ID)) %>%
  ungroup()

aca_data <- list(
  N = nrow(aca_datA),
  y = aca_datA$TotalCatch1,
  obs_days=aca_datA$TrapDays,
  Nyr = length(unique(aca_datA$Year)),
  year_id = as.integer(factor(aca_datA$Year)),
  DOYs=aca_datA$DOYs,
  DOYsqs=aca_datA$DOYsqs,
  yearc=aca_datA$yearc,
  Nplots = length(unique(aca_datA$plot_id)),
  plot_id = aca_datA$plot_id,
  DOY_sd= sd(arth_df$DOY),
  DOY_mean=mean(arth_df$DOY))

#compile model
arth_mod=cmdstan_model("arthropod_phen.stan")

#fit model
aca_mod <- arth_mod$sample(
  data = aca_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500)

aca_draws <- aca_mod$draws(format="df", variables="y_pred")
aca_pred_matrix <- as.data.frame(aca_draws)
aca_pred_matrix=aca_pred_matrix[,-1436:-1438]

aca_pred_mean <- apply(aca_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
aca_pred_lower <- apply(aca_pred_matrix, MARGIN=2, quantile, probs = 0.025)
aca_pred_upper <- apply(aca_pred_matrix, MARGIN=2, quantile, probs = 0.975)

aca_observed <- aca_datA$TotalCatch1
acaN <- length(aca_observed)

acadf_plot <- data.frame(
  pred_mean = aca_pred_mean,
  lower = aca_pred_lower,
  upper = aca_pred_upper
)%>%cbind(aca_datA)

ggplot(acadf_plot, aes(x = DOY, y = pred_mean, group = Year, col = as.factor(Year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_point(aes(y = TotalCatch1), size = 1.5) +  # points for observed data
  facet_wrap(~Plot.ID) +
  theme_classic() +labs(
    title = "Acari",
    y = "Predicted Abudance",
    color = "Year")

#extract peak DOY

acadraws_doy=aca_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(aca_datA$Year))

require(posterior)

acapeak_df <- as_draws_df(acadraws_doy)  %>%
  mutate(.draw = row_number())  %>%
  tidyr::pivot_longer(
    cols = starts_with("DOY_peak_unscaled"),
    names_pattern = "DOY_peak_unscaled\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak") %>%
  mutate(
    year_index = as.integer(year_index),
    year = years[year_index])

# peak timing summary table
acasummary_peak <- acapeak_df  %>%
  group_by(year)  %>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%filter(!(mean==150 | lower==150| upper==150))

print(acasummary_peak)

#posterior preds of phenological curves

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

aca_res=as_draws_df(aca_mod$draws())
acayr_lvls=sort(unique(aca_datA$Year))
names(acayr_lvls)=1:length(acayr_lvls)

acaalpha_mean=aca_res%>%summarise(alpha=mean(alpha))%>%as.numeric(alpha)
acabeta_DOYs_mean=aca_res%>%summarise(across(starts_with("beta_DOYs["), mean))%>%unlist()
acabeta_DOYsqs_mean=aca_res%>%summarise(across(starts_with("beta_DOYsqs["), mean))%>%unlist()

DOY_mean=mean(arth_df$DOY)
DOY_sd=sd(arth_df$DOY)
DOY_seq=seq(150, 270, by=1)
DOY_std=(DOY_seq-DOY_mean)/DOY_sd
DOY_sq_std=DOY_std^2
n_DOY=length(DOY_std)
acanyr=length(acabeta_DOYs_mean)

acafitted_curves=generate_fitted_curves(acanyr, acaalpha_mean, acabeta_DOYs_mean,
                                        acabeta_DOYsqs_mean,DOY_std, DOY_sq_std, DOY_seq, acayr_lvls)

incyears <- sort(unique(acasummary_peak$year))

acafitted_df=do.call(rbind, acafitted_curves)%>%filter(year %in% incyears)%>%
  mutate(year=as.factor(year))

acap=ggplot(acafitted_df, aes(x=DOY, y=prob, col=year))+
  geom_line(linewidth=0.6, alpha=2)+theme_classic()+
  labs(x="DOY", y="P(abundance)", title="Acari")+
  scale_color_viridis_d()+
  xlim(150,270)

#Lycosidae####
lyc_datA <- arth_df %>%filter(HoyeTaxon=="Lycosidae")%>%
  mutate(plot_id = dense_rank(Plot.ID)) %>%
  ungroup()

lyc_data <- list(
  N = nrow(lyc_datA),
  y = lyc_datA$TotalCatch1,
  obs_days=lyc_datA$TrapDays,
  Nyr = length(unique(lyc_datA$Year)),
  year_id = as.integer(factor(lyc_datA$Year)),
  DOYs=lyc_datA$DOYs,
  DOYsqs=lyc_datA$DOYsqs,
  yearc=lyc_datA$yearc,
  Nplots = length(unique(lyc_datA$plot_id)),
  plot_id = lyc_datA$plot_id,
  DOY_sd= sd(arth_df$DOY),
  DOY_mean=mean(arth_df$DOY))

#compile model
arth_mod=cmdstan_model("arthropod_phen.stan")

#fit model
lyc_mod <- arth_mod$sample(
  data = lyc_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500)

lyc_draws <- lyc_mod$draws(format="df", variables="y_pred")
lyc_pred_matrix <- as.data.frame(lyc_draws)
lyc_pred_matrix=lyc_pred_matrix[,-1348:-1350]

lyc_pred_mean <- apply(lyc_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
lyc_pred_lower <- apply(lyc_pred_matrix, MARGIN=2, quantile, probs = 0.025)
lyc_pred_upper <- apply(lyc_pred_matrix, MARGIN=2, quantile, probs = 0.975)

lyc_observed <- lyc_datA$TotalCatch1
lycN <- length(lyc_observed)

lycdf_plot <- data.frame(
  pred_mean = lyc_pred_mean,
  lower = lyc_pred_lower,
  upper = lyc_pred_upper
)%>%cbind(lyc_datA)

ggplot(lycdf_plot, aes(x = DOY, y = pred_mean, group = Year, col = as.factor(Year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_point(aes(y = TotalCatch1), size = 1.5) +  # points for observed data
  facet_wrap(~Plot.ID) +
  theme_classic() +labs(
    title = "Lycosidae",
    y = "Predicted Abudance",
    color = "Year")

#extract peak DOY

lycdraws_doy=lyc_mod$draws(variable="DOY_peak_unscaled", format="df")

years <- sort(unique(lyc_datA$Year))

require(posterior)

lycpeak_df <- as_draws_df(lycdraws_doy)  %>%
  mutate(.draw = row_number())  %>%
  tidyr::pivot_longer(
    cols = starts_with("DOY_peak_unscaled"),
    names_pattern = "DOY_peak_unscaled\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak") %>%
  mutate(
    year_index = as.integer(year_index),
    year = years[year_index])

# peak timing summary table
lycsummary_peak <- lycpeak_df  %>%
  group_by(year)  %>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95)
  )%>%filter(!(mean==150 | lower==150| upper==150))

print(lycsummary_peak)

#posterior preds of phenological curves

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

lyc_res=as_draws_df(lyc_mod$draws())
lycyr_lvls=sort(unique(lyc_datA$Year))
names(lycyr_lvls)=1:length(lycyr_lvls)

lycalpha_mean=lyc_res%>%summarise(alpha=mean(alpha))%>%as.numeric(alpha)
lycbeta_DOYs_mean=lyc_res%>%summarise(across(starts_with("beta_DOYs["), mean))%>%unlist()
lycbeta_DOYsqs_mean=lyc_res%>%summarise(across(starts_with("beta_DOYsqs["), mean))%>%unlist()

DOY_mean=mean(arth_df$DOY)
DOY_sd=sd(arth_df$DOY)
DOY_seq=seq(150, 270, by=1)
DOY_std=(DOY_seq-DOY_mean)/DOY_sd
DOY_sq_std=DOY_std^2
n_DOY=length(DOY_std)
lycnyr=length(lycbeta_DOYs_mean)

lycfitted_curves=generate_fitted_curves(lycnyr, lycalpha_mean, lycbeta_DOYs_mean,
                                        lycbeta_DOYsqs_mean,DOY_std, DOY_sq_std, DOY_seq, lycyr_lvls)

incyears <- sort(unique(lycsummary_peak$year))

lycfitted_df=do.call(rbind, lycfitted_curves)%>%filter(year %in% incyears)%>%
  mutate(year=as.factor(year))

lycp=ggplot(lycfitted_df, aes(x=DOY, y=prob, col=year))+
  geom_line(linewidth=0.6, alpha=2)+theme_classic()+
  labs(x="DOY", y="P(abundance)", title="Lycosidae")+
  scale_color_viridis_d()+
  xlim(150,270)
