#library(cmdstanr)
library(dplyr)
library(ggplot2)
library(ggpubr)


# ===================
# ====== Jarno's tests with an alternative parameterization ========

# Useful links for basic posterior plotting and checking MCMC convergence
# https://mc-stan.org/bayesplot/reference/MCMC-diagnostics.html
# https://mc-stan.org/bayesplot/index.html

library(rstan)
library(bayesplot)
library(gridExtra)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# restructure data
# # Pat's paths
# dat_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\data"
# dat_name=paste(dat_path, '\\plant_datA','.csv', sep = '')

# Jarno's paths
setwd("~/Documents_extended/Projects_data/Zackenberg_phenology/")
dat_path = paste0(getwd(), "/data")
dat_name = paste(dat_path, '/plant_datA','.csv', sep = '')

plant_datA=read.csv(dat_name, header=T, sep=',',  stringsAsFactors = F)

dry_datA <- plant_datA%>%filter(species=="Dryas")

dryas_data <- list(
  N = nrow(dry_datA),
  tot_F = dry_datA$tot_F,
  tot_NF = dry_datA$tot_NF,
  Nyr = length(unique(dry_datA$year)),
  year_id = as.integer(factor(dry_datA$year)),
  DOYs=dry_datA$DOYs,
  DOYsqs=dry_datA$DOYsqs,
  yearc=dry_datA$yearc,
  year_vec = (unique(dry_datA$yearc) - mean(unique(dry_datA$yearc)))/sd(unique(dry_datA$yearc)),
  Nplots = length(unique(dry_datA$plot_id)),
  plot_id = dry_datA$plot_id,
  DOY_sd= sd(plant_datA$DOY),
  DOY_mean=mean(plant_datA$DOY),
  ND=121)

plot(dryas_data$DOYs, dryas_data$tot_F/(dryas_data$tot_F+dryas_data$tot_NF))
# init = list(list(sigma_mu = c(runif(1)) ),
#             list(sigma_mu = c(runif(1)) ),
#             list(sigma_mu = c(runif(1)) ),
#             list(sigma_mu = c(runif(1)) ))

# Fit the model
dry_mod <- stan("plant_phen6.stan",
                 data = dryas_data,
#                 init = init,
                 seed = 123,
                 chains = 4,
                 iter = 300,
                 warmup = 200)

# Check for the convergence and mixing of sampling
# ================================================
# traceplots
traceplot(dry_mod, pars = "mu")
traceplot(dry_mod, pars = "width")
traceplot(dry_mod, pars = "alpha_year")
traceplot(dry_mod, pars = c("sigma_mu","sigma_width"))

# Rhat values
grid.arrange(mcmc_rhat_hist(rhat(dry_mod)), 
             mcmc_rhat(rhat(dry_mod)), 
             mcmc_neff_hist(neff_ratio(dry_mod)), nrow = 2)

# Regular, visual, posterior predictive checks
color_scheme_set("red")
yrep = as.matrix(dry_mod, pars= "y_pred")
yrep <- yrep[(dim(yrep)[1]-99):dim(yrep)[1],]   # replicate data, take only the last 100 samples to speed up drawing
grid.arrange(ppc_dens_overlay(y = dryas_data$tot_F,
                              yrep = yrep),
             ppc_ecdf_overlay(y = dryas_data$tot_F,
                              yrep = yrep),
             ppc_pit_ecdf(y = dryas_data$tot_F, yrep = yrep, prob = 0.99, plot_diff = FALSE),
             ppc_pit_ecdf(y = dryas_data$tot_F, yrep = yrep, prob = 0.99, plot_diff = TRUE),
             nrow = 2)
ppc_hist(y = dryas_data$tot_F, yrep[1:11,])


# Examine posterior distributions of the parameters
# =========================================
plot(dry_mod, pars="mu")     # Peak dates
plot(dry_mod, pars="width")  # width parameter of the phenology



# plot the phenology curve together with data
## ==========================================

# Vector of DoYs for plotting
DOY = seq(-2, 2, length = 100)

# Posterior samples of parameters to calculate the curve
mu = as.matrix(dry_mod, pars = "mu")
width = as.matrix(dry_mod, pars = "width")
alpha_year = as.matrix(dry_mod, pars = "alpha_year")
u_plot = as.matrix(dry_mod, pars = "u_plot")

#eta = array(nrow = dim(mu.samp)[1], ncol=length(DOY), )
eta = vector(length = length(DOY), )
plot.id = unique(dryas_data$plot_id)
cols = rainbow(6)
# Loop through the years
par(mfrow = c(3,4))
invLogit = function(x){exp(x)/(1+exp(x))}
for (y in 1:dim(mu.samp)[2]) {
  # loop through the plots
  for (pl in 1:length(plot.id)){
    # loop through DoYs
    for (nd in 1:length(DOY)){
      eta[nd] = mean( invLogit(alpha_year[,y] - ( DOY[nd] - mu[,y] )^2./width[,y]^2 + u_plot[plot.id[pl]]) )
    }
    if (pl == 1){
      plot(DOY,eta, type = "l", col = cols[pl], 
           ylim = c(0,1),
           title = paste("year", y))
    } else {
      lines(DOY,eta, type = "l", col = cols[pl])
    }
    ind = dryas_data$plot_id==pl & dryas_data$year_id==y
    points(dryas_data$DOYs[ind], dryas_data$tot_F[ind]/(dryas_data$tot_F[ind]+dryas_data$tot_NF[ind]), col=cols[pl])
  }
  
}


####### Model 6b ###################

# Fit the model
dry_mod <- stan("plant_phen6b.stan",
                data = dryas_data,
                #                 init = init,
                seed = 123,
                chains = 4,
                iter = 300,
                warmup = 200)

# Check for the convergence and mixing of sampling
# ================================================
# traceplots
traceplot(dry_mod, pars = "mu")
traceplot(dry_mod, pars = "beta_mu")
traceplot(dry_mod, pars = "width")
traceplot(dry_mod, pars = "alpha_year")
traceplot(dry_mod, pars = c("sigma_mu","sigma_width"))

# Rhat values
grid.arrange(mcmc_rhat_hist(rhat(dry_mod)), 
             mcmc_rhat(rhat(dry_mod)), 
             mcmc_neff_hist(neff_ratio(dry_mod)), nrow = 2)

# Regular, visual, posterior predictive checks
color_scheme_set("red")
yrep = as.matrix(dry_mod, pars= "y_pred")
yrep <- yrep[(dim(yrep)[1]-99):dim(yrep)[1],]   # replicate data, take only the last 100 samples to speed up drawing
grid.arrange(ppc_dens_overlay(y = dryas_data$tot_F,
                              yrep = yrep),
             ppc_ecdf_overlay(y = dryas_data$tot_F,
                              yrep = yrep),
             ppc_pit_ecdf(y = dryas_data$tot_F, yrep = yrep, prob = 0.99, plot_diff = FALSE),
             ppc_pit_ecdf(y = dryas_data$tot_F, yrep = yrep, prob = 0.99, plot_diff = TRUE),
             nrow = 2)
ppc_hist(y = dryas_data$tot_F, yrep[1:11,])


# Examine posterior distributions of the parameters
# =========================================
plot(dry_mod, pars="mu")     # Peak dates
plot(dry_mod, pars="beta_mu")
plot(dry_mod, pars="width")  # width parameter of the phenology



# plot the phenology curve together with data
## ==========================================

# Vector of DoYs for plotting
DOY = seq(-2, 2, length = 100)

# Posterior samples of parameters to calculate the curve
mu = as.matrix(dry_mod, pars = "mu")
width = as.matrix(dry_mod, pars = "width")
alpha_year = as.matrix(dry_mod, pars = "alpha_year")
u_plot = as.matrix(dry_mod, pars = "u_plot")
u_plot_mu = as.matrix(dry_mod, pars = "u_plot_mu")

#eta = array(nrow = dim(mu.samp)[1], ncol=length(DOY), )
eta = vector(length = length(DOY), )
plot.id = unique(dryas_data$plot_id)
cols = rainbow(6)
# Loop through the years
par(mfrow = c(3,4))
invLogit = function(x){exp(x)/(1+exp(x))}
for (y in 1:dim(mu.samp)[2]) {
  # loop through the plots
  for (pl in 1:length(plot.id)){
    # loop through DoYs
    for (nd in 1:length(DOY)){
      eta[nd] = mean( invLogit(alpha_year[,y] - ( DOY[nd] - (mu[,y]+u_plot_mu[plot.id[pl]]) )^2./width[,y]^2 + u_plot[plot.id[pl]]) )
    }
    if (pl == 1){
      plot(DOY,eta, type = "l", col = cols[pl], 
           ylim = c(0,1),
           title = paste("year", y))
    } else {
      lines(DOY,eta, type = "l", col = cols[pl])
    }
    ind = dryas_data$plot_id==pl & dryas_data$year_id==y
    points(dryas_data$DOYs[ind], 
           dryas_data$tot_F[ind]/(dryas_data$tot_F[ind]+dryas_data$tot_NF[ind]), 
           col=cols[pl])
  }
  
}




#########################################################################
# Below this old code 

#########################################################################

#predictions

dry_draws <- dry_mod2$draws(format="df", variables="y_pred")
y_pred_matrix <- as.data.frame(dry_draws)
y_pred_matrix=y_pred_matrix[,-1158:-1160]

y_pred_mean <- apply(y_pred_matrix, MARGIN=2, mean) #margin=2 is column mean
y_pred_lower <- apply(y_pred_matrix, MARGIN=2, quantile, probs = 0.025)
y_pred_upper <- apply(y_pred_matrix, MARGIN=2, quantile, probs = 0.975)


#plot predictions of flower totals

observed <- dry_datA$tot_F
N <- length(observed)

df_plot <- data.frame(
  pred_mean = y_pred_mean,
  lower = y_pred_lower,
  upper = y_pred_upper
)%>%cbind(dry_datA)

#to check weird years
highlight_years <- c( "1998", "2018", "1997", "2014", "2015", "2020", "2021")

ggplot(df_plot, aes(x = DOY, y = pred_mean, group = year, col = as.factor(year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  # geom_line(data = subset(df_plot, year %in% highlight_years),
  #           aes(x = DOY, y = pred_mean, group = year, color = as.factor(year)),
  #           linewidth = 1.2) +  # bold lines for highlighted years
  geom_point(aes(y = tot_F), size = 1.5) +  # points for observed data
  facet_wrap(~Plot) +
  # scale_color_manual(
  #   values = c( "1998" = "orange", "2018"="red", "1997"= "blue", "2014"="green", "2015"="pink", "2020"="black", "2021"="violet"),
  #   breaks = highlight_years,
  #   guide = guide_legend(title = "Odd Years")
  # ) +
  theme_classic() +
  labs(
    title = "Predicted Flowering Curve per Year",
    y = "Predicted Flower Totals",
    color = "Year")

#extract peak DOY
# Extract draws
draws_dry = dry_mod4$draws(variables = c("mu_year"), format="df")

beta_DOYs_dry   = as.matrix(draws_dry[, startsWith(colnames(draws_dry), "beta_DOYs[")])
beta_DOYsqs_dry = as.matrix(draws_dry[, startsWith(colnames(draws_dry), "beta_DOYsqs[")])

# Standardization constants
DOY_sd   = sd(plant_datA$DOY)
DOY_mean = mean(plant_datA$DOY)

n_draws = nrow(beta_DOYs_dry)
Nyr     = ncol(beta_DOYs_dry)

dry_peak = matrix(NA, nrow = n_draws, ncol = Nyr)

for (y in 1:Nyr) {
  DOY_peak_std = -beta_DOYs_dry[, y] / (2 * beta_DOYsqs_dry[, y])
  dry_peak[, y] = DOY_peak_std * DOY_sd + DOY_mean
}

colnames(dry_peak) = unique(dry_datA$year)
dry_peak=as.data.frame(dry_peak)

# peak timing summary table
drysummary_peak <- dry_peak%>%pivot_longer(cols=1:24, names_to="year")%>%
  rename(peak=value)%>%
  group_by(year)%>%
  summarise(
    mean  = mean(peak),
    median= median(peak),
    lower= quantile(peak, 0.05),
    upper= quantile(peak, 0.95))

print(drysummary_peak)

dryas_peak=ggplot(drysummary_peak, aes(x = year, y = median, col=as.factor(year))) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
  labs(
    title = "Dryas peak timing",
    x = "Year",
    y = "Peak Day of Year (DOY)",
    caption = "Median and 90% CI")+
  theme_classic() #1998 uncertainty (coldest year and longest flowering duration): https://www.nature.com/articles/nclimate1909

slope_samples <- peak_df |>
  group_by(.draw) |>
  summarise(
    slope = coef(lm(DOY_peak ~ year))[2]  # Extract the slope
  )

slope_summary <- slope_samples |>
  summarise(
    mean    = mean(slope),
    median  = median(slope),
    lower90 = quantile(slope, 0.05),
    upper90 = quantile(slope, 0.95),
    lower95 = quantile(slope, 0.025),
    upper95 = quantile(slope, 0.975)
  )

print(slope_summary)

library(ggplot2)

dryas_slope=ggplot(slope_samples, aes(x = slope)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "slope distribution",
    x = "Change in Peak DOY per Year",
    y = "Density"
  ) +
  theme_classic()

#curvatures####
draws_ct=dry_mod$draws(variable="beta_DOYsqs", format="df")

ct_df <- as_draws_df(draws_ct) |>
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

summary_ct <- ct_df |>
  group_by(year) |>
  summarise(
    mean  = mean(ct),
    median= median(ct),
    lower= quantile(ct, 0.05),
    upper= quantile(ct, 0.95)
  )%>%
  mutate(
    highlight_group = ifelse(year %in% highlight_years, as.character(year), "Other"))

dryas_curve=ggplot(summary_ct, aes(x = year, y = mean, col=as.factor(year))) +
  geom_line(data = subset(summary_peak, year %in% highlight_years),
            aes(x = year, y = mean, group = year, color = as.factor(year)),
            linewidth = 1.2) +  # bold lines for highlighted years
  geom_point(size = 2) +
  # stat_smooth(aes(x = year, y = mean, group = 1), method = "lm",
  #             color = "black", se = TRUE, linewidth = 1) +
 geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
  scale_color_manual(
    values = c( "1998" = "orange", "2018"="red"),
    breaks = highlight_years,
    guide = guide_legend(title = "Odd Years")
  ) +
  labs(
    title = "curvature",
    x = "Year",
    y = "curvature (beta_DOYsqs)",
    caption = "Mean and 90% CI")+
  theme_classic()+ylim(min(summary_ct$lower), max(summary_ct$upper))

dry_res=ggarrange(dryas_peak, dryas_slope, dryas_curve, nrow=1, ncol=3)
annotate_figure(dry_res, top = text_grob("Dryas", face = "bold", size = 20))

#posterior preds of phenological curves

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

dry_res=as_draws_df(dry_mod$draws())
dryyr_lvls=sort(unique(dry_datA$year))
names(dryyr_lvls)=1:length(dryyr_lvls)

dryalpha_mean=dry_res%>%summarise(alpha=mean(alpha))%>%as.numeric(alpha)
drybeta_DOYs_mean=dry_res%>%summarise(across(starts_with("beta_DOYs["), mean))%>%unlist()
drybeta_DOYsqs_mean=dry_res%>%summarise(across(starts_with("beta_DOYsqs["), mean))%>%unlist()

DOY_mean=mean(plant_datA$DOY)
DOY_sd=sd(plant_datA$DOY)
DOY_seq=seq(150, 270, by=1)
DOY_std=(DOY_seq-DOY_mean)/DOY_sd
DOY_sq_std=DOY_std^2
n_DOY=length(DOY_std)
drynyr=length(drybeta_DOYs_mean)

dryfitted_curves=generate_fitted_curves(drynyr, dryalpha_mean, drybeta_DOYs_mean,
                                        drybeta_DOYsqs_mean,DOY_std, DOY_sq_std, DOY_seq, dryyr_lvls)

incyears <- sort(unique(summary_peak$year))

dryfitted_df=do.call(rbind, dryfitted_curves)%>%filter(year %in% incyears)%>%
  mutate(year=as.factor(year))


dryp=ggplot(dryfitted_df, aes(x=DOY, y=prob, col=year))+
  geom_line(linewidth=0.6, alpha=2)+theme_classic()+
  labs(x="DOY", y="P(flower)", title="Dryas")+
  scale_color_viridis_d()+
  xlim(150,270)
