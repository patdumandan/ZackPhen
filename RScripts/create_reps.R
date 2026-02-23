#call model output names
mod_list=list.files(path="/scratch/project_2017453/spp_models",
                    pattern = "\\.rds$", full.names=TRUE)

arth_list=list.files(path="/scratch/project_2017453/arth_models",
                    pattern = "\\.rds$", full.names=TRUE)

#extract relevant global trend param
dryas_betamu <- extract_beta_mu_win(mod_list, species_name="Dryas")
silene_betamu <- extract_beta_mu_win(mod_list, species_name="Silene")
papaver_betamu <- extract_beta_mu_win(mod_list, species_name="Papaver")
salix_betamu <- extract_beta_mu_win(mod_list, species_name="Salix")
cassiope_betamu <- extract_beta_mu_win(mod_list, species_name="Cassiope")
saxifraga_betamu <- extract_beta_mu_win(mod_list, species_name="Saxifraga")

muscidae_betamu <- extract_beta_mu_win(arth_list, species_name="Muscidae")
chironomidae_betamu <- extract_beta_mu_win(arth_list, species_name="Chironomidae")
ichneumonidae_betamu <- extract_beta_mu_win(arth_list, species_name="Ichneumonidae")
sciaridae_betamu <- extract_beta_mu_win(arth_list, species_name="Sciaridae")
lycosidae_betamu <- extract_beta_mu_win(arth_list, species_name="Lycosidae")
nymphalidae_betamu <- extract_beta_mu_win(arth_list, species_name="Nymphalidae")
phoridae_betamu <- extract_beta_mu_win(arth_list, species_name="Phoridae")
acari_betamu <- extract_beta_mu_win(arth_list, species_name="Acari")
collembola_betamu <- extract_beta_mu_win(arth_list, species_name="Collembola")
linyphiidae_betamu <- extract_beta_mu_win(arth_list, species_name="Linyphiidae")
coccoidea_betamu <- extract_beta_mu_win(arth_list, species_name="Coccoidea")

mus_betamu <- lapply(muscidae_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

sci_betamu <- lapply(sciaridae_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

chi_betamu <- lapply(chironomidae_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

ich_betamu <- lapply(ichneumonidae_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

lyc_betamu <- lapply(lycosidae_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

nym_betamu <- lapply(nymphalidae_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

pho_betamu <- lapply(phoridae_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

aca_betamu <- lapply(acari_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

col_betamu <- lapply(collembola_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

lin_betamu <- lapply(linyphiidae_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

coc_betamu <- lapply(coccoidea_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

dry_betamu <- lapply(dryas_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

sax_betamu <- lapply(saxifraga_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

cas_betamu <- lapply(cassiope_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

sal_betamu <- lapply(salix_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})


sil_betamu <- lapply(silene_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

pap_betamu <- lapply(papaver_betamu, function(x) {
  data.frame(
    draw     = seq_along(x$beta_mu),
    beta_mu  = x$beta_mu,
    start_yr = x$start_yr,
    end_yr   = x$end_yr,
    tsl      = (x$end_yr-x$start_yr)+1,
    species  =x$species)
})

mus_betamu <- do.call(rbind, mus_betamu)
chi_betamu <- do.call(rbind, chi_betamu)
ich_betamu <- do.call(rbind, ich_betamu)
sci_betamu <- do.call(rbind, sci_betamu)
coc_betamu <- do.call(rbind, coc_betamu)
lin_betamu <- do.call(rbind, lin_betamu)
col_betamu <- do.call(rbind, col_betamu)
aca_betamu <- do.call(rbind, aca_betamu)
pho_betamu <- do.call(rbind, pho_betamu)
nym_betamu <- do.call(rbind, nym_betamu)
lyc_betamu <- do.call(rbind, lyc_betamu)

dry_betamu <- do.call(rbind, dry_betamu)
sil_betamu <- do.call(rbind, sil_betamu)
pap_betamu <- do.call(rbind, pap_betamu)
cas_betamu <- do.call(rbind, cas_betamu)
sal_betamu <- do.call(rbind, sal_betamu)
sax_betamu <- do.call(rbind, sax_betamu)

chironomidae_slope_df <- chi_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

sciaridae_slope_df <- sci_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

ichneumonidae_slope_df <- ich_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

muscidae_slope_df <- mus_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

dryas_slope_df <- dry_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

salix_slope_df <- sal_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

saxifraga_slope_df <- sax_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

silene_slope_df <- sil_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

cassiope_slope_df <- cas_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

papaver_slope_df <- pap_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

ichp=ggplot(ichneumonidae_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Ichneumonidae")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

scip=ggplot(sciaridae_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Sciaridae")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))
phoridae_slope_df <- pho_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

phop=ggplot(phoridae_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Phoridae")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

nymphalidae_slope_df <- nym_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

nymp=ggplot(nymphalidae_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Nymphalidae")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

lycosidae_slope_df <- lyc_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

lycp=ggplot(lycosidae_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Lycosidae")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

linyphiidae_slope_df <- lin_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

linp=ggplot(linyphiidae_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Linyphiidae")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

collembola_slope_df <- col_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

colp=ggplot(collembola_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Collembola")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

coccoidea_slope_df <- coc_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

cocp=ggplot(coccoidea_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Coccoidea")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

acari_slope_df <- aca_betamu %>%
  group_by(start_yr, end_yr, tsl) %>%
  summarise(
    slope_median = median(beta_mu),
    slope_mean   = mean(beta_mu),
    slope_sd= sd(beta_mu),
    CI_width = quantile(beta_mu, 0.975) - quantile(beta_mu, 0.025),
    .groups = "drop")

acap=ggplot(acari_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Acari")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

musp=ggplot(muscidae_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Muscidae")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

chip=ggplot(chironomidae_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Chironomidae")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

dryp=ggplot(dryas_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Dryas")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

saxp=ggplot(saxifraga_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Saxifraga")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

silp=ggplot(silene_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Silene")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

salp=ggplot(ichneumonidae_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Salix")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

papp=ggplot(papaver_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Papaver")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

casp=ggplot(papaver_slope_df,
            aes(x = tsl, y = start_yr)) +
  geom_point(aes(size = 1/CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "Trend certainty\n(1 / CI width)",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Cassiope")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

allreps_tsl=ggarrange(dryp, musp, silp, ichp, papp, chip, ncol=2, nrow=3, legend="none")
annotate_figure(allreps_tsl,
                left = text_grob("95% CI width of slope", rot = 90),
                bottom = text_grob("Time-series length (years)"))

allplants_tsl=ggarrange(dryp, casp, silp, salp, papp, saxp, ncol=2, nrow=3, legend="none")
annotate_figure(allplants_tsl,
                left = text_grob("95% CI width of slope", rot = 90),
                bottom = text_grob("Time-series length (years)"))

allarths1_tsl=ggarrange(musp, ichp, chip, acap, cocp, colp,ncol=2, nrow=3, legend="none")
annotate_figure(allarths1_tsl,
                left = text_grob("95% CI width of slope", rot = 90),
                bottom = text_grob("Time-series length (years)"))

allarths2_tsl=ggarrange(linp, lycp, nymp, phop, scip, ncol=2, nrow=3, legend="none")
annotate_figure(allarths2_tsl,
                left = text_grob("95% CI width of slope", rot = 90),
                bottom = text_grob("Time-series length (years)"))
