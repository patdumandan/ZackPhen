extract_LT_slope=function(dat, species) {
  dat%>%filter(dat$species== !!species)%>%pull(mean)
}

dryas_LT=extract_LT_slope(plant_beta_mu_summary, "Dryas")

dryas_slope_TSL <- dryas_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=dryas_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

dryas_ci_trend <- dry_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

dryas_pt <- dryas_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

dryci=ggplot(dryas_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Dryas")+
  geom_vline(xintercept = dryas_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

silene_LT=extract_LT_slope(plant_beta_mu_summary, "Silene")

silene_slope_TSL <- silene_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=silene_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

silene_ci_trend <- sil_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

silene_pt <- silene_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)


silci=ggplot(silene_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Silene")+
  geom_vline(xintercept = silene_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

papaver_LT=extract_LT_slope(plant_beta_mu_summary, "Papaver")

papaver_slope_TSL <- papaver_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=papaver_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

papaver_ci_trend <- pap_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

papaver_pt <- papaver_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

salix_LT=extract_LT_slope(plant_beta_mu_summary, "Salix")

salix_slope_TSL <- salix_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=salix_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

salix_ci_trend <- sal_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

salix_pt <- salix_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

salci=ggplot(salix_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Salix")+
  geom_vline(xintercept = salix_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

papci=ggplot(papaver_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Papaver")+
  geom_vline(xintercept = papaver_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

cassiope_LT=extract_LT_slope(plant_beta_mu_summary, "Cassiope")

cassiope_slope_TSL <- cassiope_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=cassiope_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

cassiope_ci_trend <- cas_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

cassiope_pt <- cassiope_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

casci=ggplot(cassiope_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Cassiope")+
  geom_vline(xintercept = cassiope_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

muscidae_LT=extract_LT_slope(arth_beta_mu_summary, "Muscidae")

muscidae_slope_TSL <- muscidae_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=muscidae_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

muscidae_ci_trend <- mus_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

muscidae_pt <- muscidae_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

musci=ggplot(muscidae_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Muscidae")+
  geom_vline(xintercept = muscidae_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

ichneumonidae_LT=extract_LT_slope(arth_beta_mu_summary, "Ichneumonidae")

ichneumonidae_slope_TSL <- ichneumonidae_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=ichneumonidae_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

ichneumonidae_ci_trend <- ich_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

ichneumonidae_pt <- ichneumonidae_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

ichci=ggplot(ichneumonidae_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Ichneumonidae")+
  geom_vline(xintercept = ichneumonidae_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

chironomidae_LT=extract_LT_slope(arth_beta_mu_summary, "Chironomidae")

chironomidae_slope_TSL <- chironomidae_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=chironomidae_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

chironomidae_ci_trend <- chi_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

chironomidae_pt <- chironomidae_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

chici=ggplot(chironomidae_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Chironomidae")+
  geom_vline(xintercept = chironomidae_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

sciaridae_LT=extract_LT_slope(arth_beta_mu_summary, "Sciaridae")

sciaridae_slope_TSL <- sciaridae_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=sciaridae_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

sciaridae_ci_trend <- sci_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

sciaridae_pt <- sciaridae_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

scici=ggplot(sciaridae_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Sciaridae")+
  geom_vline(xintercept = sciaridae_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

nymphalidae_LT=extract_LT_slope(arth_beta_mu_summary, "Nymphalidae")

nymphalidae_slope_TSL <- nymphalidae_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=nymphalidae_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

nymphalidae_ci_trend <- nym_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

nymphalidae_pt <- nymphalidae_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

nymci=ggplot(nymphalidae_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Nymphalidae")+
  geom_vline(xintercept = nymphalidae_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

lycosidae_LT=extract_LT_slope(arth_beta_mu_summary, "Lycosidae")

lycosidae_slope_TSL <- lycosidae_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=lycosidae_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

lycosidae_ci_trend <- lyc_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

lycosidae_pt <- lycosidae_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

lycci=ggplot(lycosidae_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Lycosidae")+
  geom_vline(xintercept = lycosidae_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

collembola_LT=extract_LT_slope(arth_beta_mu_summary, "Collembola")

collembola_slope_TSL <- collembola_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=collembola_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

collembola_ci_trend <- col_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

collembola_pt <- collembola_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

colci=ggplot(collembola_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Collembola")+
  geom_vline(xintercept = collembola_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

acari_LT=extract_LT_slope(arth_beta_mu_summary, "Acari")

acari_slope_TSL <- acari_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=acari_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

acari_ci_trend <- aca_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

acari_pt <- acari_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

acaci=ggplot(acari_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Acari")+
  geom_vline(xintercept = acari_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

coccoidea_LT=extract_LT_slope(arth_beta_mu_summary, "Coccoidea")

coccoidea_slope_TSL <- coccoidea_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=coccoidea_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

coccoidea_ci_trend <- coc_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

coccoidea_pt <- coccoidea_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

cocci=ggplot(coccoidea_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Coccoidea")+
  geom_vline(xintercept = coccoidea_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

linyphiidae_LT=extract_LT_slope(arth_beta_mu_summary, "Linyphiidae")

linyphiidae_slope_TSL <- linyphiidae_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=linyphiidae_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

linyphiidae_ci_trend <- lin_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

linyphiidae_pt <- linyphiidae_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

linci=ggplot(linyphiidae_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Linyphiidae")+
  geom_vline(xintercept = linyphiidae_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

phoridae_LT=extract_LT_slope(arth_beta_mu_summary, "Phoridae")

phoridae_slope_TSL <- phoridae_slope_df%>%
  group_by(tsl) %>%
  mutate(
    LT_slope=phoridae_LT,
    sign_check = if_else(sign(slope_mean) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(tsl)

phoridae_ci_trend <- pho_betamu %>%
  group_by(tsl, start_yr) %>%
  summarise(
    lwr_90 = quantile(beta_mu, 0.05),
    upr_90 = quantile(beta_mu, 0.95),
    .groups = "drop") %>%
  group_by(tsl) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

phoridae_pt <- phoridae_slope_TSL %>%
  group_by(tsl) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_mean) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_mean - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(tsl) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

phoci=ggplot(phoridae_ci_trend, aes(x = tsl, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Phoridae")+
  geom_vline(xintercept = phoridae_pt$tsl, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

repci=ggarrange(dryci, musci, silci, ichci, papci, chici, ncol=2, nrow=3)

annotate_figure(repci,
                left = text_grob("95% CI width of slope", rot = 90),
                bottom = text_grob("Time-series length (years)"))

plantci=ggarrange(dryci, casci, silci, salci, papci, saxci, ncol=2, nrow=3)

annotate_figure(plantci,
                left = text_grob("95% CI width of slope", rot = 90),
                bottom = text_grob("Time-series length (years)"))

arth1ci=ggarrange(musci, ichci, chici, acaci, cocci, colci, ncol=2, nrow=3)
annotate_figure(arth1ci,
                left = text_grob("95% CI width of slope", rot = 90),
                bottom = text_grob("Time-series length (years)"))

arth2ci=ggarrange(linci, lycci, nymci, phoci, scici, ncol=2, nrow=3)
annotate_figure(arth2ci,
                left = text_grob("95% CI width of slope", rot = 90),
                bottom = text_grob("Time-series length (years)"))
