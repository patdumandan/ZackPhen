#PLANTS####
dryas_ci_trend <- dryas_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

#create plant_beta_mu_summary table from plant_stanmod_PD file
dryas_LT=extract_LT_slope(plant_beta_mu_summary, "Dryas")

dryas_slope_TSL <- dryas_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=dryas_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
    #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
    arrange(TSL)

dryas_pt <- dryas_slope_TSL %>%
      group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

dryci=ggplot(dryas_ci_trend, aes(x = TSL, y = mean_CI_width)) +
      geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
      geom_line() +geom_point()+
      labs(x = NULL,y = NULL,title = "Dryas")+
      geom_vline(xintercept = dryas_pt$TSL, linetype = "dashed", color = "darkgreen") +
      theme_classic()+theme(plot.title = element_text(face = "bold"))

silene_ci_trend <- silene_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

#create plant_beta_mu_summary table from plant_stanmod_PD file
silene_LT=extract_LT_slope(plant_beta_mu_summary, "Silene")

silene_slope_TSL <- silene_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=silene_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

silene_pt <- silene_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
#  filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
 filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

silci=ggplot(silene_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Silene")+
  geom_vline(xintercept = silene_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

cassiope_ci_trend <- cassiope_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

#create plant_beta_mu_summary table from plant_stanmod_PD file
cassiope_LT=extract_LT_slope(plant_beta_mu_summary, "Cassiope")

cassiope_slope_TSL <- cassiope_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=cassiope_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

cassiope_pt <- cassiope_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
#  filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

casci=ggplot(cassiope_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Cassiope")+
  geom_vline(xintercept = cassiope_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

salix_ci_trend <- salix_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

#create plant_beta_mu_summary table from plant_stanmod_PD file
salix_LT=extract_LT_slope(plant_beta_mu_summary, "Salix")

salix_slope_TSL <- salix_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=salix_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

salix_pt <- salix_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

salci=ggplot(salix_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Salix")+
  geom_vline(xintercept = salix_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

papaver_ci_trend <- papaver_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

#create plant_beta_mu_summary table from plant_stanmod_PD file
papaver_LT=extract_LT_slope(plant_beta_mu_summary, "Papaver")

papaver_slope_TSL <- papaver_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=papaver_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

papaver_pt <- papaver_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
#  filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
 filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

papci=ggplot(papaver_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Papaver")+
  geom_vline(xintercept = papaver_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

saxifraga_ci_trend <- saxifraga_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

#create plant_beta_mu_summary table from plant_stanmod_PD file
saxifraga_LT=extract_LT_slope(plant_beta_mu_summary, "Saxifraga")

saxifraga_slope_TSL <- saxifraga_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=saxifraga_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

saxifraga_pt <- saxifraga_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
#  filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

saxci=ggplot(saxifraga_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Saxifraga")+
  geom_vline(xintercept = saxifraga_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

ggarrange(dryci, silci, papci, salci, saxci, casci, ncol=2, nrow=3)

#arthropods####
muscidae_ci_trend <- muscidae_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

#create plant_beta_mu_summary table from plant_stanmod_PD file
muscidae_LT=extract_LT_slope(arth_beta_mu_summary, "Muscidae")

muscidae_slope_TSL <- muscidae_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=muscidae_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

muscidae_pt <- muscidae_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
#  filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
 filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

musci=ggplot(muscidae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Muscidae")+
  geom_vline(xintercept = muscidae_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

ichneumonidae_ci_trend <- ichneumonidae_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

#create plant_beta_mu_summary table from plant_stanmod_PD file
ichneumonidae_LT=extract_LT_slope(arth_beta_mu_summary, "Ichneumonidae")

ichneumonidae_slope_TSL <- ichneumonidae_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=ichneumonidae_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

ichneumonidae_pt <- ichneumonidae_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
#  filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
 filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

ichci=ggplot(ichneumonidae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Ichneumonidae")+
  geom_vline(xintercept = ichneumonidae_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

chironomidae_ci_trend <- chironomidae_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

chironomidae_LT=extract_LT_slope(arth_beta_mu_summary, "Chironomidae")

chironomidae_slope_TSL <- chironomidae_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=chironomidae_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

chironomidae_pt <- chironomidae_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
#  filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
 filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

chici=ggplot(chironomidae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Chironomidae")+
  geom_vline(xintercept = chironomidae_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

lycosidae_ci_trend <- lycosidae_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

lycosidae_LT=extract_LT_slope(arth_beta_mu_summary, "Lycosidae")

lycosidae_slope_TSL <- lycosidae_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=lycosidae_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

lycosidae_pt <- lycosidae_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  # pct_reduction  = rmse_reduction / lag(median_rmse)) %>%
#  filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
 filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

lycci=ggplot(lycosidae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Lycosidae")+
  geom_vline(xintercept = lycosidae_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

phoridae_ci_trend <- phoridae_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

#create plant_beta_mu_summary table from plant_stanmod_PD file
phoridae_LT=extract_LT_slope(arth_beta_mu_summary, "Phoridae")

phoridae_slope_TSL <- phoridae_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=phoridae_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

phoridae_pt <- phoridae_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
#    filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
 filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

phoci=ggplot(phoridae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Phoridae")+
  geom_vline(xintercept = phoridae_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

sciaridae_ci_trend <- sciaridae_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

sciaridae_LT=extract_LT_slope(arth_beta_mu_summary, "Sciaridae")

sciaridae_slope_TSL <- sciaridae_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=sciaridae_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

sciaridae_pt <- sciaridae_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
#  filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
 filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

scici=ggplot(sciaridae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Sciaridae")+
  geom_vline(xintercept = sciaridae_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

acari_ci_trend <- acari_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

acari_LT=extract_LT_slope(arth_beta_mu_summary, "Acari")

acari_slope_TSL <- acari_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=acari_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

acari_pt <- acari_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
#  filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
 filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

acaci=ggplot(acari_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Acari")+
  geom_vline(xintercept = acari_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

collembola_ci_trend <- collembola_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

#create plant_beta_mu_summary table from plant_stanmod_PD file
collembola_LT=extract_LT_slope(arth_beta_mu_summary, "Collembola")

collembola_slope_TSL <- collembola_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=collembola_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

collembola_pt <- collembola_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

colci=ggplot(collembola_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Collembola")+
  geom_vline(xintercept = collembola_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

coccoidea_ci_trend <- coccoidea_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

#create plant_beta_mu_summary table from plant_stanmod_PD file
coccoidea_LT=extract_LT_slope(arth_beta_mu_summary, "Coccoidea")

coccoidea_slope_TSL <- coccoidea_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=coccoidea_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

coccoidea_pt <- coccoidea_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

cocci=ggplot(coccoidea_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Coccoidea")+
  geom_vline(xintercept = coccoidea_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

nymphalidae_ci_trend <- nymphalidae_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

#create plant_beta_mu_summary table from plant_stanmod_PD file
nymphalidae_LT=extract_LT_slope(arth_beta_mu_summary, "Nymphalidae")

nymphalidae_slope_TSL <- nymphalidae_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=nymphalidae_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

nymphalidae_pt <- nymphalidae_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

nymci=ggplot(nymphalidae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Nymphalidae")+
  geom_vline(xintercept = nymphalidae_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

linyphiidae_ci_trend <- linyphiidae_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope_draw, 0.05),
    upr_90 = quantile(slope_draw, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

#create plant_beta_mu_summary table from plant_stanmod_PD file
linyphiidae_LT=extract_LT_slope(arth_beta_mu_summary, "Linyphiidae")

linyphiidae_slope_TSL <- linyphiidae_slope_draws%>%
  group_by(TSL) %>%
  mutate(
    LT_slope=linyphiidae_LT,
    sign_check = if_else(sign(slope_draw) == sign(LT_slope),"correct","wrong"))%>%
  #  rmse_component = (slope_draw - LT_slope)^2 / n()) %>%
  arrange(TSL)

linyphiidae_pt <- linyphiidae_slope_TSL %>%
  group_by(TSL) %>%
  summarise(
    p_agree_LT =  mean(sign(slope_draw) == sign(LT_slope)),
    #   sign_certainty = pmax(p_agree_LT, 1 - p_agree_LT),
    rmse = sqrt(mean((slope_draw - LT_slope)^2, na.rm = TRUE)),
    .groups = "drop") %>%
  arrange(TSL) %>%
  mutate(
    rmse_reduction = lag(rmse) - rmse,
    pct_reduction  = rmse_reduction / lag(rmse))%>%
  #filter(!is.na(pct_reduction), pct_reduction < 0.05) %>%
  filter(p_agree_LT >= 0.95)%>% #select min.year when the mean prob of TSL slope=TS_slope >0.95
  slice(1)

linci=ggplot(linyphiidae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Linyphiidae")+
  geom_vline(xintercept = linyphiidae_pt$TSL, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

