#PLANTS####
dryas_ci_trend <- dryas_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope, 0.05),
    upr_90 = quantile(slope, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

dryci=ggplot(dryas_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Dryas")+
  #  subtitle = "Ribbon shows variability across sliding windows") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

silene_ci_trend <- silene_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope, 0.05),
    upr_90 = quantile(slope, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

silci=ggplot(silene_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Silene")+
  #  subtitle = "Ribbon shows variability across sliding windows") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

papaver_ci_trend <- papaver_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope, 0.05),
    upr_90 = quantile(slope, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

papci=ggplot(papaver_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Papaver")+
  #  subtitle = "Ribbon shows variability across sliding windows") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

#arthropods####
muscidae_ci_trend <- muscidae_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope, 0.05),
    upr_90 = quantile(slope, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

musci=ggplot(muscidae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Muscidae")+
  #  subtitle = "Ribbon shows variability across sliding windows") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

ichneumonidae_ci_trend <- ichneumonidae_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope, 0.05),
    upr_90 = quantile(slope, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

ichci=ggplot(ichneumonidae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Ichneumonidae")+
  #  subtitle = "Ribbon shows variability across sliding windows") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

chironomidae_ci_trend <- chironomidae_slope_draws %>%
  group_by(TSL, start_yr) %>%
  summarise(
    lwr_90 = quantile(slope, 0.05),
    upr_90 = quantile(slope, 0.95),
    .groups = "drop") %>%
  group_by(TSL) %>%
  summarise(
    mean_CI_width = mean(upr_90 - lwr_90),
    CI_lower = quantile(upr_90 - lwr_90, 0.05),
    CI_upper = quantile(upr_90 - lwr_90, 0.95),
    .groups = "drop")

chici=ggplot(chironomidae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Chironomidae")+
  #  subtitle = "Ribbon shows variability across sliding windows") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

allci=ggarrange(dryci,musci, silci, ichci, papci, chici,
                common.legend = T, ncol=2, nrow=3, legend = "right")
annotate_figure(allci,
                 top = text_grob("Trend uncertainty over time",
                                 face = "bold", size = 14),  # bigger title
                 left = text_grob("90% CI width of slope", rot = 90, size = 12),
                 bottom = text_grob("Time-series length (years)", size = 12))

