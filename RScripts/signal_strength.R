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

dryas_abs_slope <- dryas_slope_draws %>%
  group_by(TSL) %>%
  summarise(
    mean_abs_slope = mean(abs(slope)),
    .groups = "drop")%>%arrange(TSL) %>%
  mutate(
    slope_diff = c(NA, diff(mean_abs_slope)),          # slope change per unit TSL
    slope_diff_abs = abs(slope_diff))

dryas_threshold <- 0.05 * max(dryas_abs_slope$slope_diff_abs, na.rm = TRUE)

dryas_pt <- dryas_abs_slope %>%
  filter(slope_diff_abs < dryas_threshold & !is.na(slope_diff_abs)) %>%
  slice(1) %>%pull(TSL)

dryci=ggplot(dryas_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Dryas")+
  geom_vline(xintercept = dryas_pt, linetype = "dashed", color = "darkgreen") +
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

silene_abs_slope <- silene_slope_draws %>%
  group_by(TSL) %>%
  summarise(
    mean_abs_slope = mean(abs(slope)),
    .groups = "drop")%>%arrange(TSL) %>%
  mutate(
    slope_diff = c(NA, diff(mean_abs_slope)),          # slope change per unit TSL
    slope_diff_abs = abs(slope_diff))

silene_threshold <- 0.05 * max(silene_abs_slope$slope_diff_abs, na.rm = TRUE)

silene_pt <- silene_abs_slope %>%
  filter(slope_diff_abs < silene_threshold & !is.na(slope_diff_abs)) %>%
  slice(1) %>%pull(TSL)

silci=ggplot(silene_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Silene")+
  geom_vline(xintercept = silene_pt, linetype = "dashed", color = "darkgreen") +
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
papaver_abs_slope <- papaver_slope_draws %>%
  group_by(TSL) %>%
  summarise(
    mean_abs_slope = mean(abs(slope)),
    .groups = "drop")%>%arrange(TSL) %>%
  mutate(
    slope_diff = c(NA, diff(mean_abs_slope)),          # slope change per unit TSL
    slope_diff_abs = abs(slope_diff))

papaver_threshold <- 0.05 * max(papaver_abs_slope$slope_diff_abs, na.rm = TRUE)

papaver_pt <- papaver_abs_slope %>%
  filter(slope_diff_abs < papaver_threshold & !is.na(slope_diff_abs)) %>%
  slice(1) %>%pull(TSL)

papci=ggplot(papaver_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Papaver")+
  geom_vline(xintercept = papaver_pt , linetype = "dashed", color = "darkgreen") +
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
muscidae_abs_slope <- muscidae_slope_draws %>%
  group_by(TSL) %>%
  summarise(
    mean_abs_slope = mean(abs(slope)),
    .groups = "drop")%>%arrange(TSL) %>%
  mutate(
    slope_diff = c(NA, diff(mean_abs_slope)),          # slope change per unit TSL
    slope_diff_abs = abs(slope_diff))

muscidae_threshold <- 0.05 * max(muscidae_abs_slope$slope_diff_abs, na.rm = TRUE)

muscidae_pt <- muscidae_abs_slope %>%
  filter(slope_diff_abs < muscidae_threshold & !is.na(slope_diff_abs)) %>%
  slice(1) %>%pull(TSL)

musci=ggplot(muscidae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Muscidae")+
  geom_vline(xintercept = muscidae_pt, linetype = "dashed", color = "darkgreen") +
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

ichneumonidae_abs_slope <- ichneumonidae_slope_draws %>%
  group_by(TSL) %>%
  summarise(
    mean_abs_slope = mean(abs(slope)),
    .groups = "drop")%>%arrange(TSL) %>%
  mutate(
    slope_diff = c(NA, diff(mean_abs_slope)),          # slope change per unit TSL
    slope_diff_abs = abs(slope_diff))

ichneumonidae_threshold <- 0.05 * max(ichneumonidae_abs_slope$slope_diff_abs, na.rm = TRUE)

ichneumonidae_pt <- ichneumonidae_abs_slope %>%
  filter(slope_diff_abs < ichneumonidae_threshold & !is.na(slope_diff_abs)) %>%
  slice(1) %>%pull(TSL)

ichci=ggplot(ichneumonidae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Ichneumonidae")+
  geom_vline(xintercept = ichneumonidae_pt, linetype = "dashed", color = "darkgreen") +
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

chironomidae_abs_slope <- chironomidae_slope_draws %>%
  group_by(TSL) %>%
  summarise(
    mean_abs_slope = mean(abs(slope)),
    .groups = "drop")%>%arrange(TSL) %>%
  mutate(
    slope_diff = c(NA, diff(mean_abs_slope)),          # slope change per unit TSL
    slope_diff_abs = abs(slope_diff))

chironomidae_threshold <- 0.05 * max(chironomidae_abs_slope$slope_diff_abs, na.rm = TRUE)

chironomidae_pt <- chironomidae_abs_slope %>%
  filter(slope_diff_abs < chironomidae_threshold & !is.na(slope_diff_abs)) %>%
  slice(1) %>%pull(TSL)

chici=ggplot(chironomidae_ci_trend, aes(x = TSL, y = mean_CI_width)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) +
  geom_line() +geom_point()+
  labs(x = NULL,y = NULL,title = "Chironomidae")+
  geom_vline(xintercept = chironomidae_pt, linetype = "dashed", color = "darkgreen") +
  theme_classic()+theme(plot.title = element_text(face = "bold"))

allci=ggarrange(dryci,musci, silci, ichci, papci, chici,
                common.legend = T, ncol=2, nrow=3, legend = "right")
annotate_figure(allci,
                 left = text_grob("90% CI width of slope", rot = 90, size = 12),
                 bottom = text_grob("Time-series length (years)", size = 12))

