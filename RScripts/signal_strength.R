dryas_ci_trend <- dryas_peak_slope_ci %>%
  group_by(TSL) %>%
  summarise(
    mean_CI = mean(CI_width),
    median_CI = median(CI_width),
    sd_CI = sd(CI_width)
  )

dryci=ggplot(dryas_ci_trend, aes(x = TSL, y = mean_CI)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = mean_CI - sd_CI, ymax = mean_CI + sd_CI), alpha = 0.2) +
  theme_classic() +
  labs(x = "Time-series length (years)",
       y = "Mean 95% CI width",
       title= "Dryas",
       subtitle = "Phenological signal with increasing timeseries length")+
  theme(plot.title = element_text(face = "bold"))


dryas_ci_fit <- lm(mean_CI ~ I(1/sqrt(TSL)), data = dryas_ci_trend)
summary(dryas_ci_fit)

muscidae_ci_trend <- muscidae_peak_slope_ci %>%
  group_by(TSL) %>%
  summarise(
    mean_CI = mean(CI_width),
    median_CI = median(CI_width),
    sd_CI = sd(CI_width)
  )

musci=ggplot(muscidae_ci_trend, aes(x = TSL, y = mean_CI)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = mean_CI - sd_CI, ymax = mean_CI + sd_CI), alpha = 0.2) +
  theme_classic() +
  labs(x = "Time-series length (years)",
       y = "Mean 95% CI width",
       title= "Muscidae",
       subtitle = "Phenological signal with increasing timeseries length")+
  theme(plot.title = element_text(face = "bold"))

muscidae_ci_fit <- lm(mean_CI ~ I(1/sqrt(TSL)), data = muscidae_ci_trend)
summary(muscidae_ci_fit)

ggarrange(dryci, musci, common.legend = T)
