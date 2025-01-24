b2=ggplot(bom_dat, aes(x = DOY, y = log(TotalCatch1), col = as.factor(Year))) +
  geom_point() +
  facet_wrap(~Plot.ID) +  # Facet by Plot.ID
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "black") +  # Add spline smoothing
  ggtitle("Bombus") +
  scale_colour_viridis_d() +
  theme_classic()

c2=ggplot(col_dat, aes(x = DOY, y = log(TotalCatch1), col = as.factor(Year))) +
  geom_point() +
  facet_wrap(~Plot.ID) +  # Facet by Plot.ID
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "black") +  # Add spline smoothing
  ggtitle("Collembola") +
  scale_colour_viridis_d() +
  theme_classic()

i2=ggplot(ich_dat, aes(x = DOY, y = log(TotalCatch1), col = as.factor(Year))) +
  geom_point() +
  facet_wrap(~Plot.ID) +  # Facet by Plot.ID
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "black") +  # Add spline smoothing
  ggtitle("Ichneumonidae") +
  scale_colour_viridis_d() +
  theme_classic()

a2=ggplot(aca_dat, aes(x = DOY, y = log(TotalCatch1), col = as.factor(Year))) +
  geom_point() +
  facet_wrap(~Plot.ID) +  # Facet by Plot.ID
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "black") +  # Add spline smoothing
  ggtitle("Acari") +
  scale_colour_viridis_d() +
  theme_classic()

l2=ggplot(lyc_dat, aes(x = DOY, y = log(TotalCatch1), col = as.factor(Year))) +
  geom_point() +
  facet_wrap(~Plot.ID) +  # Facet by Plot.ID
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "black") +  # Add spline smoothing
  ggtitle("Lycosidae") +
  scale_colour_viridis_d() +
  theme_classic()

s2=ggplot(sci_dat, aes(x = DOY, y = log(TotalCatch1), col = as.factor(Year))) +
  geom_point() +
  facet_wrap(~Plot.ID) +  # Facet by Plot.ID
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "black") +  # Add spline smoothing
  ggtitle("Sciaridae") +
  scale_colour_viridis_d() +
  theme_classic()

p2=ggplot(pho_dat, aes(x = DOY, y = log(TotalCatch1), col = as.factor(Year))) +
  geom_point() +
  facet_wrap(~Plot.ID) +  # Facet by Plot.ID
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "black") +  # Add spline smoothing
  ggtitle("Phoridae") +
  scale_colour_viridis_d() +
  theme_classic()

n2=ggplot(nym_dat, aes(x = DOY, y = log(TotalCatch1), col = as.factor(Year))) +
  geom_point() +
  facet_wrap(~Plot.ID) +  # Facet by Plot.ID
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "black") +  # Add spline smoothing
  ggtitle("Nymphidalidae") +
  scale_colour_viridis_d() +
  theme_classic()

m2=ggplot(mus_dat, aes(x = DOY, y = log(TotalCatch1), col = as.factor(Year))) +
  geom_point() +
  facet_wrap(~Plot.ID) +  # Facet by Plot.ID
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "black") +  # Add spline smoothing
  ggtitle("Muscidae") +
  scale_colour_viridis_d() +
  theme_classic()

lin2=ggplot(lin_dat, aes(x = DOY, y = log(TotalCatch1), col = as.factor(Year))) +
  geom_point() +
  facet_wrap(~Plot.ID) +  # Facet by Plot.ID
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "black") +  # Add spline smoothing
  ggtitle("Linyphiidae") +
  scale_colour_viridis_d() +
  theme_classic()

cul2=ggplot(cul_dat, aes(x = DOY, y = log(TotalCatch1), col = as.factor(Year))) +
  geom_point() +
  facet_wrap(~Plot.ID) +  # Facet by Plot.ID
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "black") +  # Add spline smoothing
  ggtitle("Culicidae") +
  scale_colour_viridis_d() +
  theme_classic()

coc2=ggplot(coc_dat, aes(x = DOY, y = log(TotalCatch1), col = as.factor(Year))) +
  geom_point() +
  facet_wrap(~Plot.ID) +  # Facet by Plot.ID
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "black") +  # Add spline smoothing
  ggtitle("Coccidea") +
  scale_colour_viridis_d() +
  theme_classic()

chi2=ggplot(chi_dat, aes(x = DOY, y = log(TotalCatch1), col = as.factor(Year))) +
  geom_point() +
  facet_wrap(~Plot.ID) +  # Facet by Plot.ID
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "black") +  # Add spline smoothing
  ggtitle("Chironomidae") +
  scale_colour_viridis_d() +
  theme_classic()

ggarrange(b2, c2, i2, a2, l2, s2, p2, n2, m2, lin2, cul2, coc2, chi2,
          common.legend = T, legend="right")
