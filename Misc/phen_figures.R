dryfitted_df=dryfitted_df%>%mutate(species="Dryas")
casfitted_df=casfitted_df%>%mutate(species="Cassiope")
salfitted_df=salfitted_df%>%mutate(species="Salix")
silfitted_df=silfitted_df%>%mutate(species="Silene")
saxfitted_df=saxfitted_df%>%mutate(species="Saxifraga")
papfitted_df=papfitted_df%>%mutate(species="Papaver")

all_fits=rbind(dryfitted_df,casfitted_df,papfitted_df,salfitted_df,silfitted_df,
               saxfitted_df)

# Define all years present
years_all <- sort(unique(all_fits$year))
all_fits$year <- factor(as.numeric(as.character(all_fits$year)), levels = sort(unique(as.numeric(as.character(all_fits$year)))))

ggplot(all_fits, aes(x = DOY, y = prob, col = year, group = interaction(species, year))) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  theme_classic() +
  labs(x = "DOY", y = "P(flower)", title = "Predicted flowering", color = "Year") +
  scale_color_viridis_d(drop = FALSE) +  # prevent dropping unused levels
  facet_wrap(~species) +
  xlim(150, 270)

unique(arth_df$HoyeTaxon)

chifitted_df=chifitted_df%>%mutate(species="Chironomidae")
colfitted_df=colfitted_df%>%mutate(species="Collembola")
ichfitted_df=ichfitted_df%>%mutate(species="Ichneumonidae")
linfitted_df=linfitted_df%>%mutate(species="Linyphiidae")
musfitted_df=musfitted_df%>%mutate(species="Muscidae")
scifitted_df=scifitted_df%>%mutate(species="Sciaridae")
cocfitted_df=cocfitted_df%>%mutate(species="Coccoidea")
nymfitted_df=nymfitted_df%>%mutate(species="Nymphalidae")
phofitted_df=phofitted_df%>%mutate(species="Phoridae")
acafitted_df=acafitted_df%>%mutate(species="Acari")
lycfitted_df=lycfitted_df%>%mutate(species="Lycosidae")

arth_fits=rbind(chifitted_df,colfitted_df,ichfitted_df,linfitted_df,musfitted_df,
               scifitted_df,cocfitted_df, nymfitted_df, phofitted_df,
               acafitted_df, lycfitted_df)

# Define all years present
years_all <- sort(unique(arth_fits$year))
arth_fits$year <- factor(as.numeric(as.character(arth_fits$year)),
                         levels = sort(unique(as.numeric(as.character(arth_fits$year)))))

ggplot(arth_fits, aes(x = DOY, y = prob, col = year, group = interaction(species, year))) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  theme_classic() +
  labs(x = "DOY", y = "predicted proportions", title = "Insect Emergence", color = "Year") +
  scale_color_viridis_d(drop = FALSE) +  # prevent dropping unused levels
  facet_wrap(~species) +
  xlim(150, 270)
