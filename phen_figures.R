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
