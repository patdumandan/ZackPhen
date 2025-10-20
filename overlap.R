#compile fitted curves

musfitted_df=musfitted_df%>%mutate(taxon="Muscidae")
colfitted_df=colfitted_df%>%mutate(taxon="Collembola")
chifitted_df=chifitted_df%>%mutate(taxon="Chironomidae")
ichfitted_df=ichfitted_df%>%mutate(taxon="Ichneumonidae")
linfitted_df=linfitted_df%>%mutate(taxon="Linyphiidae")
scifitted_df=scifitted_df%>%mutate(taxon="Sciaridae")
cocfitted_df=cocfitted_df%>%mutate(taxon="Coccoidea")
nymfitted_df=nymfitted_df%>%mutate(taxon="Nymphalidae")
phofitted_df=phofitted_df%>%mutate(taxon="Phoridae")
acafitted_df=acafitted_df%>%mutate(taxon="Acari")
lycfitted_df=lycfitted_df%>%mutate(taxon="Lycosidae")


art_fc=do.call(rbind, list(musfitted_df, colfitted_df, chifitted_df, ichfitted_df,
                           linfitted_df, scifitted_df, cocfitted_df, nymfitted_df,
                           phofitted_df, acafitted_df, lycfitted_df))

casfitted_df=casfitted_df%>%mutate(taxon="Cassiope")
dryfitted_df=dryfitted_df%>%mutate(taxon="Dryas")
papfitted_df=papfitted_df%>%mutate(taxon="Papaver")
salfitted_df=salfitted_df%>%mutate(taxon="Salix")
saxfitted_df=saxfitted_df%>%mutate(taxon="Saxifraga")
silfitted_df=silfitted_df%>%mutate(taxon="Silene")

plant_fc=do.call(rbind, list(casfitted_df, dryfitted_df, papfitted_df,
                             saxfitted_df, salfitted_df, silfitted_df))

all_fc=do.call(rbind, list(art_fc, plant_fc))

dry_fc=do.call(rbind, list(dryfitted_df, musfitted_df, colfitted_df, chifitted_df, ichfitted_df,
                           linfitted_df, scifitted_df, cocfitted_df, nymfitted_df,
                           phofitted_df, acafitted_df, lycfitted_df))

dry_fc <- dry_fc %>%
  group_by(taxon, year) %>%
  mutate(prob_norm = prob / sum(prob, na.rm = TRUE)) %>%
  ungroup()

# Separate dry taxon
dry_curves <- dry_fc %>% filter(taxon == "Dryas") %>% select(DOY, year, prob_norm) %>%
  rename(prob_dry = prob_norm)

# Join with all other taxon by DOY and year
dryoverlap_df <- dry_fc %>%
  filter(taxon != "Dryas") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(dry_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_dry), na.rm = TRUE)) %>%
  ungroup()

dryoverlap_df$year <- factor(dryoverlap_df$year, levels = sort(as.numeric(levels(dryoverlap_df$year))))

dryo=ggplot(dryoverlap_df, aes(x = year, y = overlap, color = taxon, group = taxon)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_classic() +
  labs(title = "Trend of Phenological Overlap with Dryas Over Time",
       x = "Year",
       y = "Overlap (Schoener's D)",
       color = "Taxon") +
  scale_color_viridis_d() +  # colorblind-friendly palette
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  facet_wrap(~taxon)


#fitted cirves overlay
library(dplyr)
library(ggplot2)

# Assuming dry_fc contains prob_norm for all species except Dryas
# dry_curves contains prob_dry for Dryas
# We'll join the two datasets on DOY and year, for each species

# Prepare Dryas curves with appropriate columns
dry_curves_prepped <- dry_curves %>%
  select(DOY, year, prob_dry) %>%
  mutate(taxon = "Dryas")

# Prepare other species data
other_species <- dry_fc %>%
  filter(taxon != "Dryas") %>%
  select(DOY, year, taxon, prob_norm)

# Join Dryas curve to each species/year combination on DOY and year
# To plot Dryas and other species curves together, we'll gather them into long format

# First, join Dryas curves to other species curves by DOY and year
comparison_curves <- other_species %>%
  left_join(dry_curves_prepped, by = c("DOY", "year")) %>%
  rename(prob_other = prob_norm, prob_dryas = prob_dry) %>%
  select(DOY, year, taxon, prob_other, prob_dryas)

# Pivot longer to get a 'species' column with values 'Dryas' and the other species
plot_df <- comparison_curves %>%
  pivot_longer(cols = c(prob_other, prob_dryas),
               names_to = "species_curve",
               values_to = "prob") %>%
  mutate(species_curve = recode(species_curve,
                                prob_other = taxon.x,
                                prob_dryas = "Dryas"))

# Now plot:
ggplot(plot_df, aes(x = DOY, y = prob, color = species_curve, group = interaction(year, species_curve))) +
  geom_line(size = 1, alpha = 0.8) +
  facet_wrap(~taxon.x, scales = "free_y") +
  theme_classic() +
  scale_color_manual(values = c("Dryas" = "black",
                                # all other species set to gray:
                                setNames(rep("gray50", length(unique(plot_df$species_curve)) - 1),
                                         unique(plot_df$species_curve)[unique(plot_df$species_curve) != "Dryas"]))) +
  labs(title = "Dryas vs. Arthropods",
       x = "Day of Year (DOY)",
       y = "Normalized Probability",
       color = "Species") +
  theme(legend.position = "none")

#cassiop
cas_fc=do.call(rbind, list(casfitted_df, musfitted_df, colfitted_df, chifitted_df, ichfitted_df,
                           linfitted_df, scifitted_df, cocfitted_df, nymfitted_df,
                           phofitted_df, acafitted_df, lycfitted_df))

cas_fc <- cas_fc %>%
  group_by(taxon, year) %>%
  mutate(prob_norm = prob / sum(prob, na.rm = TRUE)) %>%
  ungroup()

# Separate cas taxon
cas_curves <- cas_fc %>% filter(taxon == "Cassiope") %>% select(DOY, year, prob_norm) %>%
  rename(prob_cas = prob_norm)

# Join with all other taxon by DOY and year
casoverlap_df <- cas_fc %>%
  filter(taxon != "Cassiope") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(cas_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_cas), na.rm = TRUE)) %>%
  ungroup()
casoverlap_df$year <- factor(casoverlap_df$year, levels = sort(as.numeric(levels(casoverlap_df$year))))

caso=ggplot(casoverlap_df, aes(x = year, y = overlap, color = taxon, group = taxon)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_classic() +
  labs(title = "Trend of Phenological Overlap with Cassiope Over Time",
       x = "Year",
       y = "Overlap (Schoener's D)",
       color = "Taxon") +
  scale_color_viridis_d() +  # colorblind-friendly palette
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  facet_wrap(~taxon)
# Prepare Cassiope curves with appropriate columns
cas_curves_prepped <- cas_curves %>%
  select(DOY, year, prob_cas) %>%
  mutate(taxon = "Cassiope")

# Prepare other species data
other_species <- cas_fc %>%
  filter(taxon != "Cassiope") %>%
  select(DOY, year, taxon, prob_norm)

# Join Cassiope curve to each species/year combination on DOY and year
# To plot Cassiope and other species curves together, we'll gather them into long format

# First, join Cassiope curves to other species curves by DOY and year
comparison_curves <- other_species %>%
  left_join(cas_curves_prepped, by = c("DOY", "year")) %>%
  rename(prob_other = prob_norm, prob_cas = prob_cas)

# Pivot longer to get a 'species' column with values 'Cassiope' and the other species
plot_df <- comparison_curves %>%
  pivot_longer(cols = c(prob_other, prob_cas),
               names_to = "species_curve",
               values_to = "prob") %>%
  mutate(species_curve = recode(species_curve,
                                prob_other = taxon.x,
                                prob_cas = "Cassiope"))

# Now plot:
ggplot(plot_df, aes(x = DOY, y = prob, color = species_curve, group = interaction(year, species_curve))) +
  geom_line(size = 1, alpha = 0.8) +
  facet_wrap(~taxon.x, scales = "free_y") +
  theme_classic() +
  scale_color_manual(values = c("Cassiope" = "black",
                                # all other species set to gray:
                                setNames(rep("gray50", length(unique(plot_df$species_curve)) - 1),
                                         unique(plot_df$species_curve)[unique(plot_df$species_curve) != "Cassiope"]))) +
  labs(title = "Cassiope vs. Arthropods",
       x = "Day of Year (DOY)",
       y = "Normalized Probability",
       color = "Species") +
  theme(legend.position = "none")

#papaver
pap_fc=do.call(rbind, list(papfitted_df, musfitted_df, colfitted_df, chifitted_df, ichfitted_df,
                           linfitted_df, scifitted_df, cocfitted_df, nymfitted_df,
                           phofitted_df, acafitted_df, lycfitted_df))

pap_fc <- pap_fc %>%
  group_by(taxon, year) %>%
  mutate(prob_norm = prob / sum(prob, na.rm = TRUE)) %>%
  ungroup()

# Separate pap taxon
pap_curves <- pap_fc %>% filter(taxon == "Papaver") %>% select(DOY, year, prob_norm) %>%
  rename(prob_pap = prob_norm)

# Join with all other taxon by DOY and year
papoverlap_df <- pap_fc %>%
  filter(taxon != "Papaver") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(pap_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_pap), na.rm = TRUE)) %>%
  ungroup()

papoverlap_df$year <- factor(papoverlap_df$year, levels = sort(as.numeric(levels(papoverlap_df$year))))

papo=ggplot(papoverlap_df, aes(x = year, y = overlap, color = taxon, group = taxon)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_classic() +
  labs(title = "Trend of Phenological Overlap with Papaver Over Time",
       x = "Year",
       y = "Overlap (Schoener's D)",
       color = "Taxon") +
  scale_color_viridis_d() +  # colorblind-friendly palette
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  facet_wrap(~taxon)
# Prepare Papaver curves with appropriate columns
pap_curves_prepped <- pap_curves %>%
  select(DOY, year, prob_pap) %>%
  mutate(taxon = "Papaver")

# Prepare other species data
other_species <- pap_fc %>%
  filter(taxon != "Papaver") %>%
  select(DOY, year, taxon, prob_norm)

# Join Papaver curve to each species/year combination on DOY and year
# To plot Papaver and other species curves together, we'll gather them into long format

# First, join Papaver curves to other species curves by DOY and year
comparison_curves <- other_species %>%
  left_join(pap_curves_prepped, by = c("DOY", "year")) %>%
  rename(prob_other = prob_norm, prob_pap = prob_pap)

# Pivot longer to get a 'species' column with values 'Papaver' and the other species
plot_df <- comparison_curves %>%
  pivot_longer(cols = c(prob_other, prob_pap),
               names_to = "species_curve",
               values_to = "prob") %>%
  mutate(species_curve = recode(species_curve,
                                prob_other = taxon.x,
                                prob_pap = "Papaver"))

# Now plot:
ggplot(plot_df, aes(x = DOY, y = prob, color = species_curve, group = interaction(year, species_curve))) +
  geom_line(size = 1, alpha = 0.8) +
  facet_wrap(~taxon.x, scales = "free_y") +
  theme_classic() +
  scale_color_manual(values = c("Papaver" = "black",
                                # all other species set to gray:
                                setNames(rep("gray50", length(unique(plot_df$species_curve)) - 1),
                                         unique(plot_df$species_curve)[unique(plot_df$species_curve) != "Papaver"]))) +
  labs(title = "Papaver vs. Arthropods",
       x = "Day of Year (DOY)",
       y = "Normalized Probability",
       color = "Species") +
  theme(legend.position = "none")

#salix
sal_fc=do.call(rbind, list(salfitted_df, musfitted_df, colfitted_df, chifitted_df, ichfitted_df,
                           linfitted_df, scifitted_df, cocfitted_df, nymfitted_df,
                           phofitted_df, acafitted_df, lycfitted_df))

sal_fc <- sal_fc %>%
  group_by(taxon, year) %>%
  mutate(prob_norm = prob / sum(prob, na.rm = TRUE)) %>%
  ungroup()

# Separate sal taxon
sal_curves <- sal_fc %>% filter(taxon == "Salix") %>% select(DOY, year, prob_norm) %>%
  rename(prob_sal = prob_norm)

# Join with all other taxon by DOY and year
saloverlap_df <- sal_fc %>%
  filter(taxon != "Salix") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(sal_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_sal), na.rm = TRUE)) %>%
  ungroup()

saloverlap_df$year <- factor(saloverlap_df$year, levels = sort(as.numeric(levels(saloverlap_df$year))))

salo=ggplot(saloverlap_df, aes(x = year, y = overlap, color = taxon, group = taxon)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_classic() +
  labs(title = "Trend of Phenological Overlap with Salix Over Time",
       x = "Year",
       y = "Overlap (Schoener's D)",
       color = "Taxon") +
  scale_color_viridis_d() +  # colorblind-friendly palette
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  facet_wrap(~taxon)


# Prepare Salix curves with appropriate columns
sal_curves_prepped <- sal_curves %>%
  select(DOY, year, prob_sal) %>%
  mutate(taxon = "Salix")

# Prepare other species data
other_species <- sal_fc %>%
  filter(taxon != "Salix") %>%
  select(DOY, year, taxon, prob_norm)

# Join Salix curve to each species/year combination on DOY and year
# To plot Salix and other species curves together, we'll gather them into long format

# First, join Salix curves to other species curves by DOY and year
comparison_curves <- other_species %>%
  left_join(sal_curves_prepped, by = c("DOY", "year")) %>%
  rename(prob_other = prob_norm, prob_sal = prob_sal)

# Pivot longer to get a 'species' column with values 'Salix' and the other species
plot_df <- comparison_curves %>%
  pivot_longer(cols = c(prob_other, prob_sal),
               names_to = "species_curve",
               values_to = "prob") %>%
  mutate(species_curve = recode(species_curve,
                                prob_other = taxon.x,
                                prob_sal = "Salix"))

# Now plot:
ggplot(plot_df, aes(x = DOY, y = prob, color = species_curve, group = interaction(year, species_curve))) +
  geom_line(size = 1, alpha = 0.8) +
  facet_wrap(~taxon.x, scales = "free_y") +
  theme_classic() +
  scale_color_manual(values = c("Salix" = "black",
                                # all other species set to gray:
                                setNames(rep("gray50", length(unique(plot_df$species_curve)) - 1),
                                         unique(plot_df$species_curve)[unique(plot_df$species_curve) != "Salix"]))) +
  labs(title = "Salix vs. Arthropods",
       x = "Day of Year (DOY)",
       y = "Normalized Probability",
       color = "Species") +
  theme(legend.position = "none")


#silene
sil_fc=do.call(rbind, list(silfitted_df, musfitted_df, colfitted_df, chifitted_df, ichfitted_df,
                           linfitted_df, scifitted_df, cocfitted_df, nymfitted_df,
                           phofitted_df, acafitted_df, lycfitted_df))

sil_fc <- sil_fc %>%
  group_by(taxon, year) %>%
  mutate(prob_norm = prob / sum(prob, na.rm = TRUE)) %>%
  ungroup()

# Separate sil taxon
sil_curves <- sil_fc %>% filter(taxon == "Silene") %>% select(DOY, year, prob_norm) %>%
  rename(prob_sil = prob_norm)

# Join with all other taxon by DOY and year
siloverlap_df <- sil_fc %>%
  filter(taxon != "Silene") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(sil_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_sil), na.rm = TRUE)) %>%
  ungroup()
siloverlap_df$year <- factor(siloverlap_df$year, levels = sort(as.numeric(levels(siloverlap_df$year))))

silo=ggplot(siloverlap_df, aes(x = year, y = overlap, color = taxon, group = taxon)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_classic() +
  labs(title = "Trend of Phenological Overlap with Silene Over Time",
       x = "Year",
       y = "Overlap (Schoener's D)",
       color = "Taxon") +
  scale_color_viridis_d() +  # colorblind-friendly palette
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  facet_wrap(~taxon)


# Prepare Silene curves with appropriate columns
sil_curves_prepped <- sil_curves %>%
  select(DOY, year, prob_sil) %>%
  mutate(taxon = "Silene")

# Prepare other species data
other_species <- sil_fc %>%
  filter(taxon != "Silene") %>%
  select(DOY, year, taxon, prob_norm)

# Join Silene curve to each species/year combination on DOY and year
# To plot Silene and other species curves together, we'll gather them into long format

# First, join Silene curves to other species curves by DOY and year
comparison_curves <- other_species %>%
  left_join(sil_curves_prepped, by = c("DOY", "year")) %>%
  rename(prob_other = prob_norm, prob_sil = prob_sil)

# Pivot longer to get a 'species' column with values 'Silene' and the other species
plot_df <- comparison_curves %>%
  pivot_longer(cols = c(prob_other, prob_sil),
               names_to = "species_curve",
               values_to = "prob") %>%
  mutate(species_curve = recode(species_curve,
                                prob_other = taxon.x,
                                prob_sil = "Silene"))

# Now plot:
ggplot(plot_df, aes(x = DOY, y = prob, color = species_curve, group = interaction(year, species_curve))) +
  geom_line(size = 1, alpha = 0.8) +
  facet_wrap(~taxon.x, scales = "free_y") +
  theme_classic() +
  scale_color_manual(values = c("Silene" = "black",
                                # all other species set to gray:
                                setNames(rep("gray50", length(unique(plot_df$species_curve)) - 1),
                                         unique(plot_df$species_curve)[unique(plot_df$species_curve) != "Silene"]))) +
  labs(title = "Silene vs. Arthropods",
       x = "Day of Year (DOY)",
       y = "Normalized Probability",
       color = "Species") +
  theme(legend.position = "none")

#saxifraga
sax_fc=do.call(rbind, list(saxfitted_df, musfitted_df, colfitted_df, chifitted_df, ichfitted_df,
                           linfitted_df, scifitted_df, cocfitted_df, nymfitted_df,
                           phofitted_df, acafitted_df, lycfitted_df))

sax_fc <- sax_fc %>%
  group_by(taxon, year) %>%
  mutate(prob_norm = prob / sum(prob, na.rm = TRUE)) %>%
  ungroup()

# Separate sax taxon
sax_curves <- sax_fc %>% filter(taxon == "Saxifraga") %>% select(DOY, year, prob_norm) %>%
  rename(prob_sax = prob_norm)

# Join with all other taxon by DOY and year
saxoverlap_df <- sax_fc %>%
  filter(taxon != "Saxifraga") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(sax_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_sax), na.rm = TRUE)) %>%
  ungroup()

saxo=ggplot(saxoverlap_df, aes(x = year, y = overlap, color = taxon, group = taxon)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_classic() +
  labs(title = "Trend of Phenological Overlap with Saxifraga Over Time",
       x = "Year",
       y = "Overlap (Schoener's D)",
       color = "Species") +
  scale_color_viridis_d() +  # colorblind-friendly palette
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))

# Prepare Saxifraga curves with appropriate columns
sax_curves_prepped <- sax_curves %>%
  select(DOY, year, prob_sax) %>%
  mutate(taxon = "Saxifraga")

# Prepare other species data
other_species <- sax_fc %>%
  filter(taxon != "Saxifraga") %>%
  select(DOY, year, taxon, prob_norm)

# Join Saxifraga curve to each species/year combination on DOY and year
# To plot Saxifraga and other species curves together, we'll gather them into long format

# First, join Saxifraga curves to other species curves by DOY and year
comparison_curves <- other_species %>%
  left_join(sax_curves_prepped, by = c("DOY", "year")) %>%
  rename(prob_other = prob_norm, prob_sax = prob_sax)

# Pivot longer to get a 'species' column with values 'Saxifraga' and the other species
plot_df <- comparison_curves %>%
  pivot_longer(cols = c(prob_other, prob_sax),
               names_to = "species_curve",
               values_to = "prob") %>%
  mutate(species_curve = recode(species_curve,
                                prob_other = taxon.x,
                                prob_sax = "Saxifraga"))

# Now plot:
ggplot(plot_df, aes(x = DOY, y = prob, color = species_curve, group = interaction(year, species_curve))) +
  geom_line(size = 1, alpha = 0.8) +
  facet_wrap(~taxon.x, scales = "free_y") +
  theme_classic() +
  scale_color_manual(values = c("Saxifraga" = "black",
                                # all other species set to gray:
                                setNames(rep("gray50", length(unique(plot_df$species_curve)) - 1),
                                         unique(plot_df$species_curve)[unique(plot_df$species_curve) != "Saxifraga"]))) +
  labs(title = "Saxifraga vs. Arthropods",
       x = "Day of Year (DOY)",
       y = "Normalized Probability",
       color = "Species") +
  theme(legend.position = "none")

#combine subsets
casoverlap_df <- cas_fc %>%
  filter(taxon != "Cassiope") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(cas_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_cas), na.rm = TRUE)) %>%
  ungroup()

dryoverlap_df <- dry_fc %>%
  filter(taxon != "Dryas") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(dry_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_dry), na.rm = TRUE)) %>%
  ungroup()

papoverlap_df <- pap_fc %>%
  filter(taxon != "Papaver") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(pap_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_pap), na.rm = TRUE)) %>%
  ungroup()

saloverlap_df <- sal_fc %>%
  filter(taxon != "Salix") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(sal_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_sal), na.rm = TRUE)) %>%
  ungroup()

siloverlap_df <- sil_fc %>%
  filter(taxon != "Silene") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(sil_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_sil), na.rm = TRUE)) %>%
  ungroup()

saxoverlap_df <- sax_fc %>%
  filter(taxon != "Saxifraga") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(sax_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_sax), na.rm = TRUE)) %>%
  ungroup()

all_overlap_df=do.call(rbind, list(casoverlap_df, dryoverlap_df,
                                   papoverlap_df, saloverlap_df,
                                   siloverlap_df,saxoverlap_df))

all_overlap_df=as.data.frame(all_overlap_df)



