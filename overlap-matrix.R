#Silene####
siloverlap_df <- sil_fc %>%
  filter(taxon != "Silene") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(sil_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_sil), na.rm = TRUE)) %>%
  ungroup()

siloverlap_df_summ <- sil_fc %>%
  filter(taxon != "Silene") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(sil_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_sil), na.rm = TRUE)) %>%
  ungroup()%>%
  group_by(taxon)%>%
  summarise(mean_overlap=mean(overlap))

siloverlap_df$year <- factor(siloverlap_df$year, levels = sort(as.numeric(levels(siloverlap_df$year))))

silo=ggplot(siloverlap_df, aes(x = year, y = overlap, color = taxon, group = taxon)) +
  geom_smooth(method="lm") +
  geom_point(size = 2)+
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

siloverlap_trends <- siloverlap_df %>%mutate(year = as.numeric(year))%>%
  group_by(taxon) %>%
  filter(n_distinct(year) > 1) %>%
  do({
    model <- lm(overlap ~ year, data = .)
    p_val <- summary(model)$coefficients["year", "Pr(>|t|)"]
    data.frame(
      slope = coef(model)[["year"]],
      p_value = p_val,
      r_squared = summary(model)$r.squared,
      significant = ifelse(p_val < 0.05, TRUE, FALSE)
    ) }) %>%ungroup()

silov=left_join(siloverlap_trends, siloverlap_df_summ)%>%
  mutate(resource="Silene")%>%
  rename(consumer=taxon)

#Papaver####
papoverlap_df<- pap_fc %>%
  filter(taxon != "Papaver") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(pap_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_pap), na.rm = TRUE)) %>%
  ungroup()%>%
  mutate(plant="Papaver")

papoverlap_df_summ <- pap_fc %>%
  filter(taxon != "Papaver") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(pap_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_pap), na.rm = TRUE)) %>%
  ungroup()%>%
  group_by(taxon)%>%
  summarise(mean_overlap=mean(overlap))

papoverlap_df$year <- factor(papoverlap_df$year, levels = sort(as.numeric(levels(papoverlap_df$year))))

papo=ggplot(papoverlap_df, aes(x = year, y = overlap, color = taxon, group = taxon)) +
  geom_smooth(method="lm") +
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

papoverlap_trends <- papoverlap_df %>%mutate(year = as.numeric(year))%>%
  group_by(taxon) %>%
  filter(n_distinct(year) > 1) %>%
  do({
    model <- lm(overlap ~ year, data = .)
    p_val <- summary(model)$coefficients["year", "Pr(>|t|)"]
    data.frame(
      slope = coef(model)[["year"]],
      p_value = p_val,
      r_squared = summary(model)$r.squared,
      significant = ifelse(p_val < 0.05, TRUE, FALSE))})

papov=left_join(papoverlap_trends, papoverlap_df_summ)%>%
  mutate(resource="Papaver")%>%
  rename(consumer=taxon)

#Dryas####
dryoverlap_df<- dry_fc %>%
  filter(taxon != "Dryas") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(dry_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_dry), na.rm = TRUE)) %>%
  ungroup()%>%
  mutate(plant="Dryas")

dryoverlap_df_summ <- dry_fc %>%
  filter(taxon != "Dryas") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(dry_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_dry), na.rm = TRUE)) %>%
  ungroup()%>%
  group_by(taxon)%>%
  summarise(mean_overlap=mean(overlap))

dryoverlap_df$year <- factor(dryoverlap_df$year, levels = sort(as.numeric(levels(dryoverlap_df$year))))

dryo=ggplot(dryoverlap_df, aes(x = year, y = overlap, color = taxon, group = taxon)) +
  geom_smooth(method="lm") +
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
dryoverlap_trends <- dryoverlap_df %>%mutate(year = as.numeric(year))%>%
  group_by(taxon) %>%
  filter(n_distinct(year) > 1) %>%
  do({
    model <- lm(overlap ~ year, data = .)
    p_val <- summary(model)$coefficients["year", "Pr(>|t|)"]
    data.frame(
      slope = coef(model)[["year"]],
      p_value = p_val,
      r_squared = summary(model)$r.squared,
      significant = ifelse(p_val < 0.05, TRUE, FALSE))})

dryov=left_join(dryoverlap_trends, dryoverlap_df_summ)%>%
  mutate(resource="Dryas")%>%
  rename(consumer=taxon)

#Cassiope####
casoverlap_df<- cas_fc %>%
  filter(taxon != "Cassiope") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(cas_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_cas), na.rm = TRUE)) %>%
  ungroup()%>%
  mutate(plant="Cassiope")

casoverlap_df_summ <- cas_fc %>%
  filter(taxon != "Cassiope") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(cas_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_cas), na.rm = TRUE)) %>%
  ungroup()%>%
  group_by(taxon)%>%
  summarise(mean_overlap=mean(overlap))

casoverlap_df$year <- factor(casoverlap_df$year, levels = sort(as.numeric(levels(casoverlap_df$year))))

caso=ggplot(casoverlap_df, aes(x = year, y = overlap, color = taxon, group = taxon)) +
  geom_smooth(method="lm") +
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

casoverlap_trends <- casoverlap_df %>%mutate(year = as.numeric(year))%>%
  group_by(taxon) %>%
  filter(n_distinct(year) > 1) %>%
  do({
    model <- lm(overlap ~ year, data = .)
    p_val <- summary(model)$coefficients["year", "Pr(>|t|)"]
    data.frame(
      slope = coef(model)[["year"]],
      p_value = p_val,
      r_squared = summary(model)$r.squared,
      significant = ifelse(p_val < 0.05, TRUE, FALSE))})

casov=left_join(casoverlap_trends, casoverlap_df_summ)%>%
  mutate(resource="Cassiope")%>%
  rename(consumer=taxon)

#Salix####
saloverlap_df<- sal_fc %>%
  filter(taxon != "Salix") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(sal_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_sal), na.rm = TRUE)) %>%
  ungroup()%>%
  mutate(plant="Salix")

saloverlap_df_summ <- sal_fc %>%
  filter(taxon != "Salix") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(sal_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_sal), na.rm = TRUE)) %>%
  ungroup()%>%
  group_by(taxon)%>%
  summarise(mean_overlap=mean(overlap))

saloverlap_df$year <- factor(saloverlap_df$year, levels = sort(as.numeric(levels(saloverlap_df$year))))

salo=ggplot(saloverlap_df, aes(x = year, y = overlap, color = taxon, group = taxon)) +
  geom_smooth(method="lm") +
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

saloverlap_trends <- saloverlap_df %>%mutate(year = as.numeric(year))%>%
  group_by(taxon) %>%
  filter(n_distinct(year) > 1) %>%
  do({
    model <- lm(overlap ~ year, data = .)
    p_val <- summary(model)$coefficients["year", "Pr(>|t|)"]
    data.frame(
      slope = coef(model)[["year"]],
      p_value = p_val,
      r_squared = summary(model)$r.squared,
      significant = ifelse(p_val < 0.05, TRUE, FALSE))})

salov=left_join(saloverlap_trends, saloverlap_df_summ)%>%
  mutate(resource="Salix")%>%
  rename(consumer=taxon)

#Saxifraga####
saxoverlap_df<- sax_fc %>%
  filter(taxon != "Saxifraga") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(sax_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_sax), na.rm = TRUE)) %>%
  ungroup()%>%
  mutate(plant="Saxifraga")

saxoverlap_df_summ <- sax_fc %>%
  filter(taxon != "Saxifraga") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(sax_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_sax), na.rm = TRUE)) %>%
  ungroup()%>%
  group_by(taxon)%>%
  summarise(mean_overlap=mean(overlap))

saxoverlap_df$year <- factor(saxoverlap_df$year, levels = sort(as.numeric(levels(saxoverlap_df$year))))

saxo=ggplot(saxoverlap_df, aes(x = year, y = overlap, color = taxon, group = taxon)) +
  geom_smooth(method="lm") +
  geom_point(size = 2) +facet_wrap(~taxon)+
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

saxoverlap_trends <- saxoverlap_df %>%mutate(year = as.numeric(year))%>%
  group_by(taxon) %>%
  filter(n_distinct(year) > 1) %>%
  do({
    model <- lm(overlap ~ year, data = .)
    p_val <- summary(model)$coefficients["year", "Pr(>|t|)"]
    data.frame(
      slope = coef(model)[["year"]],
      p_value = p_val,
      r_squared = summary(model)$r.squared,
      significant = ifelse(p_val < 0.05, TRUE, FALSE))})

saxov=left_join(saxoverlap_trends, saxoverlap_df_summ)%>%
  mutate(resource="Saxifraga")%>%
  rename(consumer=taxon)

#all

all_ov=do.call(rbind, list(casov,dryov,papov,silov,saxov,salov))

library(ggplot2)
library(dplyr)

# Prepare the matrix data
matrix_data <- all_ov %>%
  mutate(
    slope_color = case_when(
      significant & slope > 0 ~ "Increasing",
      significant & slope < 0 ~ "Decreasing",
      TRUE ~ NA_character_
    )
  )

# Plot: Use slope_color directly and allow NA to show as grey background
ggplot(matrix_data, aes(x = resource, y = consumer)) +
  geom_tile(aes(fill = slope_color), color = "white", size = 0.5) +
  scale_fill_manual(
    values = c(
      "Increasing" = "yellow",
      "Decreasing" = "blue"
    ),
    na.value = "grey90",  # Light grey for non-significant cells
    name = "Trend"
  ) +
  geom_text(
    data = matrix_data %>% filter(significant),
    aes(label = paste0( round(mean_overlap, 2))),
    size = 3,
    color = "black",
    vjust = 1.5
  ) +
  theme_minimal() +
  labs(
    title = "Trends in Temporal Overlap",
    x = "Resource",
    y = "Consumer"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )







