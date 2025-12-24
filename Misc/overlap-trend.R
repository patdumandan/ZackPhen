library(dplyr)
library(broom)
library(purrr)

# Filter out NAs and ensure year is numeric
clean_df <- siloverlap_df %>%
  filter(!is.na(overlap), !is.na(year)) %>%
  mutate(year = as.numeric(year))

dryoverlap_trends <- dryoverlap_df %>%
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
      significant = ifelse(p_val < 0.05, TRUE, FALSE)
    ) }) %>%ungroup()

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
      significant = ifelse(p_val < 0.05, TRUE, FALSE)
    ) }) %>%ungroup()

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
      significant = ifelse(p_val < 0.05, TRUE, FALSE)
    ) }) %>%ungroup()

#Silene###
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

siloverlap_df_summ <- sil_fc %>%
  filter(taxon != "Silene") %>%
  select(DOY, year, taxon, prob_norm) %>%
  left_join(sil_curves, by = c("DOY", "year")) %>%
  group_by(taxon, year) %>%
  summarise(overlap = sum(pmin(prob_norm, prob_sil), na.rm = TRUE)) %>%
  ungroup()%>%
  group_by(taxon)%>%
  summarise(mean_overlap=mean(overlap))

silov=left_join(siloverlap_trends, siloverlap_df_summ)%>%
  mutate(resource="Silene")%>%
  rename(consumer=taxon)

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
      significant = ifelse(p_val < 0.05, TRUE, FALSE)
    ) }) %>%ungroup()
