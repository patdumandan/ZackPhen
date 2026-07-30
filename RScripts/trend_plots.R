# Plot 1: trend PRECISION
plot_precision <- function(df, taxon_name) {

  df_taxon <- df %>%
    filter(taxon == taxon_name)

  ggplot(df_taxon, aes(x = tsl, y = start_yr, fill = trend_precision)) +
    geom_tile(color = "grey90", linewidth = 0.2) +
    scale_fill_viridis_c(
      option = "magma",
      name = "1 / CI width"
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    coord_equal() +
    theme_classic(base_size = 11) +
    labs(
      x = NULL,
      y = NULL,
      title = taxon_name
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

#Plot 2: trend certainty
plot_certainty <- function(df, taxon_name) {

  df_taxon <- df %>%
    filter(taxon == taxon_name)

  ggplot(df_taxon, aes(x = tsl, y = start_yr, fill = direction_certainty)) +
    geom_tile(color = "grey90", linewidth = 0.2) +
    scale_fill_viridis_c(
      option = "magma",
      limits = c(0.5, 1),
      name = "Certainty",
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    coord_equal() +
    theme_classic(base_size = 11) +
    labs(
      x = NULL,
      y = NULL,
      title = taxon_name) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

# Plot 2: trend directionality
plot_direction <- function(df, taxon_name) {

  df_taxon <- df %>%
    filter(taxon == taxon_name)

  max_abs_slope <- max(abs(df_taxon$slope_mean), na.rm = TRUE)

  ggplot(df_taxon, aes(x = tsl, y = start_yr, fill = slope_mean)) +
    geom_tile(color = "grey90", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-max_abs_slope, max_abs_slope),
      name = "Slope"
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    coord_equal() +
    theme_classic(base_size = 11) +
    labs(
      x = NULL,
      y = NULL,
      title = taxon_name) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

taxa <- unique(slope_df_all$taxon)

precision_plots <- map(taxa,~ plot_precision(slope_df_all, .x))
certainty_plots <- map(taxa,~ plot_certainty(slope_df_all, .x))
direction_plots <- map(taxa,~ plot_direction(slope_df_all, .x))

names(precision_plots) <- taxa
names(certainty_plots) <- taxa
names(direction_plots) <- taxa

plant_taxa <- c("Dryas", "Cassiope", "Silene", "Salix", "Papaver", "Saxifraga")
arth_taxa_1 <- c("Muscidae", "Ichneumonidae", "Chironomidae",
                 "Acari", "Coccoidea", "Collembola")

arth_taxa_2 <- c("Linyphiidae", "Lycosidae", "Nymphalidae",
                 "Phoridae", "Sciaridae")

allplants_precision <- ggarrange(
  plotlist = precision_plots[plant_taxa],
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "right"
)

allplants_certainty <- ggarrange(
  plotlist = certainty_plots[plant_taxa],
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "right"
)

allplants_direction <- ggarrange(
  plotlist = direction_plots[plant_taxa],
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "right"
)

allplants_precision
annotate_figure(
  allplants_precision,
  left = text_grob("Start year", rot = 90, size = 12),
  bottom = text_grob("Time-series length (years)", size = 12))

allplants_certainty
annotate_figure(
  allplants_certainty,
  left = text_grob("Start year", rot = 90, size = 12),
  bottom = text_grob("Time-series length (years)", size = 12))

allplants_direction
annotate_figure(
  allplants_direction,
  left = text_grob("Start year", rot = 90, size = 12),
  bottom = text_grob("Time-series length (years)", size = 12))


allarth1_precision <- ggarrange(
  plotlist = precision_plots[arth_taxa_1],
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "right")
annotate_figure(
  allarth1_precision,
  left = text_grob("Start year", rot = 90, size = 12),
  bottom = text_grob("Time-series length (years)", size = 12))


allarth2_precision <- ggarrange(
  plotlist = precision_plots[arth_taxa_2],
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "right")
annotate_figure(
  allarth2_precision,
  left = text_grob("Start year", rot = 90, size = 12),
  bottom = text_grob("Time-series length (years)", size = 12))

arth1_certainty <- ggarrange(
  plotlist = certainty_plots[arth_taxa_1],
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "right")
annotate_figure(
  arth1_certainty,
  left = text_grob("Start year", rot = 90, size = 12),
  bottom = text_grob("Time-series length (years)", size = 12))

arth2_certainty <- ggarrange(
  plotlist = certainty_plots[arth_taxa_2],
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "right")
annotate_figure(
  arth2_certainty,
  left = text_grob("Start year", rot = 90, size = 12),
  bottom = text_grob("Time-series length (years)", size = 12))

arth1_direction <- ggarrange(
  plotlist = direction_plots[arth_taxa_1],
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "right")
annotate_figure(
  arth1_direction,
  left = text_grob("Start year", rot = 90, size = 12),
  bottom = text_grob("Time-series length (years)", size = 12))

arth2_direction <- ggarrange(
  plotlist = direction_plots[arth_taxa_2],
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "right")
annotate_figure(
  arth2_direction,
  left = text_grob("Start year", rot = 90, size = 12),
  bottom = text_grob("Time-series length (years)", size = 12))

main_taxa <- c("Dryas", "Muscidae", "Silene", "Ichneumonidae", "Papaver", "Chironomidae")

main_direction <- ggarrange(
  plotlist = direction_plots[main_taxa],
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "right")
annotate_figure(
  main_direction,
  left = text_grob("Start year", rot = 90, size = 12),
  bottom = text_grob("Time-series length (years)", size = 12))

main_precision <- ggarrange(
  plotlist = precision_plots[main_taxa],
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "right")
annotate_figure(
  main_precision,
  left = text_grob("Start year", rot = 90, size = 12),
  bottom = text_grob("Time-series length (years)", size = 12))

main_certainty <- ggarrange(
  plotlist = certainty_plots[main_taxa],
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "right")
annotate_figure(
  main_certainty,
  left = text_grob("Start year", rot = 90, size = 12),
  bottom = text_grob("Time-series length (years)", size = 12))
