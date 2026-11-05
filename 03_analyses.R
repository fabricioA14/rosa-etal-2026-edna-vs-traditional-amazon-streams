# =============================================================================
# Comparative Performance of Traditional Sampling and Water Environmental DNA
# for Detecting Fish Diversity in Amazonian Streams Across Hydrological Periods
#
# Authors: Fabricio dos Anjos Santa Rosa, Derlan José Ferreira Silva,
#          Cintia Oliveira Carvalho, Silvia Britto Barreto,
#          Silvana Melo Sviggum, Anna Karolina Oliveira de Queiroz,
#          Enki Gautier, Marcia da Silva Anjos,
#          Jessica Conceição Nascimento Dergan,
#          Orlando Gabriel da Silva Solano Reis,
#          Tibério Cesar Tortola Burlamaqui, Birgitte Lisbeth Graae Thorbek,
#          João Bráullio de Luna Sales, Quentin Mauvisseau, Hugo de Boer,
#          Jonathan Stuart Ready
# Contact: Fabricio dos Anjos Santa Rosa (fabriciodosanjos14@gmail.com)
#
# Script:  03_analyses.R
# Purpose: Statistical analyses and figures (step 3): taxon counts, heatmap
#          (Figure 2), nMDS (Figure 3), ANOSIM/PERMANOVA/betadisper (Figure 4)
#          and SIMPER (Figure 5).
# Input:   SuppInfo_S1.xlsx, sheet 3_analysis_input
#          fish_icons/ (optional order illustrations for Figure 2)
# Output:  stats_82taxa.txt, simper_results_82taxa.csv and figures (JPG, 300
#          dpi)
# Notes:   Period: August = Dry, February = Wet.
# =============================================================================

# 1. Setup --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
if (length(script_dir) == 1 && nzchar(script_dir)) setwd(script_dir)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(ggrepel)
  library(stringr)
  library(purrr)
  library(patchwork)
})

data <- as.data.frame(readxl::read_excel("SuppInfo_S1.xlsx", sheet = "3_analysis_input"))

data <- data[, !colnames(data) %in% "N"]

data <- data %>%
  mutate(
    Presence = 1,
    Period = ifelse(substr(date, 1, 2) == "08", "Dry", "Wet"),
    Campaign = factor(
      case_when(
        date == "08_22" ~ "August 2022",
        date == "08_23" ~ "August 2023",
        date == "02_23" ~ "February 2023",
        date == "02_24" ~ "February 2024"
      ),
      levels = c("August 2022", "August 2023", "February 2023", "February 2024")
    )
  )

data_unique <- data[!duplicated(data[, c("Sampling_Method", "Locality", "date", "Identification")]), ]

save_fig <- function(plot, name, width, height) {
  ggsave(paste0(name, ".jpg"), plot, width = width, height = height, dpi = 300, bg = "white")
}

stats_file <- "stats_82taxa.txt"
cat("", file = stats_file)
log_out <- function(title, x) {
  cat("\n\n==== ", title, " ====\n", file = stats_file, append = TRUE, sep = "")
  capture.output(x, file = stats_file, append = TRUE)
}

method_colors <- c("eDNA" = "#377EB8", "Traditional" = "#4DAF4A")
# Method order: traditional, then eDNA
period_colors <- c(
  "Traditional_Dry" = "#238B45",
  "Traditional_Wet" = "#74C476",
  "eDNA_Dry" = "#2171B5",
  "eDNA_Wet" = "#6BAED6"
)


# 2. Taxon counts -------------------------------------------------------------

taxon_methods <- data_unique %>%
  group_by(Identification) %>%
  summarise(
    Methods = paste(sort(unique(Sampling_Method)), collapse = "_"),
    .groups = "drop"
  )

log_out("Number of taxa per method", table(data_unique %>% distinct(Sampling_Method, Identification) %>% pull(Sampling_Method)))
log_out("Exclusive/shared taxa", table(taxon_methods$Methods))
log_out("Total number of taxa", nrow(taxon_methods))


# 3. Heatmap (Figure 2) -------------------------------------------------------

desired_orders <- c(
  "Characiformes", "Siluriformes", "Gymnotiformes", "Cichliformes",
  "Cyprinodontiformes", "Beloniformes", "Synbranchiformes",
  "Tetraodontiformes", "Myliobatiformes"
)

order_colors <- c(
  "Characiformes" = "#F4A6A6", "Siluriformes" = "#F2A6E8",
  "Gymnotiformes" = "#B3B3F2", "Cichliformes" = "#F7D6A6",
  "Cyprinodontiformes" = "#EEF2A6", "Beloniformes" = "#A6D8F2",
  "Synbranchiformes" = "#A6F2B0", "Tetraodontiformes" = "#F2A6C8",
  "Myliobatiformes" = "#C8C8C8"
)

# Family order within each order
family_levels <- c(
  "Acestrorhynchidae", "Acestrorhamphidae", "Characidae", "Iguanodectidae",
  "Chilodontidae", "Erythrinidae", "Lebiasinidae", "Serrasalmidae",
  "Hemiodontidae", "Curimatidae", "Anostomidae", "Crenuchidae",
  "Loricariidae", "Callichthyidae", "Trichomycteridae", "Cetopsidae",
  "Heptapteridae", "Pimelodidae", "Pseudopimelodidae", "Aspredinidae",
  "Auchenipteridae", "Doradidae",
  "Gymnotidae", "Rhamphichthyidae", "Hypopomidae", "Apteronotidae", "Sternopygidae"
)

taxon_order <- data_unique %>%
  distinct(Identification, order, family, genus) %>%
  mutate(
    order = factor(order, levels = desired_orders),
    family_rank = match(family, family_levels, nomatch = 99),
    is_sp = grepl(" sp\\.?$", Identification)
  ) %>%
  arrange(order, family_rank, family, genus, is_sp, Identification)

heatmap_data <- data_unique %>%
  group_by(Period, Campaign, Sampling_Method, Identification, order) %>%
  summarise(Count = n_distinct(Locality), .groups = "drop") %>%
  mutate(
    Identification = factor(Identification, levels = rev(taxon_order$Identification)),
    order = factor(order, levels = desired_orders),
    Period = factor(Period, levels = c("Dry", "Wet")),
    Sampling_Method = factor(Sampling_Method, levels = c("Traditional", "eDNA"))
  )

present_orders <- desired_orders[desired_orders %in% heatmap_data$order]

p_heatmap <- ggplot(
  heatmap_data,
  aes(x = Sampling_Method, y = Identification, fill = Count)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  # White text on dark cells, dark grey on light cells
  geom_text(
    aes(label = Count, color = Count >= 7),
    size = 3.3,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey15")) +
  scale_x_discrete(labels = c("Traditional" = "Trad.", "eDNA" = "eDNA")) +
  scale_fill_gradient(
    low = "#D9F0D3",
    high = "#00441B",
    breaks = seq(0, 12, 2),
    limits = c(0, 12)
  ) +
  guides(fill = guide_colourbar(barwidth = 12, barheight = 0.8, title.position = "top", title.hjust = 0.5)) +
  labs(
    title = "Taxa Occurrences Across Seasons and Sampling Methods in Stream Fish Communities",
    x = "Sampling method",
    y = "Taxa",
    fill = "Frequency of occurrence\n(number of streams)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.08, "in"),
    panel.spacing.x = unit(c(0.08, 0.45, 0.08), "in"),
    strip.placement = "outside",
    strip.background.x = element_blank(),
    strip.text.x = element_text(size = 11, margin = margin(4, 0, 4, 0)),
    strip.text.y.left = element_blank(),
    strip.background.y = element_blank(),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.title.x = element_text(size = 11, margin = margin(t = 8)),
    axis.text.y = element_text(face = "italic", size = 9, colour = "grey10"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, margin = margin(b = 10)),
    plot.title.position = "plot",
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(colour = "grey15"),
    legend.text = element_text(colour = "grey15"),
    legend.box.spacing = unit(0, "pt")
  )

# Nested strips (period > campaign); orders shown only in the legend
p_heatmap <- p_heatmap +
  ggh4x::facet_nested(
    order ~ Period + Campaign,
    scales = "free",
    space = "free",
    switch = "y",
    nest_line = element_line(colour = "grey55", linewidth = 0.4),
    strip = ggh4x::strip_nested(
      text_x = ggh4x::elem_list_text(face = c("bold", "plain"), size = c(13, 10.5), colour = c("grey15", "grey15")),
      by_layer_x = TRUE
    )
  )

g_heatmap <- ggplotGrob(p_heatmap)
lay <- g_heatmap$layout
axis_rows <- which(grepl("^axis-l", lay$name))
axis_left <- min(lay$l[axis_rows])
axis_right <- min(lay$r[axis_rows])  # Axis of the species-name column
panel_idx <- grep("^panel", lay$name)
panel_right <- max(lay$r[panel_idx])
first_col_panels <- lay[panel_idx, ]
first_col_panels <- first_col_panels[first_col_panels$l == min(first_col_panels$l), ]
first_col_panels <- first_col_panels[order(first_col_panels$t), ]

for (k in seq_len(nrow(first_col_panels))) {
  col <- order_colors[[present_orders[k]]]
  row_t <- first_col_panels$t[k]
  row_b <- first_col_panels$b[k]

  # Colour gradient behind species names
  fade <- grid::linearGradient(
    colours = c(grDevices::adjustcolor(col, alpha.f = 0), col),
    x1 = 0, y1 = 0.5, x2 = 1, y2 = 0.5
  )
  band <- grid::rectGrob(gp = grid::gpar(col = NA, fill = fade))
  g_heatmap <- gtable::gtable_add_grob(
    g_heatmap, band,
    t = row_t, b = row_b, l = axis_left, r = axis_right,
    z = 0, clip = "off", name = paste0("order-band-", k)
  )

  # Order-coloured lines at the top and bottom of each block
  edge_gap <- unit(0.025, "in")
  edge_lwd <- unit(0.9, "pt")
  edge_y <- grid::unit.c(unit(1, "npc") + edge_gap, unit(0, "npc") - edge_gap)
  edges_fade <- grid::gTree(children = grid::gList(
    grid::rectGrob(y = edge_y[1], height = edge_lwd, gp = grid::gpar(col = NA, fill = fade)),
    grid::rectGrob(y = edge_y[2], height = edge_lwd, gp = grid::gpar(col = NA, fill = fade))
  ))
  g_heatmap <- gtable::gtable_add_grob(
    g_heatmap, edges_fade,
    t = row_t, b = row_b, l = axis_left, r = axis_right,
    z = 20, clip = "off", name = paste0("order-edges-fade-", k)
  )
  edges_solid <- grid::segmentsGrob(
    x0 = unit(c(0, 0), "npc"), x1 = unit(c(1, 1), "npc"),
    y0 = edge_y, y1 = edge_y,
    gp = grid::gpar(col = col, lwd = 0.9, lineend = "butt")
  )
  g_heatmap <- gtable::gtable_add_grob(
    g_heatmap, edges_solid,
    t = row_t, b = row_b, l = axis_right + 1, r = panel_right,
    z = 20, clip = "off", name = paste0("order-edges-", k)
  )
}

# Frequency scale and order legend in one row
guide_i <- which(g_heatmap$layout$name == "guide-box-bottom")
freq_guide <- g_heatmap$grobs[[guide_i]]
guide_row <- g_heatmap$layout$t[guide_i]
g_heatmap$grobs[[guide_i]] <- ggplot2::zeroGrob()
g_heatmap$heights[guide_row] <- unit(0, "pt")

leg_ncol <- 5
leg_cell_w <- 1.55
leg_cell_h <- 0.3
leg_title_w <- 0.75
leg_nrow <- ceiling(length(present_orders) / leg_ncol)
orders_w <- leg_title_w + leg_ncol * leg_cell_w
leg_h <- max(leg_nrow * leg_cell_h, 0.7) + 0.3

leg_children <- list(grid::textGrob(
  "Orders", x = unit(0, "npc"), y = unit(0.5, "npc") + unit((leg_nrow - 1) * leg_cell_h / 2, "in"),
  just = c("left", "centre"), gp = grid::gpar(fontsize = 12, fontface = "bold", col = "grey15")
))
for (k in seq_along(present_orders)) {
  cc <- (k - 1) %% leg_ncol
  rr <- (k - 1) %/% leg_ncol
  x0 <- unit(leg_title_w + cc * leg_cell_w, "in")
  yy <- unit(0.5, "npc") + unit(((leg_nrow - 1) / 2 - rr) * leg_cell_h, "in")
  leg_children[[length(leg_children) + 1]] <- grid::rectGrob(
    x = x0, y = yy, width = unit(0.2, "in"), height = unit(0.2, "in"), just = c("left", "centre"),
    gp = grid::gpar(fill = order_colors[[present_orders[k]]], col = NA)
  )
  leg_children[[length(leg_children) + 1]] <- grid::textGrob(
    present_orders[k], x = x0 + unit(0.27, "in"), y = yy, just = c("left", "centre"),
    gp = grid::gpar(fontsize = 10.5, col = "grey15")
  )
}
orders_legend <- grid::gTree(children = do.call(grid::gList, leg_children))

legend_row <- gtable::gtable(
  widths = grid::unit.c(unit(1, "null"), gtable::gtable_width(freq_guide), unit(0.5, "in"),
                        unit(orders_w, "in"), unit(1, "null")),
  heights = unit(1, "null")
)
legend_row <- gtable::gtable_add_grob(legend_row, freq_guide, t = 1, l = 2, clip = "off", name = "freq")
legend_row <- gtable::gtable_add_grob(legend_row, orders_legend, t = 1, l = 4, clip = "off", name = "orders")

g_heatmap <- gtable::gtable_add_rows(g_heatmap, unit(leg_h, "in"), pos = -1)
g_heatmap <- gtable::gtable_add_grob(
  g_heatmap, legend_row,
  t = nrow(g_heatmap), l = 1, r = ncol(g_heatmap), clip = "off", name = "legends"
)

# Order illustrations (fish_icons/) linked to each order block
heat_w <- 10
fish_col_w <- 1.9
heat_h <- 20

g_heatmap <- gtable::gtable_add_cols(g_heatmap, unit(fish_col_w, "in"), pos = panel_right)
fish_col <- panel_right + 1
panel_idx <- grep("^panel", g_heatmap$layout$name)
first_col_panels <- g_heatmap$layout[panel_idx, ]
first_col_panels <- first_col_panels[first_col_panels$l == min(first_col_panels$l), ]
first_col_panels <- first_col_panels[order(first_col_panels$t), ]
rows_t <- min(first_col_panels$t)
rows_b <- max(first_col_panels$b)
g_heatmap <- gtable::gtable_add_grob(
  g_heatmap, grid::rectGrob(gp = grid::gpar(col = NA, fill = NA)), t = rows_t, b = rows_b, l = fish_col,
  name = "fish-column", clip = "off"
)

# Positions (inches) of order blocks and illustration column
vp_extent <- function(vp_name) {
  grid::seekViewport(vp_name)
  lo <- grid::deviceLoc(unit(0, "npc"), unit(0, "npc"))
  hi <- grid::deviceLoc(unit(1, "npc"), unit(1, "npc"))
  c(x0 = as.numeric(lo$x), y0 = as.numeric(lo$y), x1 = as.numeric(hi$x), y1 = as.numeric(hi$y))
}
vp_name_of <- function(i) {
  l <- g_heatmap$layout[i, ]
  paste0(l$name, ".", paste(l$t, l$l, l$b, l$r, sep = "-"))
}
grDevices::pdf(NULL, width = heat_w + fish_col_w, height = heat_h + leg_h)
grid::grid.newpage()
grid::grid.draw(g_heatmap)
col_ext <- vp_extent(vp_name_of(which(g_heatmap$layout$name == "fish-column")))
block_ext <- t(sapply(rownames(first_col_panels), function(i) vp_extent(vp_name_of(as.integer(i)))))
grDevices::dev.off()

block_df <- data.frame(
  order = present_orders,
  ytop = block_ext[, "y1"] - col_ext["y0"],
  ybot = block_ext[, "y0"] - col_ext["y0"]
)

icon_dir <- "fish_icons"
icon_w <- 1.25
icon_gap <- 0.06
# Illustrations are optional: an order without a PNG in fish_icons/ is drawn without it
icons <- lapply(block_df$order, function(o) {
  f <- file.path(icon_dir, paste0(o, ".png"))
  if (file.exists(f)) png::readPNG(f) else NULL
})
if (any(sapply(icons, is.null))) {
  message("fish_icons/: no illustration for ", paste(block_df$order[sapply(icons, is.null)], collapse = ", "))
}
block_df$icon_h <- sapply(icons, function(im) if (is.null(im)) 0.62 else min(icon_w * nrow(im) / ncol(im), 0.62))
block_df$icon_w <- sapply(seq_along(icons), function(k) {
  if (is.null(icons[[k]])) icon_w else block_df$icon_h[k] * ncol(icons[[k]]) / nrow(icons[[k]])
})

# Centre each illustration on its block without overlaps
target <- (block_df$ytop + block_df$ybot) / 2
col_h <- col_ext["y1"] - col_ext["y0"]
yc <- target
for (k in seq_along(yc)) {
  top_limit <- if (k == 1) col_h else yc[k - 1] - block_df$icon_h[k - 1] / 2 - icon_gap
  yc[k] <- min(yc[k], top_limit - block_df$icon_h[k] / 2)
}
for (k in rev(seq_along(yc))) {
  bottom_limit <- if (k == length(yc)) 0 else yc[k + 1] + block_df$icon_h[k + 1] / 2 + icon_gap
  yc[k] <- max(yc[k], bottom_limit + block_df$icon_h[k] / 2)
}
block_df$yc <- yc

fan_len <- fish_col_w - icon_w - 0.08
fish_grobs <- list()
for (k in seq_len(nrow(block_df))) {
  b <- block_df[k, ]
  col <- order_colors[[b$order]]
  fish_grobs[[length(fish_grobs) + 1]] <- grid::polygonGrob(
    x = unit(c(0, fan_len, fan_len, 0), "in"),
    y = unit(c(b$ytop, b$yc + b$icon_h * 0.35, b$yc - b$icon_h * 0.35, b$ybot), "in"),
    gp = grid::gpar(col = NA, fill = grid::linearGradient(
      colours = c(col, grDevices::adjustcolor(col, alpha.f = 0.15)),
      x1 = 0, y1 = 0.5, x2 = 1, y2 = 0.5
    ))
  )
  if (is.null(icons[[k]])) next
  fish_grobs[[length(fish_grobs) + 1]] <- grid::rasterGrob(
    icons[[k]],
    x = unit(fan_len + 0.04, "in"), y = unit(b$yc, "in"),
    width = unit(b$icon_w, "in"), height = unit(b$icon_h, "in"),
    just = c("left", "centre"), interpolate = TRUE
  )
}
fish_tree <- grid::gTree(children = do.call(grid::gList, fish_grobs))
g_heatmap$grobs[[which(g_heatmap$layout$name == "fish-column")]] <- fish_tree

ggsave("heatmap_82taxa.jpg", g_heatmap, width = heat_w + fish_col_w, height = heat_h + leg_h, dpi = 300, bg = "white")


# 4. nMDS (Figure 3) ----------------------------------------------------------

nmds_data <- data_unique %>%
  select(Sampling_Method, Locality, date, Identification, Presence) %>%
  pivot_wider(
    names_from = Identification,
    values_from = Presence,
    values_fill = 0
  )

nmds_metadata <- nmds_data[, 1:3]
nmds_matrix <- nmds_data[, 4:ncol(nmds_data)]

set.seed(123)
nmds_result <- metaMDS(nmds_matrix, distance = "jaccard", trymax = 50, trace = 0)
nmds_stress <- nmds_result$stress
log_out("NMDS stress", nmds_stress)

nmds_sites <- data.frame(
  nmds_metadata,
  MDS1 = nmds_result$points[, 1],
  MDS2 = nmds_result$points[, 2]
) %>%
  mutate(
    Period = ifelse(substr(date, 1, 2) == "08", "Dry", "Wet"),
    Method_Period = paste(Sampling_Method, Period, sep = "_")
  )

abbreviate_taxon <- function(x) {
  tokens <- str_split(x, "\\s+")[[1]]
  if (any(str_detect(tokens, "^sp\\.?$")) || length(tokens) < 2) return(x)
  paste0(substr(tokens[1], 1, 1), ". ", tokens[2])
}

nmds_taxa <- as.data.frame(scores(nmds_result, display = "species"))
colnames(nmds_taxa)[1:2] <- c("MDS1", "MDS2")
nmds_taxa$Taxon <- rownames(nmds_taxa)
nmds_taxa <- nmds_taxa %>%
  mutate(Taxon_Label = map_chr(Taxon, abbreviate_taxon)) %>%
  left_join(taxon_methods, by = c("Taxon" = "Identification")) %>%
  mutate(
    Method_Category = case_when(
      Methods == "eDNA" ~ "Exclusive_eDNA",
      Methods == "Traditional" ~ "Exclusive_Traditional",
      TRUE ~ "Shared"
    )
  )

n_excl_edna <- sum(nmds_taxa$Method_Category == "Exclusive_eDNA")
n_excl_trad <- sum(nmds_taxa$Method_Category == "Exclusive_Traditional")
n_shared <- sum(nmds_taxa$Method_Category == "Shared")

nmds_hulls_period <- nmds_sites %>%
  group_by(Method_Period) %>%
  slice(chull(MDS1, MDS2)) %>%
  ungroup()

all_scores <- rbind(nmds_result$points, as.matrix(nmds_taxa[, c("MDS1", "MDS2")]))
x_limits <- range(all_scores[, 1])
y_limits <- range(all_scores[, 2])

nmds_panel <- function(taxa, subtitle, tag) {
  ggplot() +
    geom_polygon(
      data = nmds_hulls_period,
      aes(MDS1, MDS2, group = Method_Period, fill = Method_Period),
      alpha = 0.25, color = NA
    ) +
    geom_point(
      data = nmds_sites,
      aes(MDS1, MDS2, fill = Method_Period),
      shape = 21, color = "grey20", size = 3.5, stroke = 1
    ) +
    geom_point(data = taxa, aes(MDS1, MDS2), shape = 17, color = "grey45", size = 1.8, alpha = 0.7) +
    geom_text_repel(
      data = taxa,
      aes(MDS1, MDS2, label = Taxon_Label),
      fontface = "italic", size = 3.2, color = "grey40",
      max.overlaps = Inf, box.padding = 0.4, point.padding = 0.1,
      min.segment.length = 0, segment.color = "grey70", force = 4, seed = 1,
      max.time = Inf  # no time limit, so label positions are reproducible
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    scale_fill_manual(values = period_colors, breaks = names(period_colors)) +
    labs(subtitle = subtitle, tag = tag, x = "nMDS1", y = "nMDS2", fill = "Method + period") +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey50"),
      plot.subtitle = element_text(hjust = 0.5, size = 13, color = "grey30"),
      plot.tag = element_text(size = 20),
      plot.tag.position = c(0.06, 0.93),
      legend.position = "bottom"
    )
}

p_excl <- nmds_panel(
  filter(nmds_taxa, Method_Category != "Shared"),
  paste0("Exclusive taxa: Traditional = ", n_excl_trad, " | eDNA = ", n_excl_edna),
  "A"
)
p_shared <- nmds_panel(
  filter(nmds_taxa, Method_Category == "Shared"),
  paste0("Shared taxa: ", n_shared),
  "B"
)

p_nmds_vertical <- (p_excl / p_shared) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "nMDS – Grouping by Period (Dry/Wet) and Sampling Method",
    caption = paste0("Jaccard distance; stress = ", round(nmds_stress, 3)),
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 17))
  ) &
  theme(legend.position = "bottom")

save_fig(p_nmds_vertical, "nmds_vertical_82taxa", 12, 15.6)
save_fig(p_excl + labs(tag = NULL), "nmds_exclusive_82taxa", 12, 8)
save_fig(p_shared + labs(tag = NULL), "nmds_shared_82taxa", 12, 8)


# 5. ANOSIM, PERMANOVA and betadisper (Figure 4) ------------------------------

community_metadata <- nmds_sites %>%
  transmute(
    Sampling_Method = factor(Sampling_Method),
    Locality = Locality,
    Stream = substr(Locality, 1, 3),
    Period = factor(Period),
    Method_Period = interaction(Sampling_Method, Period, sep = ".")
  )
community_matrix <- as.matrix(nmds_matrix)
community_distance <- vegdist(community_matrix, method = "jaccard")

set.seed(123)
anosim_method <- anosim(community_distance, community_metadata$Sampling_Method, permutations = 999)
set.seed(123)
anosim_method_period <- anosim(community_distance, community_metadata$Method_Period, permutations = 999)
log_out("ANOSIM by method", anosim_method)
log_out("ANOSIM by method x period", anosim_method_period)

anosim_df <- data.frame(
  Dissimilarity = anosim_method_period$dis.rank,
  Comparison = factor(
    as.character(anosim_method_period$class.vec),
    levels = c("Between", "Traditional.Dry", "Traditional.Wet", "eDNA.Dry", "eDNA.Wet")
  )
)

p_anosim <- ggplot(anosim_df, aes(Comparison, Dissimilarity, fill = Comparison)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = "grey80", linewidth = 1) +
  geom_boxplot(width = 0.10, outlier.shape = NA, alpha = 0.8, color = "grey40") +
  scale_fill_manual(values = c(
    "Between" = "#E67155",
    "eDNA.Dry" = "#2171B5", "eDNA.Wet" = "#6BAED6",
    "Traditional.Dry" = "#238B45", "Traditional.Wet" = "#74C476"
  )) +
  labs(
    title = "Differences in Fish Community Composition Across Methods and Periods",
    subtitle = paste0(
      "ANOSIM: R = ", round(anosim_method_period$statistic, 3),
      " | p = ", anosim_method_period$signif
    ),
    x = "Comparison",
    y = "Ranked dissimilarity"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 13, color = "grey30")
  )

save_fig(p_anosim, "anosim_82taxa", 9.9, 6.1)

set.seed(123)
permanova_result <- adonis2(
  community_matrix ~ Sampling_Method * Period,
  data = community_metadata,
  method = "jaccard",
  by = "terms",
  permutations = 999
)
log_out("PERMANOVA (adonis2, Sampling_Method * Period, sequential terms)", permanova_result)

set.seed(123)
permanova_strata <- adonis2(
  community_matrix ~ Sampling_Method * Period,
  data = community_metadata,
  method = "jaccard",
  by = "terms",
  permutations = how(blocks = community_metadata$Stream, nperm = 999)
)
log_out("PERMANOVA with stream as block (strata)", permanova_strata)

betadisper_result <- betadisper(community_distance, community_metadata$Method_Period)
set.seed(123)
log_out("BETADISPER", betadisper_result)
log_out("BETADISPER ANOVA", anova(betadisper_result))
log_out("BETADISPER permutest", permutest(betadisper_result, permutations = 999))
log_out("BETADISPER TukeyHSD", TukeyHSD(betadisper_result))



# 6. SIMPER (Figure 5) --------------------------------------------------------

set.seed(123)
simper_method_period <- simper(community_matrix, group = community_metadata$Method_Period, permutations = 999)

simper_results <- bind_rows(lapply(names(simper_method_period), function(cmp) {
  df <- as.data.frame(simper_method_period[[cmp]])
  df$Taxon <- rownames(df)
  df$Comparison <- cmp
  df
})) %>%
  select(Comparison, Taxon, everything())

write.csv(simper_results, "simper_results_82taxa.csv", row.names = FALSE)

simper_comparison_order <- c(
  "Traditional.Dry_Traditional.Wet",
  "eDNA.Wet_eDNA.Dry",
  "Traditional.Dry_eDNA.Dry",
  "Traditional.Wet_eDNA.Wet"
)

# Comparison names follow the factor-level order
stopifnot(all(simper_comparison_order %in% simper_results$Comparison))

simper_plot_data <- simper_results %>%
  filter(!is.na(p) & p < 0.05, Comparison %in% simper_comparison_order) %>%
  left_join(taxon_methods, by = c("Taxon" = "Identification")) %>%
  mutate(
    Category = factor(recode(Methods, "eDNA_Traditional" = "Shared"), levels = c("Traditional", "Shared", "eDNA")),
    Comparison = factor(Comparison, levels = simper_comparison_order)
  )

log_out("SIMPER: number of significant taxa per comparison and category",
        table(simper_plot_data$Comparison, simper_plot_data$Category))

# Sum contributions per category (one block per bar)
simper_bar_data <- simper_plot_data %>%
  group_by(Comparison, Category) %>%
  summarise(average = sum(average), .groups = "drop")

p_simper <- ggplot(simper_bar_data, aes(Comparison, average, fill = Category)) +
  geom_col(color = NA) +
  scale_fill_manual(values = c("Traditional" = "#74C476", "Shared" = "#E68373", "eDNA" = "#6BAED6")) +
  scale_x_discrete(labels = function(x) gsub("_", "\nvs\n", x)) +
  labs(
    title = "Taxon Contributions to Dissimilarity Between Groups",
    subtitle = "SIMPER analysis by traditional sampling, shared and eDNA taxa (p < 0.05)",
    x = "Comparison",
    y = "SIMPER contribution",
    fill = "Category"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30")
  )

save_fig(p_simper, "simper_82taxa", 9, 6.7)

cat("Done. Results in", getwd(), "\n")
