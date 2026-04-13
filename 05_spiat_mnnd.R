# MNND Analysis — D-ESMEL Study
# Computes mean nearest-neighbour distance (MNND) between:
#   1. Melanoma cells and each immune cell type
#   2. Immune cell pairs (immune–immune distances)
#   3. Ratio: MNND(Melanoma–CD8) / MNND(Melanoma–Macrophage)
#
# Input:  SPE objects loaded in the R environment (named SPE_*no)
# Output: Excel tables and boxplot figures saved to the working directory
#
# Set your output directory here — all files will be saved to this location
# setwd("your/output/folder")

library(spatstat.geom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)

# =============================================================================
# SHARED FUNCTION: compute MNND between two cell types in a SPE object
# =============================================================================

compute_mnnd <- function(spe_object, ref_type, target_type) {
  coords     <- spatialCoords(spe_object)
  cell_types <- colData(spe_object)$Cell.Type
  ref_coords    <- coords[cell_types == ref_type,    , drop = FALSE]
  target_coords <- coords[cell_types == target_type, , drop = FALSE]
  if (nrow(ref_coords) == 0 || nrow(target_coords) == 0) return(NA_real_)
  win <- owin(range(coords[,1]), range(coords[,2]))
  pp_ref    <- ppp(ref_coords[,1],    ref_coords[,2],    window = win)
  pp_target <- ppp(target_coords[,1], target_coords[,2], window = win)
  mean(nncross(pp_ref, pp_target)$dist, na.rm = TRUE)
}

# Get all SPE objects in the environment
spe_objects <- ls(pattern = "^SPE.*no$")

# =============================================================================
# PART 1: MELANOMA–IMMUNE MNND
# =============================================================================

mnnd_results <- lapply(spe_objects, function(spe_name) {
  spe <- get(spe_name)
  tibble(
    Sample               = spe_name,
    Pair_ID              = sub("^SPE_(\\d+)[P][0-9A-Za-z]*no$", "\\1", spe_name),
    Group                = ifelse(grepl("P0", spe_name), "Control", "Case"),
    Melanoma_Cytotoxic_T = compute_mnnd(spe, "Melanoma", "Cytotoxic T-cell"),
    Melanoma_T_cell      = compute_mnnd(spe, "Melanoma", "T-cell"),
    Melanoma_Macrophage  = compute_mnnd(spe, "Melanoma", "Macrophage"),
    Melanoma_B_cell      = compute_mnnd(spe, "Melanoma", "B-cell lineage")
  )
}) %>% bind_rows()

# Keep only complete 1:1 pairs
valid_pairs <- mnnd_results %>%
  distinct(Pair_ID, Group) %>%
  count(Pair_ID, Group) %>%
  pivot_wider(names_from = Group, values_from = n, values_fill = 0) %>%
  filter(Control == 1, Case == 1) %>%
  pull(Pair_ID)

mnnd_results <- mnnd_results %>%
  filter(Pair_ID %in% valid_pairs) %>%
  mutate(Group = factor(Group, levels = c("Control", "Case")))

# Long format for plotting
mnnd_long <- mnnd_results %>%
  pivot_longer(
    cols = c(Melanoma_Cytotoxic_T, Melanoma_T_cell, Melanoma_Macrophage, Melanoma_B_cell),
    names_to = "Cell_Pair", values_to = "MNND"
  ) %>%
  mutate(Cell_Pair = recode(Cell_Pair,
    Melanoma_Cytotoxic_T = "Melanoma\u2013Cytotoxic T",
    Melanoma_T_cell      = "Melanoma\u2013T cell",
    Melanoma_Macrophage  = "Melanoma\u2013Macrophage",
    Melanoma_B_cell      = "Melanoma\u2013B cell"
  ))

# Paired Wilcoxon per Cell_Pair + BH correction
stats_melanoma <- mnnd_long %>%
  group_by(Cell_Pair) %>%
  summarise(
    p = {
      wide <- pivot_wider(cur_data_all(), id_cols = Pair_ID,
                          names_from = Group, values_from = MNND)
      if (all(c("Case", "Control") %in% names(wide)) && nrow(wide) > 0)
        wilcox.test(wide$Case, wide$Control, paired = TRUE)$p.value
      else NA_real_
    },
    n_pairs = n_distinct(Pair_ID),
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    label = paste0("p=", signif(p, 3), ifelse(!is.na(p.adj) & p.adj < 0.05, "*", ""))
  )

# Summary table and Excel export
mnnd_summary <- mnnd_long %>%
  group_by(Cell_Pair, Group) %>%
  summarise(
    mean_mnnd   = mean(MNND, na.rm = TRUE),
    median_mnnd = median(MNND, na.rm = TRUE),
    n_samples   = n(),
    .groups = "drop"
  )

write.xlsx(
  left_join(mnnd_summary, stats_melanoma, by = "Cell_Pair") %>% arrange(Cell_Pair, Group),
  file = "mnnd_melanoma_immune_stats.xlsx"
)

# Boxplot
p_mnnd <- ggplot(mnnd_long, aes(x = Group, y = MNND, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5) +
  geom_text(data = stats_melanoma, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.3, size = 3) +
  scale_fill_manual(values = c("Case" = "#E69F00", "Control" = "#56B4E9")) +
  facet_wrap(~ Cell_Pair, scales = "free_y", ncol = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = NULL, y = "Mean nearest-neighbour distance (\u00b5m)") +
  theme_bw(base_size = 10) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_blank(),
    strip.text         = element_text(face = "bold", size = 10),
    strip.background   = element_rect(fill = "white"),
    panel.spacing      = unit(1.5, "lines"),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, size = 0.5),
    plot.margin        = unit(c(0.8, 0.5, 0.5, 0.5), "cm")
  )

ggsave("mnnd_melanoma_immune_boxplots.png", p_mnnd, width = 8, height = 6, dpi = 300)
ggsave("mnnd_melanoma_immune_boxplots.pdf", p_mnnd, width = 8, height = 6)

# =============================================================================
# PART 2: IMMUNE–IMMUNE MNND
# =============================================================================

immune_cell_types <- c("Cytotoxic T-cell", "T-cell", "Macrophage", "B-cell lineage")

mnnd_immune_results <- list()
for (spe_name in spe_objects) {
  spe <- get(spe_name)
  for (i in 1:(length(immune_cell_types) - 1)) {
    for (j in (i + 1):length(immune_cell_types)) {
      mnnd_immune_results[[length(mnnd_immune_results) + 1]] <- tibble(
        Sample    = spe_name,
        Pair_ID   = sub("^SPE_(\\d+)[P][0-9A-Za-z]*no$", "\\1", spe_name),
        Group     = ifelse(grepl("P0", spe_name), "Control", "Case"),
        Cell_Pair = paste(immune_cell_types[i], immune_cell_types[j], sep = " \u2013 "),
        MNND      = compute_mnnd(spe, immune_cell_types[i], immune_cell_types[j])
      )
    }
  }
}
mnnd_immune_results <- bind_rows(mnnd_immune_results)

# Keep only complete 1:1 pairs
valid_pairs_immune <- mnnd_immune_results %>%
  distinct(Pair_ID, Group) %>%
  count(Pair_ID, Group) %>%
  pivot_wider(names_from = Group, values_from = n, values_fill = 0) %>%
  filter(Control == 1, Case == 1) %>%
  pull(Pair_ID)

mnnd_immune_results <- mnnd_immune_results %>%
  filter(Pair_ID %in% valid_pairs_immune) %>%
  mutate(Group = factor(Group, levels = c("Control", "Case")))

# Paired Wilcoxon per Cell_Pair + BH correction
stats_immune <- mnnd_immune_results %>%
  group_by(Cell_Pair) %>%
  summarise(
    p = {
      wide <- pivot_wider(cur_data_all(), id_cols = Pair_ID,
                          names_from = Group, values_from = MNND)
      if (all(c("Case", "Control") %in% names(wide)) && nrow(wide) > 0)
        wilcox.test(wide$Case, wide$Control, paired = TRUE)$p.value
      else NA_real_
    },
    n_pairs = n_distinct(Pair_ID),
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    label = paste0("p=", signif(p, 3), ifelse(!is.na(p.adj) & p.adj < 0.05, "*", ""))
  )

# Summary table and Excel export
summary_immune <- mnnd_immune_results %>%
  group_by(Cell_Pair, Group) %>%
  summarise(
    mean_mnnd   = mean(MNND, na.rm = TRUE),
    median_mnnd = median(MNND, na.rm = TRUE),
    n_samples   = n(),
    .groups = "drop"
  )

write.xlsx(
  left_join(summary_immune, stats_immune, by = "Cell_Pair") %>% arrange(Cell_Pair, Group),
  file = "mnnd_immune_immune_stats.xlsx"
)

# Boxplot
p_mnnd_immune <- ggplot(mnnd_immune_results, aes(x = Group, y = MNND, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5) +
  geom_text(data = stats_immune, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.3, size = 3) +
  scale_fill_manual(values = c("Case" = "#E69F00", "Control" = "#56B4E9")) +
  facet_wrap(~ Cell_Pair, scales = "free_y", ncol = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = NULL, y = "Mean nearest-neighbour distance (\u00b5m)") +
  theme_bw(base_size = 10) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_blank(),
    strip.text         = element_text(face = "bold", size = 10),
    strip.background   = element_rect(fill = "white"),
    panel.spacing      = unit(1.5, "lines"),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, size = 0.5),
    plot.margin        = unit(c(0.8, 0.5, 0.5, 0.5), "cm")
  )

ggsave("mnnd_immune_immune_boxplots.png", p_mnnd_immune, width = 8, height = 6, dpi = 300)
ggsave("mnnd_immune_immune_boxplots.pdf", p_mnnd_immune, width = 8, height = 6)

# =============================================================================
# PART 3: RATIO MNND(Melanoma–CD8) / MNND(Melanoma–Macrophage)
# =============================================================================

# Compute ratio per sample (reuses mnnd_results from Part 1)
ratio_df <- mnnd_results %>%
  mutate(ratio_cd8_over_mac = ifelse(Melanoma_Macrophage > 0,
                                     Melanoma_Cytotoxic_T / Melanoma_Macrophage,
                                     NA_real_))

# Keep complete 1:1 pairs with non-NA ratios
valid_pairs_ratio <- ratio_df %>%
  filter(!is.na(ratio_cd8_over_mac)) %>%
  distinct(Pair_ID, Group) %>%
  count(Pair_ID, Group) %>%
  pivot_wider(names_from = Group, values_from = n, values_fill = 0) %>%
  filter(Control == 1, Case == 1) %>%
  pull(Pair_ID)

ratio_df <- ratio_df %>%
  filter(Pair_ID %in% valid_pairs_ratio, !is.na(ratio_cd8_over_mac))

# Paired Wilcoxon
wide_ratio <- ratio_df %>%
  select(Pair_ID, Group, ratio_cd8_over_mac) %>%
  pivot_wider(names_from = Group, values_from = ratio_cd8_over_mac)

p_val_ratio <- if (nrow(wide_ratio) > 0)
  wilcox.test(wide_ratio$Case, wide_ratio$Control, paired = TRUE)$p.value
else NA_real_

stats_ratio <- tibble(
  Metric = "MNND(Mel\u2013CD8)/MNND(Mel\u2013Mac)",
  p      = p_val_ratio,
  p.adj  = p.adjust(p_val_ratio, method = "BH"),
  label  = paste0("p=", signif(p_val_ratio, 3),
                  ifelse(!is.na(p.adjust(p_val_ratio, method = "BH")) &
                           p.adjust(p_val_ratio, method = "BH") < 0.05, "*", ""))
)

# Boxplot
p_ratio <- ggplot(ratio_df, aes(x = Group, y = ratio_cd8_over_mac, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5) +
  geom_text(data = stats_ratio, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.3, size = 3) +
  scale_fill_manual(values = c("Case" = "#E69F00", "Control" = "#56B4E9")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = NULL, y = "MNND(Mel\u2013CD8) / MNND(Mel\u2013Mac)") +
  theme_bw(base_size = 10) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, size = 0.5)
  )

ggsave("mnnd_cd8_mac_ratio_boxplot.png", p_ratio, width = 6, height = 5, dpi = 300)
ggsave("mnnd_cd8_mac_ratio_boxplot.pdf", p_ratio, width = 6, height = 5)

# Summary table and Excel export
summary_ratio <- ratio_df %>%
  group_by(Group) %>%
  summarise(
    mean_ratio   = mean(ratio_cd8_over_mac, na.rm = TRUE),
    median_ratio = median(ratio_cd8_over_mac, na.rm = TRUE),
    n            = n(),
    .groups      = "drop"
  ) %>%
  mutate(Metric = "MNND(Mel\u2013CD8)/MNND(Mel\u2013Mac)") %>%
  left_join(stats_ratio, by = "Metric") %>%
  arrange(Group)

write.xlsx(summary_ratio, "mnnd_cd8_mac_ratio_stats.xlsx")
