# Load required libraries
library(SPIAT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)
library(purrr)

# Initialize an empty dataframe to store all p_cells
combined_p_cells <- data.frame()

# Get all spatial objects starting with "SPE" and ending with "no" (cases and controls)
spe_objects <- ls(pattern = "^SPE.*no$")

# Loop through each SPE object, calculate proportions, and combine the results
for (spe_obj in spe_objects) {
  spe_data <- get(spe_obj)
  
  p_cells <- calculate_cell_proportions(
    spe_object = spe_data,
    reference_celltypes = NULL,
    feature_colname = "Cell.Type",
    celltypes_to_exclude = "Others",
    plot.image = FALSE
  )
  
  p_cells$Sample <- spe_obj
  combined_p_cells <- bind_rows(combined_p_cells, p_cells)
}

# Save the combined p_cells to an Excel file
excel_file_path <- "cell_proportions_results.xlsx"
write.xlsx(combined_p_cells, file = excel_file_path)

# Group and Pair_ID
combined_p_cells <- combined_p_cells %>%
  mutate(
    Group = ifelse(grepl("P0", Sample), "Control", "Case"),
    Group = factor(Group, levels = c("Control", "Case")),
    Pair_ID = sub("^SPE_(\\d+)[P][0-9A-Za-z]*no$", "\\1", Sample)
  )

# Keep only complete 1:1 pairs (exactly one Control and one Case)
valid_pairs <- combined_p_cells %>%
  distinct(Sample, Group, Pair_ID) %>%
  count(Pair_ID, Group) %>%
  pivot_wider(names_from = Group, values_from = n, values_fill = 0) %>%
  filter(Control == 1, Case == 1) %>%
  pull(Pair_ID)

combined_p_cells <- combined_p_cells %>%
  filter(Pair_ID %in% valid_pairs)

# Keep immune cells (exclude melanoma)
immune_cells <- combined_p_cells %>%
  filter(Cell_type != "Melanoma")

# ensure group order
immune_cells <- immune_cells %>%
  mutate(Group = factor(Group, levels = c("Control", "Case")))

# complete missing combinations (Pair_ID x Group x Cell_type); absent types -> proportion 0
immune_cells_complete <- immune_cells %>%
  group_by(Pair_ID, Group, Cell_type) %>%
  summarise(Proportion = sum(Proportion), .groups = "drop") %>%   # in case of duplicates
  group_by(Pair_ID, Cell_type) %>%
  complete(Group = factor(c("Control","Case"), levels = c("Control","Case")),
           fill = list(Proportion = 0)) %>%
  ungroup()

# wide for paired test
wide <- immune_cells_complete %>%
  pivot_wider(names_from = Group, values_from = Proportion, values_fill = 0)

# paired Wilcoxon per cell type, FDR adjust, make labels
stats <- wide %>%
  group_by(Cell_type) %>%
  summarise(
    p = {
      x <- Case; y <- Control
      if (length(x) == length(y) && length(x) > 0) {
        wilcox.test(x, y, paired = TRUE)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "fdr"),
    label = paste0("p=", signif(p, 3), ifelse(!is.na(p.adj) & p.adj < 0.05, "*", ""))
  )

# plot (you can keep your existing plot code; just point geom_text to `stats`)
ggplot(immune_cells, aes(x = Group, y = Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5) +
  geom_text(
    data = stats,
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE, vjust = 1.3, size = 3
  ) +
  scale_fill_manual(values = c("Case" = "#E69F00", "Control" = "#56B4E9")) +
  facet_wrap(~ Cell_type, scales = "free_y", ncol = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = NULL, y = "Cell Proportion") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "white"),
    panel.spacing = unit(1.5, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.8, 0.5, 0.5, 0.5), "cm"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# Save the plot to file
plot_file <- "immune_cell_proportions.png"

# Set your output directory here — plots and tables will be saved to this location
# setwd("your/output/folder")
ggsave(
  filename = plot_file,
  plot = last_plot(),   # or replace with your plot object if you saved it as p <- ggplot(...)
  width = 8,
  height = 6,
  dpi = 300
)


# Tumor vs Non-Tumor immune proportions with paired case–control tests and BH correction

# ---- collect per-area proportions ----
all_region_cells <- data.frame()
spe_objects <- ls(pattern = "^SPE.*no$")

for (spe_obj in spe_objects) {
  spe <- get(spe_obj)
  
  if ("Region" %in% colnames(colData(spe))) {
    regions <- unique(colData(spe)$Region)
    for (region in regions) {
      region_spe <- spe[, colData(spe)$Region == region]
      if (ncol(region_spe) == 0) next
      
      pc <- calculate_cell_proportions(
        spe_object = region_spe,
        reference_celltypes = NULL,
        feature_colname = "Cell.Type",
        celltypes_to_exclude = "Others",
        plot.image = FALSE
      )
      pc$Sample <- spe_obj
      pc$Region <- region
      all_region_cells <- bind_rows(all_region_cells, pc)
    }
  } else {
    pc <- calculate_cell_proportions(
      spe_object = spe,
      reference_celltypes = NULL,
      feature_colname = "Cell.Type",
      celltypes_to_exclude = "Others",
      plot.image = FALSE
    )
    pc$Sample <- spe_obj
    pc$Region <- "Whole_Slide"
    all_region_cells <- bind_rows(all_region_cells, pc)
  }
}

# standardize and keep areas of interest
std_area <- function(x) {
  lx <- tolower(x)
  ifelse(lx %in% c("tumor", "tumour", "tum"), "Tumor",
         ifelse(lx %in% c("nontumor", "non_tumor", "non-tumor", "nontum", "non tumor"), "Non_Tumor", NA))
}
all_region_cells <- all_region_cells %>%
  mutate(
    Area = std_area(Region),
    Group = factor(ifelse(grepl("P0", Sample), "Control", "Case"), levels = c("Control","Case")),
    Pair_ID = sub("^SPE_(\\d+)[P][0-9A-Za-z]*no$", "\\1", Sample)
  ) %>%
  filter(!is.na(Area))

# save raw table
write.xlsx(all_region_cells, file = "cell_proportions_by_region.xlsx")

# ---- restrict to complete 1:1 pairs per Area ----
valid_pairs_per_area <- all_region_cells %>%
  distinct(Pair_ID, Group, Area) %>%
  count(Pair_ID, Area, Group) %>%
  pivot_wider(names_from = Group, values_from = n, values_fill = 0) %>%
  filter(Control == 1, Case == 1) %>%
  select(Pair_ID, Area)

region_cells_completepairs <- all_region_cells %>%
  inner_join(valid_pairs_per_area, by = c("Pair_ID","Area"))

# ---- immune cells only and fill missing types with 0 within each Pair x Group x Area ----
immune_cells <- region_cells_completepairs %>%
  filter(Cell_type != "Melanoma") %>%
  group_by(Pair_ID, Group, Area, Cell_type) %>%
  summarise(Proportion = sum(Proportion), .groups = "drop") %>%
  group_by(Pair_ID, Area, Cell_type) %>%
  complete(Group = factor(c("Control","Case"), levels = c("Control","Case")),
           fill = list(Proportion = 0)) %>%
  ungroup()

# ---- paired Wilcoxon per (Area, Cell_type); BH adjust across all tests ----
wide <- immune_cells %>%
  pivot_wider(names_from = Group, values_from = Proportion, values_fill = 0)

stats <- wide %>%
  group_by(Area, Cell_type) %>%
  summarise(
    p = {
      x <- Case; y <- Control
      if (length(x) == length(y) && length(x) > 0) wilcox.test(x, y, paired = TRUE)$p.value else NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    label = paste0("p=", signif(p, 3), ifelse(!is.na(p.adj) & p.adj < 0.05, "*", ""))
  )

# ---- plot across areas and cell types with p-values ----
area_plot <- ggplot(immune_cells, aes(x = Group, y = Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5) +
  geom_text(
    data = stats,
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE, vjust = 1.3, size = 3
  ) +
  scale_fill_manual(values = c("Case" = "#E69F00", "Control" = "#56B4E9")) +
  facet_grid(Area ~ Cell_type, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = NULL, y = "Cell Proportion") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "white"),
    panel.spacing = unit(1.5, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.8, 0.5, 0.5, 0.5), "cm"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# save
ggsave("immune_infiltration_by_area.png", plot = area_plot, width = 12, height = 8, dpi = 300)
ggsave("immune_infiltration_by_area.pdf", plot = area_plot, width = 12, height = 8)


# make summary table of raw + adjusted p-values per (Area, Cell_type)
# Summarise mean and median proportions per group, area, and cell type
prop_summary <- immune_cells %>%
  group_by(Area, Cell_type, Group) %>%
  summarise(
    mean_prop   = mean(Proportion, na.rm = TRUE),
    median_prop = median(Proportion, na.rm = TRUE),
    .groups = "drop"
  )

# Combine with stats (p-values)
summary_table <- stats %>%
  left_join(prop_summary, by = c("Area","Cell_type")) %>%
  arrange(Area, Cell_type, Group)

# Save to Excel
write.xlsx(summary_table, file = "immune_proportions_stats_with_values.xlsx")

# Preview
print(summary_table, n = 40)




