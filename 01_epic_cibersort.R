# Publication-ready analysis: Total Immune Cells by Case/Control Status
# EPIC deconvolution data - D-ESMEL study

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
# Load EPIC data
# Update the path below to point to your local copy of the input file
data <- read_excel("EPIC_output_withprognosis.xlsx")

# Compute total immune cell percentage (including CAFs and Endothelial cells)
data$total_immune_cells_2 <- rowSums(data[, c("Bcells", "CD4_Tcells", "CD8_Tcells", 
                                              "Macrophages", "NKcells", "CAFs", "Endothelial")])

# Recode prognosis to Case/Control
data$group <- factor(ifelse(data$prognosis == "poor", "Case", "Control"),
                     levels = c("Control", "Case"))

# Wilcoxon signed-rank test (paired)
# Arrange by pair to ensure correct pairing
data_ordered <- data %>% arrange(pair, prognosis)
control_samples <- data_ordered %>% filter(group == "Control") %>% pull(total_immune_cells_2)
case_samples <- data_ordered %>% filter(group == "Case") %>% pull(total_immune_cells_2)

wilcox_result <- wilcox.test(control_samples, case_samples, paired = TRUE)
p_val <- wilcox_result$p.value

# Summary statistics for table
summary_stats <- data %>%
  group_by(group) %>%
  summarise(
    n = n(),
    Median = median(total_immune_cells_2, na.rm = TRUE),
    IQR_25 = quantile(total_immune_cells_2, 0.25, na.rm = TRUE),
    IQR_75 = quantile(total_immune_cells_2, 0.75, na.rm = TRUE),
    Mean = mean(total_immune_cells_2, na.rm = TRUE),
    SD = sd(total_immune_cells_2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Median (IQR)` = sprintf("%.3f (%.3f–%.3f)", Median, IQR_25, IQR_75),
    `Mean ± SD` = sprintf("%.3f ± %.3f", Mean, SD)
  )

# Print results table
cat("\n=== RESULTS TABLE FOR MANUSCRIPT ===\n\n")
cat("Table: Total immune cell fraction by case-control status (EPIC deconvolution)\n\n")
print(summary_stats %>% select(group, n, `Median (IQR)`, `Mean ± SD`))
cat("\nWilcoxon signed-rank test (paired): V =", wilcox_result$statistic, 
    ", p =", format(p_val, digits = 3), "\n")

# Publication-ready plot
p <- ggplot(data, aes(x = group, y = total_immune_cells_2, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.2, color = "black") +
  scale_fill_manual(values = c("Control" = "#4DAF4A", "Case" = "#E41A1C")) +
  labs(
    x = NULL,
    y = "Total immune cell fraction EPIC"
  ) +
  annotate("segment", x = 1, xend = 2, y = max(data$total_immune_cells_2) + 0.02, 
           yend = max(data$total_immune_cells_2) + 0.02, linewidth = 0.4) +
  annotate("text", x = 1.5, y = max(data$total_immune_cells_2) + 0.035, 
           label = ifelse(p_val < 0.001, "p < 0.001", 
                          paste0("p = ", format(round(p_val, 3), nsmall = 3))),
           size = 3.5, family = "sans") +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.title.y = element_text(size = 11, margin = margin(r = 10)),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.position = "none",
    plot.margin = margin(10, 15, 10, 10)
  ) +
  coord_cartesian(ylim = c(min(data$total_immune_cells_2) - 0.02, 
                           max(data$total_immune_cells_2) + 0.06))

# Save plot
# Set your working directory here — plots will be saved to this location
# setwd("your/output/folder")
ggsave("Figure_EPIC_TotalImmuneCells.pdf", p, width = 3.5, height = 4, units = "in", dpi = 300)
png("Figure_EPIC_TotalImmuneCells.png", width = 3.5, height = 4, units = "in", res = 300)
print(p)
dev.off()

cat("\nPlot saved as:\n  - Figure_EPIC_TotalImmuneCells.pdf\n  - Figure_EPIC_TotalImmuneCells.png\n")

## Individual cell types
# Split and order by pair to ensure correct pairing
good_samples <- data[data$prognosis == "good", ]
poor_samples <- data[data$prognosis == "poor", ]
good_samples <- good_samples[order(good_samples$pair), ]
poor_samples <- poor_samples[order(poor_samples$pair), ]

# Define immune cell columns explicitly
immune_columns <- c("Bcells", "CAFs", "CD4_Tcells", "CD8_Tcells", 
                    "Endothelial", "Macrophages", "NKcells", "otherCells")

# Clean names for display
cell_names <- c("B cells", "CAFs", "CD4+ T cells", "CD8+ T cells",
                "Endothelial cells", "Macrophages", "NK cells", "Other cells")
names(cell_names) <- immune_columns

# Perform Wilcoxon signed-rank test for each immune cell type
results_list <- lapply(immune_columns, function(cell) {
  good_values <- good_samples[[cell]]
  poor_values <- poor_samples[[cell]]
  
  test_result <- wilcox.test(good_values, poor_values, paired = TRUE)
  
  data.frame(
    cell_type = cell,
    cell_name = cell_names[cell],
    control_median = median(good_values, na.rm = TRUE),
    control_iqr_25 = quantile(good_values, 0.25, na.rm = TRUE),
    control_iqr_75 = quantile(good_values, 0.75, na.rm = TRUE),
    case_median = median(poor_values, na.rm = TRUE),
    case_iqr_25 = quantile(poor_values, 0.25, na.rm = TRUE),
    case_iqr_75 = quantile(poor_values, 0.75, na.rm = TRUE),
    W_statistic = test_result$statistic,
    p_value = test_result$p.value
  )
})

results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

# Benjamini-Hochberg correction
results_df$p_adjusted <- p.adjust(results_df$p_value, method = "BH")

# Format for manuscript table
results_df$control_median_iqr <- sprintf("%.3f (%.3f–%.3f)", 
                                         results_df$control_median,
                                         results_df$control_iqr_25,
                                         results_df$control_iqr_75)
results_df$case_median_iqr <- sprintf("%.3f (%.3f–%.3f)",
                                      results_df$case_median,
                                      results_df$case_iqr_25,
                                      results_df$case_iqr_75)

# Print results table
cat("\n", paste(rep("=", 90), collapse = ""), "\n")
cat("RESULTS TABLE FOR MANUSCRIPT\n")
cat(paste(rep("=", 90), collapse = ""), "\n\n")
cat("Table: Immune cell fractions by case-control status (EPIC deconvolution)\n")
cat(paste(rep("-", 90), collapse = ""), "\n")
cat(sprintf("%-18s %-28s %-28s %-10s %-10s\n", 
            "Cell type", "Control median (IQR)", "Case median (IQR)", "p-value", "p (adj)"))
cat(paste(rep("-", 90), collapse = ""), "\n")

for (i in 1:nrow(results_df)) {
  p_raw <- ifelse(results_df$p_value[i] < 0.001, "<0.001", sprintf("%.3f", results_df$p_value[i]))
  p_adj <- ifelse(results_df$p_adjusted[i] < 0.001, "<0.001", sprintf("%.3f", results_df$p_adjusted[i]))
  cat(sprintf("%-18s %-28s %-28s %-10s %-10s\n",
              results_df$cell_name[i],
              results_df$control_median_iqr[i],
              results_df$case_median_iqr[i],
              p_raw, p_adj))
}

cat(paste(rep("-", 90), collapse = ""), "\n")
cat("p-value: Wilcoxon signed-rank test (paired); p (adj): Benjamini-Hochberg adjusted\n")
cat(paste(rep("=", 90), collapse = ""), "\n")

# Transform data to long format for plotting
data_long <- data %>%
  pivot_longer(
    cols = all_of(immune_columns),
    names_to = "cell_type",
    values_to = "value"
  ) %>%
  mutate(cell_name = factor(cell_names[cell_type], levels = cell_names))

# Merge adjusted p-values for annotation
pval_df <- results_df %>%
  select(cell_type, p_adjusted) %>%
  mutate(
    cell_name = factor(cell_names[cell_type], levels = cell_names),
    sig_label = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01 ~ "**",
      p_adjusted < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Calculate y positions for significance brackets
bracket_positions <- data_long %>%
  group_by(cell_name) %>%
  summarise(
    y_max = max(value, na.rm = TRUE),
    y_min = min(value, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  mutate(
    bracket_y = y_max + 0.05 * y_range,
    label_y = y_max + 0.08 * y_range
  ) %>%
  left_join(pval_df %>% select(cell_name, sig_label), by = "cell_name")

# Publication-ready plot
p <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, 
               color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8, color = "black") +
  geom_segment(data = bracket_positions,
               aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = bracket_positions,
               aes(x = 1, xend = 1, y = bracket_y - 0.02 * y_range, yend = bracket_y),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = bracket_positions,
               aes(x = 2, xend = 2, y = bracket_y - 0.02 * y_range, yend = bracket_y),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = bracket_positions,
            aes(x = 1.5, y = label_y, label = sig_label),
            inherit.aes = FALSE, size = 3) +
  facet_wrap(~ cell_name, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("Control" = "#4DAF4A", "Case" = "#E41A1C")) +
  labs(x = NULL, y = "Cell fraction") +
  theme_classic(base_size = 10) +
  theme(
    axis.title.y = element_text(size = 11, margin = margin(r = 10)),
    axis.text = element_text(size = 9, color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    strip.text = element_text(size = 10, face = "plain"),
    strip.background = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.8, "lines")
  )

# Save plot
png("Figure_EPIC_ImmuneCells_Individual.png", width = 10, height = 5.5, units = "in", res = 300)
print(p)
dev.off()

ggsave("Figure_EPIC_ImmuneCells_Individual.pdf", p, width = 10, height = 5.5, units = "in")

cat("\nPlot saved as:\n")
cat("  - Figure_EPIC_ImmuneCells_Individual.png\n")
cat("  - Figure_EPIC_ImmuneCells_Individual.pdf\n")

# Save results table as CSV
output_table <- results_df %>%
  select(cell_name, control_median_iqr, case_median_iqr, p_value, p_adjusted) %>%
  rename(
    `Cell type` = cell_name,
    `Control median (IQR)` = control_median_iqr,
    `Case median (IQR)` = case_median_iqr,
    `p-value (raw)` = p_value,
    `p-value (BH adjusted)` = p_adjusted
  )

write.csv(output_table, "Table_EPIC_ImmuneCells_Individual.csv", row.names = FALSE)
cat("  - Table_EPIC_ImmuneCells_Individual.csv\n")

###CIBERSORT
# Load data
# Load CIBERSORT data
# Update the path below to point to your local copy of the input file
data <- read_excel("CIBERSORT_output_withprognosis.xlsx")

# Define immune cell columns (22 cell types from LM22 signature)
immune_columns <- c("Bcellsnaive", "Bcellsmemory", "Plasmacells", 
                    "T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting",
                    "T cells CD4 memory activated", "T cells follicular helper",
                    "T cells regulatory (Tregs)", "T cells gamma delta",
                    "NK cells resting", "NK cells activated",
                    "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
                    "Dendritic cells resting", "Dendritic cells activated",
                    "Mast cells resting", "Mast cells activated",
                    "Eosinophils", "Neutrophils")

# Clean names for display
cell_names <- c("B cells naive", "B cells memory", "Plasma cells",
                "CD8+ T cells", "CD4+ T cells naive", "CD4+ T cells memory resting",
                "CD4+ T cells memory activated", "T follicular helper",
                "Tregs", "γδ T cells",
                "NK cells resting", "NK cells activated",
                "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
                "Dendritic cells resting", "Dendritic cells activated",
                "Mast cells resting", "Mast cells activated",
                "Eosinophils", "Neutrophils")
names(cell_names) <- immune_columns

# Recode prognosis to Case/Control
data$group <- factor(ifelse(data$prognosis == "poor", "Case", "Control"),
                     levels = c("Control", "Case"))

# =============================================================================
# PART 1: TOTAL IMMUNE CELLS
# =============================================================================

# Compute total immune cell fraction
data$total_immune_cells <- rowSums(data[, immune_columns], na.rm = TRUE)

# Split and order by pair for correct pairing
control_samples <- data[data$group == "Control", ]
case_samples <- data[data$group == "Case", ]
control_samples <- control_samples[order(control_samples$pair), ]
case_samples <- case_samples[order(case_samples$pair), ]

# Wilcoxon signed-rank test (paired) for total immune cells
wilcox_total <- wilcox.test(control_samples$total_immune_cells, 
                            case_samples$total_immune_cells, 
                            paired = TRUE)

# Summary statistics for total immune cells
total_summary <- data %>%
  group_by(group) %>%
  summarise(
    n = n(),
    Median = median(total_immune_cells, na.rm = TRUE),
    IQR_25 = quantile(total_immune_cells, 0.25, na.rm = TRUE),
    IQR_75 = quantile(total_immune_cells, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("CIBERSORT: TOTAL IMMUNE CELL FRACTION\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
cat(sprintf("%-10s %-6s %-30s\n", "Group", "n", "Median (IQR)"))
cat(paste(rep("-", 50), collapse = ""), "\n")
for (i in 1:nrow(total_summary)) {
  cat(sprintf("%-10s %-6d %.3f (%.3f–%.3f)\n",
              total_summary$group[i], total_summary$n[i],
              total_summary$Median[i], total_summary$IQR_25[i], total_summary$IQR_75[i]))
}
cat(paste(rep("-", 50), collapse = ""), "\n")
cat(sprintf("Wilcoxon signed-rank test: V = %.0f, p = %.3f\n", 
            wilcox_total$statistic, wilcox_total$p.value))
cat(paste(rep("=", 70), collapse = ""), "\n")

# Plot for total immune cells
p_total <- ggplot(data, aes(x = group, y = total_immune_cells, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, 
               color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.2, color = "black") +
  scale_fill_manual(values = c("Control" = "#4DAF4A", "Case" = "#E41A1C")) +
  labs(x = NULL, y = "Total immune cell fraction") +
  annotate("segment", x = 1, xend = 2, 
           y = max(data$total_immune_cells) + 0.02, 
           yend = max(data$total_immune_cells) + 0.02, linewidth = 0.4) +
  annotate("text", x = 1.5, y = max(data$total_immune_cells) + 0.04,
           label = sprintf("p = %.3f", wilcox_total$p.value), size = 3.5) +
  theme_classic(base_size = 11) +
  theme(
    axis.title.y = element_text(size = 11, margin = margin(r = 10)),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(min(data$total_immune_cells) - 0.02, 
                           max(data$total_immune_cells) + 0.08))

# Save total immune cells plot
png("Figure_CIBERSORT_TotalImmuneCells.png", width = 3.5, height = 4, units = "in", res = 300)
print(p_total)
dev.off()

ggsave("Figure_CIBERSORT_TotalImmuneCells.pdf", p_total, width = 3.5, height = 4, units = "in")

# =============================================================================
# PART 2: INDIVIDUAL IMMUNE CELL TYPES
# =============================================================================

# Perform Wilcoxon signed-rank test for each immune cell type
results_list <- lapply(immune_columns, function(cell) {
  control_values <- control_samples[[cell]]
  case_values <- case_samples[[cell]]
  
  test_result <- wilcox.test(control_values, case_values, paired = TRUE)
  
  data.frame(
    cell_type = cell,
    cell_name = cell_names[cell],
    control_median = median(control_values, na.rm = TRUE),
    control_iqr_25 = quantile(control_values, 0.25, na.rm = TRUE),
    control_iqr_75 = quantile(control_values, 0.75, na.rm = TRUE),
    case_median = median(case_values, na.rm = TRUE),
    case_iqr_25 = quantile(case_values, 0.25, na.rm = TRUE),
    case_iqr_75 = quantile(case_values, 0.75, na.rm = TRUE),
    W_statistic = test_result$statistic,
    p_value = test_result$p.value
  )
})

results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

# Benjamini-Hochberg correction
results_df$p_adjusted <- p.adjust(results_df$p_value, method = "BH")

# Format for manuscript table
results_df$control_median_iqr <- sprintf("%.3f (%.3f–%.3f)", 
                                         results_df$control_median,
                                         results_df$control_iqr_25,
                                         results_df$control_iqr_75)
results_df$case_median_iqr <- sprintf("%.3f (%.3f–%.3f)",
                                      results_df$case_median,
                                      results_df$case_iqr_25,
                                      results_df$case_iqr_75)

# Print results table
cat("\n", paste(rep("=", 100), collapse = ""), "\n")
cat("CIBERSORT: INDIVIDUAL IMMUNE CELL FRACTIONS\n")
cat(paste(rep("=", 100), collapse = ""), "\n\n")
cat(sprintf("%-30s %-25s %-25s %-10s %-10s\n", 
            "Cell type", "Control median (IQR)", "Case median (IQR)", "p-value", "p (adj)"))
cat(paste(rep("-", 100), collapse = ""), "\n")

for (i in 1:nrow(results_df)) {
  p_raw <- ifelse(results_df$p_value[i] < 0.001, "<0.001", sprintf("%.3f", results_df$p_value[i]))
  p_adj <- ifelse(results_df$p_adjusted[i] < 0.001, "<0.001", sprintf("%.3f", results_df$p_adjusted[i]))
  cat(sprintf("%-30s %-25s %-25s %-10s %-10s\n",
              results_df$cell_name[i],
              results_df$control_median_iqr[i],
              results_df$case_median_iqr[i],
              p_raw, p_adj))
}

cat(paste(rep("-", 100), collapse = ""), "\n")
cat("p-value: Wilcoxon signed-rank test (paired); p (adj): Benjamini-Hochberg adjusted\n")
cat(paste(rep("=", 100), collapse = ""), "\n")

# Transform data to long format for plotting
data_long <- data %>%
  pivot_longer(
    cols = all_of(immune_columns),
    names_to = "cell_type",
    values_to = "value"
  ) %>%
  mutate(cell_name = factor(cell_names[cell_type], levels = cell_names))

# Merge adjusted p-values for annotation
pval_df <- results_df %>%
  select(cell_type, p_adjusted) %>%
  mutate(
    cell_name = factor(cell_names[cell_type], levels = cell_names),
    sig_label = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01 ~ "**",
      p_adjusted < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Calculate y positions for significance brackets
bracket_positions <- data_long %>%
  group_by(cell_name) %>%
  summarise(
    y_max = max(value, na.rm = TRUE),
    y_min = min(value, na.rm = TRUE),
    y_range = max(y_max - y_min, 0.01),
    .groups = "drop"
  ) %>%
  mutate(
    bracket_y = y_max + 0.08 * y_range,
    label_y = y_max + 0.12 * y_range
  ) %>%
  left_join(pval_df %>% select(cell_name, sig_label), by = "cell_name")

# Publication-ready plot for individual cells (6 rows x 4 columns = 24 panels, using 22)
p_individual <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, 
               color = "black", linewidth = 0.3) +
  geom_jitter(width = 0.12, alpha = 0.25, size = 0.5, color = "black") +
  geom_segment(data = bracket_positions,
               aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_segment(data = bracket_positions,
               aes(x = 1, xend = 1, y = bracket_y - 0.03 * y_range, yend = bracket_y),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_segment(data = bracket_positions,
               aes(x = 2, xend = 2, y = bracket_y - 0.03 * y_range, yend = bracket_y),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = bracket_positions,
            aes(x = 1.5, y = label_y, label = sig_label),
            inherit.aes = FALSE, size = 2.5) +
  facet_wrap(~ cell_name, scales = "free_y", nrow = 4) +
  scale_fill_manual(values = c("Control" = "#4DAF4A", "Case" = "#E41A1C")) +
  labs(x = NULL, y = "Cell fraction") +
  theme_classic(base_size = 9) +
  theme(
    axis.title.y = element_text(size = 10, margin = margin(r = 8)),
    axis.text = element_text(size = 7, color = "black"),
    axis.text.x = element_text(size = 7),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    strip.text = element_text(size = 8, face = "plain"),
    strip.background = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines")
  )

# Save individual cells plot
png("Figure_CIBERSORT_ImmuneCells_Individual.png", width = 10, height = 10, units = "in", res = 300)
print(p_individual)
dev.off()

ggsave("Figure_CIBERSORT_ImmuneCells_Individual.pdf", p_individual, width = 10, height = 10, units = "in")

cat("\nPlots saved as:\n")
cat("  - Figure_CIBERSORT_TotalImmuneCells.png/.pdf\n")
cat("  - Figure_CIBERSORT_ImmuneCells_Individual.png/.pdf\n")

# =============================================================================
# SAVE RESULTS TABLES
# =============================================================================

# Total immune cells table
total_table <- data.frame(
  Group = c("Control", "Case", "Wilcoxon signed-rank test"),
  n = c(total_summary$n[1], total_summary$n[2], ""),
  `Median (IQR)` = c(
    sprintf("%.3f (%.3f–%.3f)", total_summary$Median[1], total_summary$IQR_25[1], total_summary$IQR_75[1]),
    sprintf("%.3f (%.3f–%.3f)", total_summary$Median[2], total_summary$IQR_25[2], total_summary$IQR_75[2]),
    sprintf("V = %.0f, p = %.3f", wilcox_total$statistic, wilcox_total$p.value)
  ),
  check.names = FALSE
)
write.csv(total_table, "Table_CIBERSORT_TotalImmuneCells.csv", row.names = FALSE)

# Individual cells table
output_table <- results_df %>%
  select(cell_name, control_median_iqr, case_median_iqr, p_value, p_adjusted) %>%
  rename(
    `Cell type` = cell_name,
    `Control median (IQR)` = control_median_iqr,
    `Case median (IQR)` = case_median_iqr,
    `p-value (raw)` = p_value,
    `p-value (BH adjusted)` = p_adjusted
  )
write.csv(output_table, "Table_CIBERSORT_ImmuneCells_Individual.csv", row.names = FALSE)

cat("  - Table_CIBERSORT_TotalImmuneCells.csv\n")
cat("  - Table_CIBERSORT_ImmuneCells_Individual.csv\n")

