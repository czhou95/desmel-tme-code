# Set your working directory here
# setwd("your/working/directory")

# Publication-ready analysis: TCR Repertoire Diversity Metrics
# Case-Control Comparison - D-ESMEL study

# install.packages("readxl")  # uncomment to install if needed
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)  # for diversity calculations

# Load data
all_tcr_data <- read_excel("all_tcr_data.xlsx")
# Update the path below to point to your local copy of the key file
key_file <- read_excel("THREAT_sample_key.xlsx")

# Process key file to extract pair and case/control status
key_file <- key_file %>%
  mutate(
    pair_id = as.numeric(str_extract(`ID customer`, "^\\d+")),
    case_control = ifelse(grepl("-P-0$", `ID customer`), "Control", "Case")
  ) %>%
  rename(sample_id = `SkylineDx Sample ID`)

key_file %>% filter(is.na(pair_id)) %>% select(`ID customer`)

# Calculate TCR diversity metrics per sample
metrics_df <- all_tcr_data %>%
  group_by(sample_id) %>%
  summarize(
    richness = n_distinct(cloneId),
    shannon_entropy = diversity(readCount, index = "shannon"),
    evenness = diversity(readCount, index = "simpson"),
    .groups = "drop"
  )

# Sort by pair for correct pairing
merged_metrics <- merged_metrics %>% arrange(pair_id, case_control)

# Split by group
control_samples <- merged_metrics %>% filter(group == "Control")
case_samples <- merged_metrics %>% filter(group == "Case")

# Find pairs that have BOTH case and control
complete_pairs <- intersect(control_samples$pair_id, case_samples$pair_id)
cat("Complete matched pairs:", length(complete_pairs), "\n")

# Keep only complete pairs
control_samples <- control_samples %>% 
  filter(pair_id %in% complete_pairs) %>% 
  arrange(pair_id)
case_samples <- case_samples %>% 
  filter(pair_id %in% complete_pairs) %>% 
  arrange(pair_id)

# Verify pairing
stopifnot(nrow(control_samples) == nrow(case_samples))
stopifnot(all(control_samples$pair_id == case_samples$pair_id))

cat("Analysis includes", nrow(control_samples), "matched pairs\n")


# =============================================================================
# STATISTICAL ANALYSIS
# =============================================================================

# Define metrics to analyze
metrics <- c("richness", "shannon_entropy", "evenness")
metric_names <- c("Richness", "Shannon entropy", "Evenness")
names(metric_names) <- metrics

# Perform Wilcoxon signed-rank test for each metric
results_list <- lapply(metrics, function(metric) {
  control_values <- control_samples[[metric]]
  case_values <- case_samples[[metric]]
  
  test_result <- wilcox.test(control_values, case_values, paired = TRUE)
  
  data.frame(
    metric = metric,
    metric_name = metric_names[metric],
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
results_df$control_median_iqr <- sprintf("%.2f (%.2f–%.2f)", 
                                         results_df$control_median,
                                         results_df$control_iqr_25,
                                         results_df$control_iqr_75)
results_df$case_median_iqr <- sprintf("%.2f (%.2f–%.2f)",
                                      results_df$case_median,
                                      results_df$case_iqr_25,
                                      results_df$case_iqr_75)

# Print results table
cat("\n", paste(rep("=", 100), collapse = ""), "\n")
cat("TCR REPERTOIRE DIVERSITY METRICS\n")
cat(paste(rep("=", 100), collapse = ""), "\n\n")
cat(sprintf("%-20s %-25s %-25s %-12s %-12s\n", 
            "Metric", "Control median (IQR)", "Case median (IQR)", "p-value", "p (adj)"))
cat(paste(rep("-", 100), collapse = ""), "\n")

for (i in 1:nrow(results_df)) {
  p_raw <- ifelse(results_df$p_value[i] < 0.001, "<0.001", sprintf("%.3f", results_df$p_value[i]))
  p_adj <- ifelse(results_df$p_adjusted[i] < 0.001, "<0.001", sprintf("%.3f", results_df$p_adjusted[i]))
  cat(sprintf("%-20s %-25s %-25s %-12s %-12s\n",
              results_df$metric_name[i],
              results_df$control_median_iqr[i],
              results_df$case_median_iqr[i],
              p_raw, p_adj))
}

cat(paste(rep("-", 100), collapse = ""), "\n")
cat("p-value: Wilcoxon signed-rank test (paired); p (adj): Benjamini-Hochberg adjusted\n")
cat("n = ", nrow(control_samples), " matched pairs\n")
cat(paste(rep("=", 100), collapse = ""), "\n")

# =============================================================================
# PUBLICATION-READY FIGURE
# =============================================================================

# Transform data to long format
data_long <- merged_metrics %>%
  pivot_longer(
    cols = all_of(metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric_name = factor(metric_names[metric], levels = metric_names))

# Merge p-values for annotation
pval_df <- results_df %>%
  select(metric, p_adjusted) %>%
  mutate(
    metric_name = factor(metric_names[metric], levels = metric_names),
    sig_label = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01 ~ "**",
      p_adjusted < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Calculate bracket positions
bracket_positions <- data_long %>%
  group_by(metric_name) %>%
  summarise(
    y_max = max(value, na.rm = TRUE),
    y_min = min(value, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  mutate(
    bracket_y = y_max + 0.06 * y_range,
    label_y = y_max + 0.10 * y_range
  ) %>%
  left_join(pval_df %>% select(metric_name, sig_label), by = "metric_name")

# Create plot
p <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, 
               color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1, color = "black") +
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
            inherit.aes = FALSE, size = 4) +
  facet_wrap(~ metric_name, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Control" = "#4DAF4A", "Case" = "#E41A1C")) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    strip.text = element_text(size = 11, face = "plain"),
    strip.background = element_blank(),
    legend.position = "none",
    panel.spacing = unit(1, "lines")
  )

# Save plot
# Set your output directory here — plots will be saved to this location
# setwd("your/output/folder")

png("Figure_TCR_Diversity.png", width = 9, height = 4, units = "in", res = 300)
print(p)
dev.off()

ggsave("Figure_TCR_Diversity.pdf", p, width = 9, height = 4, units = "in")

cat("\nPlot saved as:\n")
cat("  - Figure_TCR_Diversity.png\n")
cat("  - Figure_TCR_Diversity.pdf\n")

# =============================================================================
# SAVE RESULTS TABLE
# =============================================================================

output_table <- results_df %>%
  select(metric_name, control_median_iqr, case_median_iqr, p_value, p_adjusted) %>%
  rename(
    `Metric` = metric_name,
    `Control median (IQR)` = control_median_iqr,
    `Case median (IQR)` = case_median_iqr,
    `p-value (raw)` = p_value,
    `p-value (BH adjusted)` = p_adjusted
  )

write.csv(output_table, "Table_TCR_Diversity.csv", row.names = FALSE)
cat("  - Table_TCR_Diversity.csv\n")


### BCR ###

# =============================================================================
# BCR Repertoire Analysis — Publication-Ready Code
# D-ESMEL Study
# =============================================================================
#
# This script performs B-cell receptor (BCR) repertoire analysis from MiXCR
# output, computing richness, Shannon entropy, and evenness (Gini-Simpson)
# for matched case-control pairs. It also tests whether the absence of
# detectable B-cell infiltration differs between groups (McNemar's test)
# and performs a sensitivity analysis including samples with zero richness.
# =============================================================================

library(tidyverse)
library(vegan)

# =============================================================================
# 1. READ AND PROCESS BCR DATA
# =============================================================================

# Update the path below to point to your local MiXCR BCR output directory
main_dir <- "path/to/mixcr_bcr_output"
sample_dirs <- list.dirs(main_dir, full.names = TRUE, recursive = FALSE)

read_bcr_file <- function(file_path) {
  read_tsv(file_path, col_types = cols(.default = "c"), show_col_types = FALSE)
}

process_sample_dir <- function(sample_dir) {
  bcr_files <- list.files(sample_dir, pattern = "\\.clones_IG[HKL]\\.tsv$", full.names = TRUE)
  if (length(bcr_files) == 0) return(NULL)
  
  sample_bcr_data <- bind_rows(lapply(bcr_files, read_bcr_file))
  sample_bcr_data$readCount <- as.numeric(sample_bcr_data$readCount)
  sample_bcr_data$cloneId <- as.numeric(sample_bcr_data$cloneId)
  sample_bcr_data$sample_id <- basename(sample_dir)
  
  return(sample_bcr_data)
}

all_bcr_data <- bind_rows(lapply(sample_dirs, process_sample_dir))

# setwd("your/working/directory")
write_tsv(all_bcr_data, "all_bcr_data.tsv")

# =============================================================================
# 2. CALCULATE REPERTOIRE METRICS
# =============================================================================

metrics_df <- all_bcr_data %>%
  group_by(sample_id) %>%
  summarize(
    richness = n_distinct(cloneId),
    shannon_entropy = diversity(readCount, index = "shannon"),
    evenness = diversity(readCount, index = "simpson"),
    .groups = "drop"
  ) %>%
  as.data.frame()

# Merge with metadata
merged_metrics <- inner_join(
  metrics_df, key_EMC_SkylineDx,
  by = c("sample_id" = "SampleID_Skyline")
)

# Create group label for plotting
merged_metrics$group <- factor(
  ifelse(merged_metrics$case_control == 0, "Control", "Case"),
  levels = c("Control", "Case")
)

# =============================================================================
# 3. DETECTABLE B-CELL INFILTRATION ANALYSIS (McNemar's test)
# =============================================================================

samples_with_bcr <- unique(all_bcr_data$sample_id)

key_EMC_SkylineDx <- key_EMC_SkylineDx %>%
  mutate(bcr_detected = ifelse(SampleID_Skyline %in% samples_with_bcr, 1, 0))

paired_bcr <- key_EMC_SkylineDx %>%
  select(set, case_control, bcr_detected) %>%
  pivot_wider(
    names_from = case_control,
    values_from = bcr_detected,
    names_prefix = "bcr_"
  )

mcnemar_table <- table(
  Control = paired_bcr$bcr_0,
  Case = paired_bcr$bcr_1
)
mcnemar_result <- mcnemar.test(mcnemar_table)

# Exact binomial test for small discordant cell counts
n_discordant_10 <- mcnemar_table[2, 1]
n_discordant_01 <- mcnemar_table[1, 2]
n_discordant_total <- n_discordant_10 + n_discordant_01

if (n_discordant_total < 25) {
  exact_result <- binom.test(n_discordant_10, n_discordant_total, p = 0.5)
}

# =============================================================================
# 4. IDENTIFY COMPLETE PAIRS
# =============================================================================

richness_wide_all <- merged_metrics %>%
  select(set, case_control, richness) %>%
  pivot_wider(
    names_from = case_control,
    values_from = richness,
    names_prefix = "richness_"
  )

complete_pairs <- richness_wide_all %>%
  filter(!is.na(richness_1) & !is.na(richness_0))

complete_pair_sets <- complete_pairs$set
n_complete_pairs <- nrow(complete_pairs)
cat("Complete BCR pairs:", n_complete_pairs, "\n")

# Filter to complete pairs
bcr_complete <- merged_metrics %>%
  filter(set %in% complete_pair_sets)

control_samples <- bcr_complete %>% filter(case_control == 0) %>% arrange(set)
case_samples <- bcr_complete %>% filter(case_control == 1) %>% arrange(set)

stopifnot(nrow(control_samples) == nrow(case_samples))
stopifnot(all(control_samples$set == case_samples$set))

# =============================================================================
# 5. STATISTICAL ANALYSIS
# =============================================================================

metrics <- c("richness", "shannon_entropy", "evenness")
metric_names <- c("Richness", "Shannon entropy", "Evenness")
names(metric_names) <- metrics

results_list <- lapply(metrics, function(metric) {
  control_values <- control_samples[[metric]]
  case_values <- case_samples[[metric]]
  
  test_result <- wilcox.test(control_values, case_values, paired = TRUE)
  
  data.frame(
    metric = metric,
    metric_name = metric_names[metric],
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

# Format for manuscript
results_df$control_median_iqr <- sprintf("%.2f (%.2f\u2013%.2f)",
                                         results_df$control_median,
                                         results_df$control_iqr_25,
                                         results_df$control_iqr_75)
results_df$case_median_iqr <- sprintf("%.2f (%.2f\u2013%.2f)",
                                      results_df$case_median,
                                      results_df$case_iqr_25,
                                      results_df$case_iqr_75)

# Print results
cat("\n", paste(rep("=", 100), collapse = ""), "\n")
cat("BCR REPERTOIRE DIVERSITY METRICS (n =", n_complete_pairs, "complete pairs)\n")
cat(paste(rep("=", 100), collapse = ""), "\n\n")
cat(sprintf("%-20s %-25s %-25s %-12s %-12s\n",
            "Metric", "Control median (IQR)", "Case median (IQR)", "p-value", "p (adj)"))
cat(paste(rep("-", 100), collapse = ""), "\n")

for (i in 1:nrow(results_df)) {
  p_raw <- ifelse(results_df$p_value[i] < 0.001, "<0.001", sprintf("%.3f", results_df$p_value[i]))
  p_adj <- ifelse(results_df$p_adjusted[i] < 0.001, "<0.001", sprintf("%.3f", results_df$p_adjusted[i]))
  cat(sprintf("%-20s %-25s %-25s %-12s %-12s\n",
              results_df$metric_name[i],
              results_df$control_median_iqr[i],
              results_df$case_median_iqr[i],
              p_raw, p_adj))
}
cat(paste(rep("-", 100), collapse = ""), "\n")
cat("p-value: Wilcoxon signed-rank test (paired); p (adj): Benjamini-Hochberg adjusted\n")

# =============================================================================
# 6. SENSITIVITY ANALYSIS: RICHNESS WITH ZEROS (all 178 pairs)
# =============================================================================

bcr_richness_all <- key_EMC_SkylineDx %>%
  select(SampleID_Skyline, case_control, set) %>%
  left_join(
    metrics_df %>% select(sample_id, richness),
    by = c("SampleID_Skyline" = "sample_id")
  ) %>%
  mutate(richness_with_zeros = replace_na(richness, 0))

richness_wide_zeros <- bcr_richness_all %>%
  select(set, case_control, richness_with_zeros) %>%
  pivot_wider(
    names_from = case_control,
    values_from = richness_with_zeros,
    names_prefix = "richness_"
  )

wilcox_richness_zeros <- wilcox.test(
  richness_wide_zeros$richness_1,
  richness_wide_zeros$richness_0,
  paired = TRUE
)

cat("\n", paste(rep("=", 100), collapse = ""), "\n")
cat("McNEMAR TEST \u2014 Detectable B-cell infiltration\n")
cat(paste(rep("=", 100), collapse = ""), "\n")
cat("McNemar p-value:", round(mcnemar_result$p.value, 2), "\n")
if (exists("exact_result")) {
  cat("Exact binomial p-value:", round(exact_result$p.value, 2), "\n")
  cat("Discordant pairs:", n_discordant_total, "\n")
}
cat("\nSENSITIVITY \u2014 Richness with zeros (n = 178 pairs)\n")
cat("Paired Wilcoxon p-value:", round(wilcox_richness_zeros$p.value, 2), "\n")
cat(paste(rep("=", 100), collapse = ""), "\n")

# =============================================================================
# 7. PUBLICATION-READY FIGURE
# =============================================================================

data_long <- bcr_complete %>%
  pivot_longer(
    cols = all_of(metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric_name = factor(metric_names[metric], levels = metric_names))

pval_df <- results_df %>%
  select(metric, p_adjusted) %>%
  mutate(
    metric_name = factor(metric_names[metric], levels = metric_names),
    sig_label = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

bracket_positions <- data_long %>%
  group_by(metric_name) %>%
  summarise(
    y_max = max(value, na.rm = TRUE),
    y_min = min(value, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  mutate(
    bracket_y = y_max + 0.06 * y_range,
    label_y = y_max + 0.10 * y_range
  ) %>%
  left_join(pval_df %>% select(metric_name, sig_label), by = "metric_name")

p <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7,
               color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1, color = "black") +
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
            inherit.aes = FALSE, size = 4) +
  facet_wrap(~ metric_name, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Control" = "#4DAF4A", "Case" = "#E41A1C")) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    strip.text = element_text(size = 11, face = "plain"),
    strip.background = element_blank(),
    legend.position = "none",
    panel.spacing = unit(1, "lines")
  )

# Save figure
# setwd("your/working/directory")
png("Figure_BCR_Diversity.png", width = 9, height = 4, units = "in", res = 300)
print(p)
dev.off()

ggsave("Figure_BCR_Diversity.pdf", p, width = 9, height = 4, units = "in")

# =============================================================================
# 8. SAVE TABLE
# =============================================================================

output_table <- results_df %>%
  select(metric_name, control_median_iqr, case_median_iqr, p_value, p_adjusted) %>%
  rename(
    `Metric` = metric_name,
    `Control median (IQR)` = control_median_iqr,
    `Case median (IQR)` = case_median_iqr,
    `p-value (raw)` = p_value,
    `p-value (BH adjusted)` = p_adjusted
  )

write.csv(output_table, "Table_BCR_Diversity.csv", row.names = FALSE)

cat("\nFiles saved:\n")
cat("  - Figure_BCR_Diversity.png\n")
cat("  - Figure_BCR_Diversity.pdf\n")
cat("  - Table_BCR_Diversity.csv\n")

