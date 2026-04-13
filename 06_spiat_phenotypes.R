################################################################################
# SPATIAL IMMUNE PHENOTYPE CLASSIFICATION - COMPLETE SCRIPT
# D-ESMEL TIME Study
# 
# Classifies tumors into: Infiltrated, Excluded, or Desert
# Based on CD8+ T-cell density and proximity to melanoma cells
#
# PREREQUISITE: Run density and MNND analyses first to have:
#   - complete_pairs_df (from density analysis)
#   - mnnd_results (from MNND analysis)
################################################################################

library(tidyverse)

# =============================================================================
# PART 1: PREPARE DATA FOR CLASSIFICATION
# =============================================================================

# Get CD8 density per sample
cd8_density <- complete_pairs_df %>%
  filter(Cell.Type == "Cytotoxic T-cell") %>%
  select(Sample, Patient, Group, Density = Density_per_mm2)

# Get MNND (Melanoma to CD8) per sample
cd8_mnnd <- mnnd_results %>%
  select(Sample, Pair_ID, Group, MNND = Melanoma_Cytotoxic_T)

# Merge density and MNND
phenotype_data <- cd8_density %>%
  left_join(cd8_mnnd %>% select(Sample, MNND), by = "Sample") %>%
  filter(!is.na(Density) & !is.na(MNND))

cat("Samples with both density and MNND data:", nrow(phenotype_data), "\n")


# =============================================================================
# PART 2: CLASSIFY INTO IMMUNE PHENOTYPES
# =============================================================================

# Use median cutoffs (you could also use tertiles or literature-based cutoffs)
density_median <- median(phenotype_data$Density, na.rm = TRUE)
mnnd_median <- median(phenotype_data$MNND, na.rm = TRUE)

cat("\nCutoffs (median-based):\n")
cat(sprintf("  CD8 density: %.1f cells/mm²\n", density_median))
cat(sprintf("  MNND (Mel-CD8): %.1f µm\n", mnnd_median))

# Classification:
# - Infiltrated: High density + Low MNND (CD8 present AND close)
# - Excluded: High density + High MNND (CD8 present BUT far)
# - Desert: Low density (regardless of MNND)

phenotype_data <- phenotype_data %>%
  mutate(
    Density_cat = ifelse(Density >= density_median, "High", "Low"),
    MNND_cat = ifelse(MNND >= mnnd_median, "High", "Low"),
    Phenotype = case_when(
      Density_cat == "Low" ~ "Desert",
      Density_cat == "High" & MNND_cat == "Low" ~ "Infiltrated",
      Density_cat == "High" & MNND_cat == "High" ~ "Excluded"
    ),
    Phenotype = factor(Phenotype, levels = c("Infiltrated", "Excluded", "Desert"))
  )


# =============================================================================
# PART 3: STATISTICAL ANALYSIS
# =============================================================================

# Contingency table
contingency <- table(phenotype_data$Phenotype, phenotype_data$Group)
cat("\n========== CONTINGENCY TABLE ==========\n")
print(contingency)

# Proportions per group
cat("\n========== PROPORTIONS BY GROUP ==========\n")
prop_table <- phenotype_data %>%
  count(Group, Phenotype) %>%
  group_by(Group) %>%
  mutate(
    Total = sum(n),
    Percentage = n / Total * 100
  ) %>%
  ungroup()
print(prop_table)

# Chi-square test (or Fisher's exact if expected counts < 5)
# NOTE: These treat samples as independent - see paired analysis below
chi_test <- chisq.test(contingency)
fisher_test <- fisher.test(contingency)

cat("\n========== UNPAIRED STATISTICAL TESTS (for reference) ==========\n")
cat(sprintf("Chi-square test: X² = %.2f, df = %d, p = %.4f\n", 
            chi_test$statistic, chi_test$parameter, chi_test$p.value))
cat(sprintf("Fisher's exact test: p = %.4f\n", fisher_test$p.value))

# Check expected counts
cat("\nExpected counts (Chi-square):\n")
print(round(chi_test$expected, 1))


# =============================================================================
# PAIRED ANALYSIS - Stuart-Maxwell Test (appropriate for matched design)
# =============================================================================

cat("\n========== PAIRED ANALYSIS (matched case-control) ==========\n")

# Create paired contingency table: rows = Control phenotype, cols = Case phenotype
paired_data <- phenotype_data %>%
  select(Patient, Group, Phenotype) %>%
  pivot_wider(names_from = Group, values_from = Phenotype) %>%
  filter(!is.na(Control) & !is.na(Case))  # Keep only complete pairs

cat(sprintf("\nComplete matched pairs: %d\n", nrow(paired_data)))

# Cross-tabulation of paired phenotypes
paired_table <- table(Control = paired_data$Control, Case = paired_data$Case)
cat("\nPaired contingency table (Control phenotype vs Case phenotype):\n")
print(paired_table)

# Calculate discordant pairs
cat("\nDiscordant pairs:\n")
cat(sprintf("  Control=Infiltrated, Case=Desert: %d\n", 
            sum(paired_data$Control == "Infiltrated" & paired_data$Case == "Desert")))
cat(sprintf("  Control=Desert, Case=Infiltrated: %d\n", 
            sum(paired_data$Control == "Desert" & paired_data$Case == "Infiltrated")))
cat(sprintf("  Control=Infiltrated, Case=Excluded: %d\n", 
            sum(paired_data$Control == "Infiltrated" & paired_data$Case == "Excluded")))
cat(sprintf("  Control=Excluded, Case=Infiltrated: %d\n", 
            sum(paired_data$Control == "Excluded" & paired_data$Case == "Infiltrated")))

# Stuart-Maxwell test (generalization of McNemar for k×k tables)
# Using the coin package or manual calculation
# install.packages("coin")  # uncomment if not yet installed
library(coin)

# Alternative: McNemar-Bowker test of symmetry
mcnemar_bowker <- function(x) {
  # Test of symmetry for k×k table
  k <- nrow(x)
  if (k != ncol(x)) stop("Table must be square")
  
  # Calculate test statistic
  stat <- 0
  df <- 0
  for (i in 1:(k-1)) {
    for (j in (i+1):k) {
      if ((x[i,j] + x[j,i]) > 0) {
        stat <- stat + (x[i,j] - x[j,i])^2 / (x[i,j] + x[j,i])
        df <- df + 1
      }
    }
  }
  p_value <- pchisq(stat, df, lower.tail = FALSE)
  
  list(statistic = stat, df = df, p.value = p_value)
}

# Run McNemar-Bowker test
mb_test <- mcnemar_bowker(paired_table)
cat(sprintf("\nMcNemar-Bowker test of symmetry: X² = %.2f, df = %d, p = %.4f\n",
            mb_test$statistic, mb_test$df, mb_test$p.value))

# Also calculate marginal homogeneity (Stuart-Maxwell) using exact approach
# This tests whether the marginal distributions differ
marginal_control <- table(paired_data$Control)
marginal_case <- table(paired_data$Case)

cat("\nMarginal distributions:\n")
cat("Control: "); print(marginal_control)
cat("Case:    "); print(marginal_case)

# Simple sign test approach for direction: 
# Among pairs where phenotype differs, which direction is more common?
discordant <- paired_data %>%
  filter(Control != Case) %>%
  mutate(
    direction = case_when(
      # "Worse" = going from Infiltrated to Desert/Excluded, or Excluded to Desert
      (Control == "Infiltrated" & Case %in% c("Excluded", "Desert")) ~ "Worse",
      (Control == "Excluded" & Case == "Desert") ~ "Worse",
      # "Better" = opposite direction
      (Case == "Infiltrated" & Control %in% c("Excluded", "Desert")) ~ "Better",
      (Case == "Excluded" & Control == "Desert") ~ "Better",
      TRUE ~ "Lateral"  # Excluded <-> other at same "level"
    )
  )

cat("\nAmong discordant pairs:\n")
print(table(discordant$direction))

# Binomial test on Infiltrated specifically
infiltrated_control_not_case <- sum(paired_data$Control == "Infiltrated" & paired_data$Case != "Infiltrated")
infiltrated_case_not_control <- sum(paired_data$Case == "Infiltrated" & paired_data$Control != "Infiltrated")
cat(sprintf("\nInfiltrated in Control but not Case: %d\n", infiltrated_control_not_case))
cat(sprintf("Infiltrated in Case but not Control: %d\n", infiltrated_case_not_control))

if ((infiltrated_control_not_case + infiltrated_case_not_control) > 0) {
  binom_infiltrated <- binom.test(infiltrated_control_not_case, 
                                  infiltrated_control_not_case + infiltrated_case_not_control)
  cat(sprintf("McNemar test for Infiltrated: p = %.4f\n", binom_infiltrated$p.value))
}

# Desert specifically
desert_case_not_control <- sum(paired_data$Case == "Desert" & paired_data$Control != "Desert")
desert_control_not_case <- sum(paired_data$Control == "Desert" & paired_data$Case != "Desert")
cat(sprintf("\nDesert in Case but not Control: %d\n", desert_case_not_control))
cat(sprintf("Desert in Control but not Case: %d\n", desert_control_not_case))

if ((desert_case_not_control + desert_control_not_case) > 0) {
  binom_desert <- binom.test(desert_case_not_control, 
                             desert_case_not_control + desert_control_not_case)
  cat(sprintf("McNemar test for Desert: p = %.4f\n", binom_desert$p.value))
}

# Use paired test p-value for final result
final_p <- mb_test$p.value
final_p_infiltrated <- if(exists("binom_infiltrated")) binom_infiltrated$p.value else NA
final_p_desert <- if(exists("binom_desert")) binom_desert$p.value else NA


# =============================================================================
# PART 5: PUBLICATION-READY FIGURES
# =============================================================================

# --- Figure A: Stacked bar chart ---
plot_data <- phenotype_data %>%
  count(Group, Phenotype) %>%
  group_by(Group) %>%
  mutate(Percentage = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    Group = factor(Group, levels = c("Control", "Case"))
  )

p_bar <- ggplot(plot_data, aes(x = Group, y = Percentage, fill = Phenotype)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", Percentage)),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(
    values = c("Infiltrated" = "#2E86AB", "Excluded" = "#F6AE2D", "Desert" = "#E94F37"),
    name = "Immune phenotype"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 107)) +
  labs(
    x = NULL,
    y = "Percentage of tumors",
    title = NULL
  ) +
  annotate("text", x = 1.5, y = 103, 
           label = sprintf("p = %.2f", final_p),
           size = 3.5) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.4),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save bar chart
png("Figure_Immune_Phenotypes_Bar.png", width = 5, height = 4, units = "in", res = 300)
print(p_bar)
dev.off()

ggsave("Figure_Immune_Phenotypes_Bar.pdf", p_bar, width = 5, height = 4)


# --- Figure B: Scatter plot showing classification ---
p_scatter <- ggplot(phenotype_data, aes(x = Density, y = MNND)) +
  geom_point(aes(fill = Phenotype, shape = Group), size = 3, alpha = 0.7, stroke = 0.5) +
  geom_hline(yintercept = mnnd_median, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = density_median, linetype = "dashed", color = "grey50") +
  scale_fill_manual(
    values = c("Infiltrated" = "#2E86AB", "Excluded" = "#F6AE2D", "Desert" = "#E94F37"),
    name = "Phenotype"
  ) +
  scale_shape_manual(
    values = c("Control" = 21, "Case" = 24),
    name = "Group"
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(
    x = expression("CD8+ T-cell density (cells/mm"^2*")"),
    y = "MNND Melanoma–CD8 (µm)"
  ) +
  # Add quadrant labels
  annotate("text", x = max(phenotype_data$Density) * 0.95, y = min(phenotype_data$MNND) * 1.1,
           label = "Infiltrated", color = "#2E86AB", fontface = "bold", hjust = 1, size = 3.5) +
  annotate("text", x = max(phenotype_data$Density) * 0.95, y = max(phenotype_data$MNND) * 0.95,
           label = "Excluded", color = "#F6AE2D", fontface = "bold", hjust = 1, size = 3.5) +
  annotate("text", x = min(phenotype_data$Density) * 1.1, y = max(phenotype_data$MNND) * 0.95,
           label = "Desert", color = "#E94F37", fontface = "bold", hjust = 0, size = 3.5) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.4),
    legend.position = "right"
  )

# Save scatter plot
png("Figure_Immune_Phenotypes_Scatter.png", width = 7, height = 5, units = "in", res = 300)
print(p_scatter)
dev.off()

ggsave("Figure_Immune_Phenotypes_Scatter.pdf", p_scatter, width = 7, height = 5)


cat("\n\nFigures saved:\n")
cat("  - Figure_Immune_Phenotypes_Bar.png/pdf\n")
cat("  - Figure_Immune_Phenotypes_Scatter.png/pdf\n")
