library(dplyr)
library(ggplot2)
library(ggsignif)
library(tidyr)
library(survival)  # for clogit

# Function to calculate cell density
calculate_cell_density <- function(spe_object, 
                                   cell_types = NULL,  
                                   unit_conversion = 1000) {  
  
  # Get coordinates and Cell.Type
  coords <- SpatialExperiment::spatialCoords(spe_object)
  cells_df <- data.frame(
    Cell.X.Position = coords[,1],
    Cell.Y.Position = coords[,2],
    Cell.Type = colData(spe_object)$Cell.Type
  )
  
  # Get image boundaries
  x_max <- max(cells_df$Cell.X.Position)
  x_min <- min(cells_df$Cell.X.Position)
  y_max <- max(cells_df$Cell.Y.Position)
  y_min <- min(cells_df$Cell.Y.Position)
  
  # Convert μm² to mm²
  area_mm2 <- (x_max - x_min) * (y_max - y_min) / (unit_conversion^2)
  
  if (is.null(cell_types)) {
    # Calculate density for all cells
    total_cells <- nrow(cells_df)
    density <- total_cells / area_mm2
    
    result <- data.frame(
      Cell.Type = "All",
      Number_of_Cells = total_cells,
      Area_mm2 = area_mm2,
      Density_per_mm2 = density
    )
    
  } else {
    # Calculate density for specific cell types
    densities <- lapply(cell_types, function(cell_type) {
      n_cells <- sum(cells_df$Cell.Type == cell_type)
      density <- n_cells / area_mm2
      
      data.frame(
        Cell.Type = cell_type,
        Number_of_Cells = n_cells,
        Area_mm2 = area_mm2,
        Density_per_mm2 = density
      )
    })
    result <- do.call(rbind, densities)
  }
  
  return(result)
}

# Define cell types of interest
cell_types_of_interest <- c("T-cell", "Cytotoxic T-cell", "B-cell lineage", "Macrophage")

# Get all SPE objects
spe_objects <- ls(pattern = "^SPE.*no$")

# Create list for results
all_densities <- list()

# Calculate densities for each SPE object
for(obj_name in spe_objects) {
  spe_obj <- get(obj_name)
  densities <- calculate_cell_density(spe_obj, cell_types = cell_types_of_interest)
  densities$Sample <- obj_name
  densities$Group <- ifelse(grepl("P0", obj_name), "Control", "Case")
  densities$Patient <- sub("(SPE_\\d+)P[01].*", "\\1", obj_name)
  all_densities[[obj_name]] <- densities
}

# Combine all results
all_densities_df <- do.call(rbind, all_densities)

# Get complete pairs
complete_patients <- all_densities_df %>%
  group_by(Patient) %>%
  summarize(n_groups = n_distinct(Group)) %>%
  filter(n_groups == 2) %>%
  pull(Patient)

# Filter for complete pairs
complete_pairs_df <- all_densities_df %>%
  filter(Patient %in% complete_patients)

# Perform statistical comparison (paired Wilcoxon) with BH correction and plot-friendly labels
statistical_comparison <- lapply(cell_types_of_interest, function(cell_type) {
  subset_data <- complete_pairs_df[complete_pairs_df$Cell.Type == cell_type, ]
  
  paired_data <- subset_data %>%
    select(Patient, Group, Density_per_mm2) %>%
    pivot_wider(names_from = Group, values_from = Density_per_mm2) %>%
    drop_na()
  
  if (nrow(paired_data) < 3) {
    return(data.frame(
      Cell.Type = cell_type,
      p_value = NA_real_,
      n_pairs = nrow(paired_data)
    ))
  }
  
  test_result <- tryCatch({
    wilcox.test(paired_data$Control, paired_data$Case, paired = TRUE)
  }, error = function(e) {
    message("Error in wilcox test for ", cell_type, ": ", e$message)
    return(list(p.value = NA_real_))
  })
  
  data.frame(
    Cell.Type = cell_type,
    p_value = test_result$p.value,
    n_pairs = nrow(paired_data)
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    label = paste0("p=", signif(p_value, 3),
                   ifelse(!is.na(p_adjusted) & p_adjusted < 0.05, "*", ""))
  )

# Summary stats per group (mean/median densities)
density_summary <- complete_pairs_df %>%
  dplyr::mutate(Group = factor(Group, levels = c("Control","Case"))) %>%
  dplyr::group_by(Cell.Type, Group) %>%
  dplyr::summarise(
    mean_density   = mean(Density_per_mm2, na.rm = TRUE),
    median_density = median(Density_per_mm2, na.rm = TRUE),
    n_samples      = dplyr::n(),
    .groups = "drop"
  )

# Excel table combining stats + group summaries (one row per Cell.Type x Group)
summary_table <- density_summary %>%
  dplyr::left_join(
    statistical_comparison %>%
      dplyr::select(Cell.Type, p_value, p_adjusted, label, n_pairs),
    by = "Cell.Type"
  ) %>%
  dplyr::arrange(Cell.Type, Group)

openxlsx::write.xlsx(summary_table, file = "cell_density_stats_with_values.xlsx")

# Compute y-positions for facet labels and prepare annotation data
max_y_df <- complete_pairs_df %>%
  dplyr::group_by(Cell.Type) %>%
  dplyr::summarise(max_y = max(Density_per_mm2, na.rm = TRUE) * 1.05, .groups = "drop")

plot_annotations <- dplyr::left_join(
  statistical_comparison %>% dplyr::select(Cell.Type, label),
  max_y_df, by = "Cell.Type"
)

# Boxplot in same style as proportions plots, with p-values and BH-significance asterisk
p_density <- ggplot(
  complete_pairs_df %>% dplyr::mutate(Group = factor(Group, levels = c("Control","Case"))),
  aes(x = Group, y = Density_per_mm2, fill = Group)
) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5) +
  geom_text(
    data = plot_annotations,
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.3,
    size = 3
  ) +
  scale_fill_manual(values = c("Case" = "#E69F00", "Control" = "#56B4E9")) +
  facet_wrap(~ Cell.Type, scales = "free_y", ncol = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    x = NULL,
    y = "Density (cells/mm²)"
  ) +
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
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.margin = unit(c(0.8, 0.5, 0.5, 0.5), "cm")
  )

print(p_density)

# Save plots
ggsave("cell_density_boxplots.png", p_density, width = 8, height = 6, dpi = 300)
ggsave("cell_density_boxplots.pdf", p_density, width = 8, height = 6)

# 1. Extract cytotoxic T-cell density per sample
cd8_df <- complete_pairs_df %>%
  filter(Cell.Type == "Cytotoxic T-cell") %>%
  mutate(Group = factor(Group, levels = c("Control","Case")))

# 2. Prepare data for clogit
# Outcome: Case=1, Control=0; stratified by Pair_ID
clogit_df <- cd8_df %>%
  mutate(
    status = ifelse(Group == "Case", 1, 0),
    Pair_ID = as.factor(Patient)   # or Pair_ID depending on column names
  ) %>%
  dplyr::select(Pair_ID, status, Density_per_mm2)


# 3. Fit conditional logistic regression
m <- clogit(status ~ Density_per_mm2 + strata(Pair_ID), data = clogit_df)

summary(m)

#Coefficient: –0.0035
#Odds ratio (exp(coef)): 0.9965 per 1 cell/mm² increase
#95% CI: 0.9935–0.9995
#p = 0.024 (significant)

#Per 100 cells/mm2 increase 

clogit_df <- clogit_df %>%
  mutate(density_per100 = Density_per_mm2 / 100)

m100 <- clogit(status ~ density_per100 + strata(Pair_ID), data = clogit_df)
summary(m100)

## Per doubling of density 
clogit_df <- clogit_df %>%
  mutate(log2_density = log2(Density_per_mm2 + 1))

m_log2 <- clogit(status ~ log2_density + strata(Pair_ID), data = clogit_df)
summary(m_log2)

###
