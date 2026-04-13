# ================================
# MXIF clustering + metrics + plots (DBSCAN3)
# ================================

# ---- Packages ----

library(tidyverse)
library(dbscan)
library(writexl)
library(rstatix)
library(scales)

# ================================
# 0) Parameters / paths
# ================================
# Set your output directory — RDS files, Excel tables, and figures will be saved here
output_dir <- "your/output/directory"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

min_pts  <- 70
n_tiles  <- 6
feature_colname <- "Cell.Type"
cell_types_of_interest <- c("T-cell","B-cell lineage","Cytotoxic T-cell","Macrophage")
overwrite_rds <- TRUE  # set FALSE to skip files that already exist

# ================================
# 1) Function: process_sample_hdbscan_plus
# ================================
process_sample_hdbscan_plus <- function(spe_object,
                                        cell_types_of_interest = c("T-cell","B-cell lineage","Cytotoxic T-cell","Macrophage"),
                                        min_pts = 70,
                                        feature_colname = "Cell.Type",
                                        n_tiles = 6) {
  library(dplyr); library(tidyr); library(tibble)
  # --- Extract immune cells of interest
  cell_data <- as.data.frame(cbind(colData(spe_object), spatialCoords(spe_object)))
  cells_of_interest <- cell_data[cell_data[[feature_colname]] %in% cell_types_of_interest, ]
  if (nrow(cells_of_interest) == 0) {
    return(tibble(
      # original 26
      total_clusters = 0L, total_immune_cells = 0L, mean_cluster_size = NA_real_,
      CD8_dominant_clusters = 0L, frac_CD8_dominant_clusters = NA_real_,
      CD3_dominant_clusters = 0L, frac_CD3_dominant_clusters = NA_real_,
      Macro_dominant_clusters = 0L, frac_Macro_dominant_clusters = NA_real_,
      Bcell_dominant_clusters = 0L, frac_Bcell_dominant_clusters = NA_real_,
      mean_CD3_CD8_ratio = NA_real_,
      clusters_CD3_CD8_ratio_gt2 = 0L, clusters_CD3_CD8_ratio_lt1 = 0L,
      frac_clusters_CD3_CD8_gt2 = NA_real_, frac_clusters_CD3_CD8_lt1 = NA_real_,
      CD8_dominant_clusters_50 = 0L, frac_CD8_dominant_clusters_50 = NA_real_,
      CD8_dominant_clusters_40 = 0L, frac_CD8_dominant_clusters_40 = NA_real_,
      CD8_top_dominant_min30 = 0L, frac_CD8_top_dominant_min30 = NA_real_,
      has_mixed_cluster = FALSE, frac_clustered_cells = NA_real_,
      has_tls_candidate = FALSE, has_entropy_gt_1.2 = FALSE,
      # new
      pct_cells_in_largest_cluster = NA_real_,
      pct_cells_in_top3_clusters = NA_real_,
      largest_cluster_frac_CD8 = NA_real_,
      largest_cluster_frac_Macro = NA_real_,
      mean_entropy_across_clusters = NA_real_,
      frac_CD8_cells_clustered = NA_real_,
      frac_CD3_cells_clustered = NA_real_,
      frac_B_cells_clustered = NA_real_,
      frac_Macro_cells_clustered = NA_real_,
      frac_clusters_with_CD8_and_Macro_20_20 = NA_real_,
      frac_clusters_with_CD8_and_Macro_30_30 = NA_real_,
      # looser variants
      has_tls_candidate_loose = FALSE,
      CD8_top_dominant_min25 = 0L,
      frac_CD8_top_dominant_min25 = NA_real_,
      CD8_min30_no_top = 0L,
      frac_CD8_min30_no_top = NA_real_
    ))
  }
  
  # --- Tiled HDBSCAN (memory-safe)
  x_breaks <- stats::quantile(cells_of_interest$Cell.X.Position, probs = seq(0, 1, length.out = n_tiles + 1))
  cells_of_interest$tile <- cut(cells_of_interest$Cell.X.Position, breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
  
  full_cluster_vec <- rep(NA_character_, nrow(cells_of_interest))
  cluster_id_counter <- 0L
  
  for (tile_num in seq_len(n_tiles)) {
    tile_cells <- cells_of_interest[cells_of_interest$tile == tile_num, ]
    if (nrow(tile_cells) < min_pts) next
    coords <- tile_cells[, c("Cell.X.Position", "Cell.Y.Position")]
    result <- try(dbscan::hdbscan(coords, minPts = min_pts), silent = TRUE)
    if (inherits(result, "try-error")) next
    ids_in_tile <- which(cells_of_interest$tile == tile_num)
    shifted <- ifelse(result$cluster == 0, NA_integer_, result$cluster + cluster_id_counter)
    full_cluster_vec[ids_in_tile] <- ifelse(is.na(shifted), NA_character_, paste0("Cluster_", shifted))
    max_c <- suppressWarnings(max(shifted, na.rm = TRUE))
    if (!is.finite(max_c)) max_c <- cluster_id_counter
    cluster_id_counter <- max(cluster_id_counter, max_c)
  }
  
  # --- Attach cluster labels & extract clustered cells
  cells_of_interest$Neighborhood <- full_cluster_vec
  colData(spe_object)$Neighborhood <- NA_character_
  match_idx <- match(cells_of_interest$Cell.ID, colData(spe_object)$Cell.ID)
  colData(spe_object)$Neighborhood[match_idx] <- cells_of_interest$Neighborhood
  
  idx_keep <- !is.na(colData(spe_object)$Neighborhood)
  hdbscan_cells <- as.data.frame(cbind(colData(spe_object)[idx_keep, ], spatialCoords(spe_object)[idx_keep, ]))
  
  base_na_return <- function(frac_clustered_cells_val) {
    tibble(
      total_clusters = 0L, total_immune_cells = 0L, mean_cluster_size = NA_real_,
      CD8_dominant_clusters = 0L, frac_CD8_dominant_clusters = NA_real_,
      CD3_dominant_clusters = 0L, frac_CD3_dominant_clusters = NA_real_,
      Macro_dominant_clusters = 0L, frac_Macro_dominant_clusters = NA_real_,
      Bcell_dominant_clusters = 0L, frac_Bcell_dominant_clusters = NA_real_,
      mean_CD3_CD8_ratio = NA_real_,
      clusters_CD3_CD8_ratio_gt2 = 0L, clusters_CD3_CD8_ratio_lt1 = 0L,
      frac_clusters_CD3_CD8_gt2 = NA_real_, frac_clusters_CD3_CD8_lt1 = NA_real_,
      CD8_dominant_clusters_50 = 0L, frac_CD8_dominant_clusters_50 = NA_real_,
      CD8_dominant_clusters_40 = 0L, frac_CD8_dominant_clusters_40 = NA_real_,
      CD8_top_dominant_min30 = 0L, frac_CD8_top_dominant_min30 = NA_real_,
      has_mixed_cluster = FALSE, frac_clustered_cells = frac_clustered_cells_val,
      has_tls_candidate = FALSE, has_entropy_gt_1.2 = FALSE,
      pct_cells_in_largest_cluster = NA_real_,
      pct_cells_in_top3_clusters = NA_real_,
      largest_cluster_frac_CD8 = NA_real_,
      largest_cluster_frac_Macro = NA_real_,
      mean_entropy_across_clusters = NA_real_,
      frac_CD8_cells_clustered = NA_real_,
      frac_CD3_cells_clustered = NA_real_,
      frac_B_cells_clustered = NA_real_,
      frac_Macro_cells_clustered = NA_real_,
      frac_clusters_with_CD8_and_Macro_20_20 = NA_real_,
      frac_clusters_with_CD8_and_Macro_30_30 = NA_real_,
      has_tls_candidate_loose = FALSE,
      CD8_top_dominant_min25 = 0L,
      frac_CD8_top_dominant_min25 = NA_real_,
      CD8_min30_no_top = 0L,
      frac_CD8_min30_no_top = NA_real_
    )
  }
  
  if (nrow(hdbscan_cells) == 0) {
    return(base_na_return(frac_clustered_cells_val = nrow(hdbscan_cells) / nrow(cells_of_interest)))
  }
  
  # --- Per-cluster composition
  required_types <- c("Cytotoxic T-cell","T-cell","B-cell lineage","Macrophage")
  
  counts_wide <- hdbscan_cells %>%
    group_by(Neighborhood, Cell.Type) %>%
    summarise(n = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Cell.Type, values_from = n, values_fill = 0, names_repair = "minimal")
  
  for (ct in required_types) if (!ct %in% names(counts_wide)) counts_wide[[ct]] <- 0
  
  counts_wide <- counts_wide %>%
    mutate(
      total_cells = `Cytotoxic T-cell` + `T-cell` + `B-cell lineage` + `Macrophage`,
      frac_CD8    = if_else(total_cells == 0, NA_real_, `Cytotoxic T-cell` / total_cells),
      frac_CD3    = if_else(total_cells == 0, NA_real_, `T-cell` / total_cells),
      frac_B      = if_else(total_cells == 0, NA_real_, `B-cell lineage` / total_cells),
      frac_Macro  = if_else(total_cells == 0, NA_real_, `Macrophage` / total_cells),
      entropy     = apply(cbind(`Cytotoxic T-cell`, `T-cell`, `B-cell lineage`, `Macrophage`), 1, function(x) {
        s <- sum(x); if (s == 0) return(NA_real_)
        p <- x / s; -sum(ifelse(p > 0, p * log2(p), 0))
      }),
      nB = `B-cell lineage`,
      nT = `T-cell` + `Cytotoxic T-cell`,
      fT = ifelse(total_cells > 0, nT/total_cells, NA_real_),
      tls_loose = (frac_B >= 0.10) & (fT >= 0.20),
      CD8_top_min25 = (frac_CD8 >= 0.25) & (frac_CD8 >= pmax(frac_CD3, frac_B, frac_Macro, na.rm = TRUE)),
      CD8_min30_no_top = (frac_CD8 > 0.30)
    )
  
  # --- Original 26 metrics
  total_clusters     <- nrow(counts_wide)
  total_immune_cells <- sum(counts_wide$total_cells)
  mean_cluster_size  <- mean(counts_wide$total_cells)
  
  CD8_dom_60   <- sum(counts_wide$frac_CD8   > 0.60, na.rm = TRUE)
  CD3_dom_60   <- sum(counts_wide$frac_CD3   > 0.60, na.rm = TRUE)
  Macro_dom_50 <- sum(counts_wide$frac_Macro > 0.50, na.rm = TRUE)
  Bcell_dom_10 <- sum(counts_wide$frac_B     > 0.10, na.rm = TRUE)
  
  cd3_cd8_ratio <- (counts_wide$`T-cell` + 1) / (counts_wide$`Cytotoxic T-cell` + 1)
  mean_ratio <- mean(cd3_cd8_ratio, na.rm = TRUE)
  ratio_gt2  <- sum(cd3_cd8_ratio > 2, na.rm = TRUE)
  ratio_lt1  <- sum(cd3_cd8_ratio < 1, na.rm = TRUE)
  
  CD8_dom_50 <- sum(counts_wide$frac_CD8 > 0.50, na.rm = TRUE)
  CD8_dom_40 <- sum(counts_wide$frac_CD8 > 0.40, na.rm = TRUE)
  
  top_flag <- (counts_wide$frac_CD8 >= pmax(counts_wide$frac_CD8, counts_wide$frac_CD3, counts_wide$frac_Macro, counts_wide$frac_B)) &
    (counts_wide$frac_CD8 > 0.30)
  CD8_top_min30 <- sum(top_flag, na.rm = TRUE)
  
  denom <- if (total_clusters > 0) total_clusters else NA_real_
  
  has_mixed_cluster <- any(rowSums(cbind(counts_wide$frac_CD8, counts_wide$frac_CD3, counts_wide$frac_B, counts_wide$frac_Macro) > 0.25, na.rm = TRUE) >= 3)
  frac_clustered_cells <- nrow(hdbscan_cells) / nrow(cells_of_interest)
  has_tls_candidate <- any((counts_wide$frac_CD3 > 0.20 | counts_wide$frac_CD8 > 0.20) & counts_wide$frac_B > 0.10, na.rm = TRUE)
  has_entropy_gt_1.2 <- any(counts_wide$entropy > 1.2, na.rm = TRUE)
  
  # --- New derived metrics
  ord <- order(counts_wide$total_cells, decreasing = TRUE)
  top1_total <- if (length(ord) >= 1) counts_wide$total_cells[ord[1]] else NA_real_
  top3_total <- if (length(ord) >= 3) sum(counts_wide$total_cells[ord[1:3]]) else sum(counts_wide$total_cells[ord], na.rm = TRUE)
  pct_cells_in_largest_cluster <- if (total_immune_cells > 0) top1_total / total_immune_cells else NA_real_
  pct_cells_in_top3_clusters   <- if (total_immune_cells > 0) top3_total / total_immune_cells else NA_real_
  
  largest_cluster_frac_CD8   <- if (length(ord) >= 1) counts_wide$frac_CD8[ord[1]] else NA_real_
  largest_cluster_frac_Macro <- if (length(ord) >= 1) counts_wide$frac_Macro[ord[1]] else NA_real_
  
  mean_entropy_across_clusters <- mean(counts_wide$entropy, na.rm = TRUE)
  
  type_total_all <- cells_of_interest %>%
    count(Cell.Type, name = "total_all") %>%
    tidyr::complete(Cell.Type = required_types, fill = list(total_all = 0))
  type_total_clustered <- hdbscan_cells %>%
    count(Cell.Type, name = "total_clust") %>%
    tidyr::complete(Cell.Type = required_types, fill = list(total_clust = 0))
  per_type <- type_total_all %>%
    left_join(type_total_clustered, by = "Cell.Type") %>%
    mutate(frac_clustered = ifelse(total_all > 0, total_clust / total_all, NA_real_))
  frac_CD8_cells_clustered   <- per_type$frac_clustered[per_type$Cell.Type == "Cytotoxic T-cell"]
  frac_CD3_cells_clustered   <- per_type$frac_clustered[per_type$Cell.Type == "T-cell"]
  frac_B_cells_clustered     <- per_type$frac_clustered[per_type$Cell.Type == "B-cell lineage"]
  frac_Macro_cells_clustered <- per_type$frac_clustered[per_type$Cell.Type == "Macrophage"]
  
  both_20_20 <- mean((counts_wide$frac_CD8 > 0.20) & (counts_wide$frac_Macro > 0.20), na.rm = TRUE)
  both_30_30 <- mean((counts_wide$frac_CD8 > 0.30) & (counts_wide$frac_Macro > 0.30), na.rm = TRUE)
  
  # Looser variants
  has_tls_candidate_loose <- any(counts_wide$tls_loose, na.rm = TRUE)
  CD8_top_min25_cnt <- sum(counts_wide$CD8_top_min25, na.rm = TRUE)
  CD8_top_min25_frac <- CD8_top_min25_cnt / denom
  CD8_min30_no_top_cnt <- sum(counts_wide$CD8_min30_no_top, na.rm = TRUE)
  CD8_min30_no_top_frac <- CD8_min30_no_top_cnt / denom
  
  # --- Return combined tibble
  tibble(
    total_clusters = total_clusters,
    total_immune_cells = total_immune_cells,
    mean_cluster_size = mean_cluster_size,
    CD8_dominant_clusters = CD8_dom_60,
    frac_CD8_dominant_clusters = CD8_dom_60 / denom,
    CD3_dominant_clusters = CD3_dom_60,
    frac_CD3_dominant_clusters = CD3_dom_60 / denom,
    Macro_dominant_clusters = Macro_dom_50,
    frac_Macro_dominant_clusters = Macro_dom_50 / denom,
    Bcell_dominant_clusters = Bcell_dom_10,
    frac_Bcell_dominant_clusters = Bcell_dom_10 / denom,
    mean_CD3_CD8_ratio = mean_ratio,
    clusters_CD3_CD8_ratio_gt2 = ratio_gt2,
    clusters_CD3_CD8_ratio_lt1 = ratio_lt1,
    frac_clusters_CD3_CD8_gt2 = ratio_gt2 / denom,
    frac_clusters_CD3_CD8_lt1 = ratio_lt1 / denom,
    CD8_dominant_clusters_50 = CD8_dom_50,
    frac_CD8_dominant_clusters_50 = CD8_dom_50 / denom,
    CD8_dominant_clusters_40 = CD8_dom_40,
    frac_CD8_dominant_clusters_40 = CD8_dom_40 / denom,
    CD8_top_dominant_min30 = CD8_top_min30,
    frac_CD8_top_dominant_min30 = CD8_top_min30 / denom,
    has_mixed_cluster = has_mixed_cluster,
    frac_clustered_cells = frac_clustered_cells,
    has_tls_candidate = has_tls_candidate,
    has_entropy_gt_1.2 = has_entropy_gt_1.2,
    pct_cells_in_largest_cluster = pct_cells_in_largest_cluster,
    pct_cells_in_top3_clusters = pct_cells_in_top3_clusters,
    largest_cluster_frac_CD8 = largest_cluster_frac_CD8,
    largest_cluster_frac_Macro = largest_cluster_frac_Macro,
    mean_entropy_across_clusters = mean_entropy_across_clusters,
    frac_CD8_cells_clustered = frac_CD8_cells_clustered,
    frac_CD3_cells_clustered = frac_CD3_cells_clustered,
    frac_B_cells_clustered = frac_B_cells_clustered,
    frac_Macro_cells_clustered = frac_Macro_cells_clustered,
    frac_clusters_with_CD8_and_Macro_20_20 = both_20_20,
    frac_clusters_with_CD8_and_Macro_30_30 = both_30_30,
    has_tls_candidate_loose = has_tls_candidate_loose,
    CD8_top_dominant_min25 = CD8_top_min25_cnt,
    frac_CD8_top_dominant_min25 = CD8_top_min25_frac,
    CD8_min30_no_top = CD8_min30_no_top_cnt,
    frac_CD8_min30_no_top = CD8_min30_no_top_frac
  )
}

# ================================
# 2) Create RDS files in dbscan3
# ================================
# ================================
# Process ONLY complete matched sets, memory-efficiently (DBSCAN3)
# ================================

# Assumes you already defined:
# - process_sample_hdbscan_plus()
# - output_dir, cell_types_of_interest, min_pts, feature_colname, n_tiles
# If not, define them before running this block.


# --- 0) Helper that processes exactly ONE sample name and frees memory
process_one_sample <- function(spe_name) {
  sample_id <- sub("^SPE_", "", spe_name)
  out_file  <- file.path(output_dir, paste0(sample_id, ".rds"))
  if (!overwrite_rds && file.exists(out_file)) {
    message("Skipping (exists): ", sample_id)
    return(invisible(TRUE))
  }
  
  if (!exists(spe_name, inherits = TRUE)) {
    warning("Missing in memory: ", spe_name)
    return(invisible(FALSE))
  }
  
  message("Processing: ", sample_id)
  # Load a reference; avoid extra copies
  spe <- get(spe_name, inherits = TRUE)
  
  # Do work inside a local scope so temporaries die as soon as we exit
  res <- try({
    process_sample_hdbscan_plus(
      spe_object = spe,
      cell_types_of_interest = cell_types_of_interest,
      min_pts = min_pts,
      feature_colname = feature_colname,
      n_tiles = n_tiles
    )
  }, silent = TRUE)
  
  # Immediately drop SPE (largest object) before any further work
  rm(spe); gc()
  
  if (inherits(res, "try-error") || is.null(res)) {
    warning("Failed on sample: ", sample_id)
    rm(list = c("res")); gc()
    return(invisible(FALSE))
  }
  
  # Save and clear
  res$sample_id <- sample_id
  saveRDS(res, out_file)
  rm(res); gc()
  
  invisible(TRUE)
}

# --- 1) Find only COMPLETE pairs among currently loaded SPE_* objects
all_spe_names <- ls(pattern = "^SPE_")
stopifnot(length(all_spe_names) > 0)

loaded_map <- tibble::tibble(
  spe_name  = all_spe_names,
  sample_id = sub("^SPE_", "", all_spe_names)
) %>%
  dplyr::mutate(
    SetID = sub("P[01].*$", "", sample_id),
    Group = dplyr::case_when(
      grepl("P0", sample_id) ~ "control",
      grepl("P1", sample_id) ~ "case",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Group))

complete_setids <- loaded_map %>%
  dplyr::group_by(SetID) %>%
  dplyr::filter(dplyr::n_distinct(Group) == 2) %>%
  dplyr::ungroup() %>%
  dplyr::pull(SetID) %>%
  unique()

if (length(complete_setids) == 0) {
  stop("No complete matched sets among currently loaded SPE_* objects.")
}

# Only those SPEs that belong to a complete matched SetID
spe_to_process <- loaded_map %>%
  dplyr::filter(SetID %in% complete_setids) %>%
  dplyr::arrange(SetID, Group) %>%     # (optional) control then case per set
  dplyr::pull(spe_name)

message("Complete matched sets: ", length(complete_setids))
message("SPE objects to process: ", length(spe_to_process))

# --- 2) Process one-by-one, freeing memory every iteration
ok  <- logical(length(spe_to_process))
for (i in seq_along(spe_to_process)) {
  ok[i] <- process_one_sample(spe_to_process[i])
  # Defensive: drop any stray large objects created by other code
  # Keep only functions/params and minimal environment stuff
  # (commented out by default to avoid removing things you need)
  # keep <- c("process_one_sample", "process_sample_hdbscan_plus",
  #           "output_dir","cell_types_of_interest","min_pts","feature_colname","n_tiles",
  #           "overwrite_rds","all_spe_names","loaded_map","complete_setids","spe_to_process","ok")
  # rm(list = setdiff(ls(), keep)); gc()
}

message("Done writing RDS for complete sets to: ", output_dir)
message("Success: ", sum(ok), " / ", length(ok))

# ================================
# 3) Load ONLY dbscan3 + build metrics table
# ================================
files <- list.files(output_dir, pattern = "\\.rds$", full.names = TRUE)

stopifnot(length(files) > 0)

metrics <- files %>%
  setNames(nm = basename(.)) %>%
  purrr::map_dfr(
    ~ {
      obj <- readRDS(.x)
      # ensure a tibble row-bindable object
      if (is.data.frame(obj)) tibble::as_tibble(obj) else tibble::tibble(obj = list(obj))
    },
    .id = "filename"
  ) %>%
  dplyr::mutate(
    sample_id = sub("\\.rds$", "", filename),
    SetID     = sub("P[01].*$", "", sample_id),
    Group     = dplyr::if_else(grepl("P0", sample_id), "Control", "Case")
  ) %>%
  dplyr::select(-filename) %>%
  dplyr::relocate(sample_id, SetID, Group)

names(metrics)

# ================================
# 4) Plot EVERYTHING (2 plots/metric) — safe column names, no jitter
# ================================


# Where to save plots
plots_dir <- file.path(output_dir, "plots_plus")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Ensure IDs & metric columns exist
id_cols <- c("sample_id","SetID","Group")
all_metric_cols <- setdiff(names(metrics), id_cols)

# Make Group/SetID factors once
metrics <- metrics %>%
  dplyr::mutate(
    Group = factor(Group, levels = c("control","case")),
    SetID = factor(SetID)
  )

# ---------- helpers ----------
coerce_to_numeric <- function(v) {
  if (is.numeric(v)) return(v)
  if (is.logical(v)) return(as.numeric(v))
  if (is.factor(v))  return(coerce_to_numeric(as.character(v)))
  if (is.character(v)) {
    vv <- trimws(tolower(v))
    out <- suppressWarnings(as.numeric(v))
    tf_mask <- vv %in% c("true","false","1","0")
    if (any(tf_mask) && all(vv[tf_mask] %in% c("true","false","1","0"))) {
      out2 <- rep(NA_real_, length(vv))
      out2[vv == "true" | vv == "1"]  <- 1
      out2[vv == "false" | vv == "0"] <- 0
      if (sum(!is.na(out2)) >= sum(!is.na(out))) return(out2)
    }
    return(out)
  }
  suppressWarnings(as.numeric(v))
}

is_binary_vec <- function(x) {
  xn <- coerce_to_numeric(x)
  ux <- unique(na.omit(xn))
  length(ux) > 0 && all(ux %in% c(0,1))
}

binary_metrics <- all_metric_cols[vapply(metrics[all_metric_cols], is_binary_vec, logical(1))]
log_scale_candidates <- c("mean_CD3_CD8_ratio")
log_scale_metrics <- intersect(log_scale_candidates, all_metric_cols)

plot_box <- function(df, y, is_binary = FALSE) {
  ynum <- coerce_to_numeric(df[[y]])
  if (all(is.na(ynum))) return(NULL)
  df2 <- df %>% dplyr::mutate(.Y = ynum, Group = factor(Group, levels = c("control","case")))
  
  p <- ggplot(df2, aes(x = Group, y = .Y, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.55) +
    geom_point(position = position_dodge(width = 0.45), alpha = 0.6, shape = 21, stroke = 0.2)
  
  if (is_binary) {
    p <- p + scale_y_continuous(breaks = c(0,1), labels = c("FALSE","TRUE"), limits = c(-0.1, 1.1))
  }
  
  p + labs(title = paste0(y, " (boxplot)"), x = NULL, y = y) + theme_bw()
}

plot_pairs <- function(df, y, logy = FALSE) {
  ynum <- coerce_to_numeric(df[[y]])
  if (all(is.na(ynum))) return(NULL)
  df2 <- df %>% dplyr::mutate(.Y = ynum, Group = factor(Group, levels = c("control","case")))
  
  p <- ggplot(df2, aes(x = Group, y = .Y, group = SetID)) +
    geom_line(alpha = 0.35) +
    geom_point(aes(color = Group), size = 2, position = position_dodge(width = 0.35)) +
    labs(title = paste0(y, " (matched pairs)"), x = NULL, y = y) +
    theme_bw()
  
  if (logy && all(df2$.Y > 0, na.rm = TRUE)) p + scale_y_log10() else p
}

plot_prop <- function(df, y) {
  ynum <- coerce_to_numeric(df[[y]])
  if (all(is.na(ynum))) return(NULL)
  df2 <- df %>% dplyr::mutate(.Y = ynum, Group = factor(Group, levels = c("control","case")))
  
  df2 %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(prop = mean(.Y, na.rm = TRUE), .groups = "drop") %>%
    ggplot(aes(x = Group, y = prop, fill = Group)) +
    geom_col(width = 0.55) +
    scale_y_continuous(labels = percent, limits = c(0, 1)) +
    labs(title = paste0(y, " (group proportion)"), x = NULL, y = "% samples") +
    theme_bw()
}

# ---------- generate & save ----------
failed <- list()
saved  <- character(0)

for (m in all_metric_cols) {
  is_bin <- m %in% binary_metrics
  logy   <- m %in% log_scale_metrics
  
  p1 <- try(plot_box(metrics, m, is_binary = is_bin), silent = TRUE)
  p2 <- try(if (is_bin) plot_prop(metrics, m) else plot_pairs(metrics, m, logy = logy), silent = TRUE)
  
  if (!inherits(p1, "try-error") && inherits(p1, "ggplot")) {
    ggsave(filename = file.path(plots_dir, paste0(m, "_box.png")), plot = p1,
           width = 6.5, height = 4.5, dpi = 300)
  } else {
    failed[[m]] <- c(failed[[m]], "box")
  }
  
  if (!inherits(p2, "try-error") && inherits(p2, "ggplot")) {
    outname <- if (is_bin) paste0(m, "_prop.png") else paste0(m, "_pairs.png")
    ggsave(filename = file.path(plots_dir, outname), plot = p2,
           width = 6.5, height = 4.5, dpi = 300)
  } else {
    failed[[m]] <- c(failed[[m]], if (is_bin) "prop" else "pairs")
  }
  
  if (!length(failed[[m]])) saved <- c(saved, m)
}

cat("Saved plots for", length(saved), "metrics into:\n  ", plots_dir, "\n")

# 5) Excel summary + matched-pair stats (raw + BH-adjusted p-values)
# ================================

# Safety: helper present

# Identify binary metrics (logical or 0/1 after coercion)
binary_metrics <- all_metric_cols[vapply(metrics[all_metric_cols], is_binary_vec, logical(1))]
cont_like      <- setdiff(all_metric_cols, binary_metrics)

# ---- Group summaries (mean, sd, median, IQR, n) ----
metrics_summary <- metrics %>%
  group_by(Group) %>%
  summarise(
    across(
      .cols = all_of(all_metric_cols),
      .fns = list(
        mean   = ~mean(coerce_to_numeric(.x), na.rm = TRUE),
        sd     = ~sd(coerce_to_numeric(.x),   na.rm = TRUE),
        median = ~median(coerce_to_numeric(.x), na.rm = TRUE),
        IQR    = ~IQR(coerce_to_numeric(.x),    na.rm = TRUE),
        n      = ~sum(!is.na(coerce_to_numeric(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# ---- Restrict to complete matched pairs ----
pairs_df <- metrics %>%
  group_by(SetID) %>% filter(n() == 2) %>% ungroup()

# ---- Paired Wilcoxon for continuous/prop metrics ----
wilcox_results <- map_df(cont_like, function(met) {
  tmp <- pairs_df %>%
    mutate(value = coerce_to_numeric(.data[[met]])) %>%
    filter(!is.na(value))
  
  # Need both groups per pair
  tmp2 <- tmp %>%
    group_by(SetID) %>% filter(n_distinct(Group) == 2) %>% ungroup()
  
  # Summaries per group
  gsum <- tmp2 %>%
    group_by(Group) %>%
    summarise(
      mean  = mean(value, na.rm = TRUE),
      median= median(value, na.rm = TRUE),
      n     = sum(!is.na(value)),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = Group,
                       values_from = c(mean, median, n),
                       names_sep = ".")
  
  # Paired test
  p_raw <- tryCatch({
    tmp2 %>% wilcox_test(value ~ Group, paired = TRUE) %>% pull(p)
  }, error = function(e) NA_real_)
  
  # Effect: median difference (case - control) on paired diffs
  diffs <- tmp2 %>% arrange(SetID, Group) %>%
    group_by(SetID) %>%
    summarise(diff = value[Group == "case"] - value[Group == "control"],
              .groups = "drop")
  med_diff <- median(diffs$diff, na.rm = TRUE)
  n_pairs  <- sum(!is.na(diffs$diff))
  
  tibble(
    metric      = met,
    test        = "wilcoxon_paired",
    n_pairs     = n_pairs,
    mean_control  = gsum$`mean.control`,
    mean_case     = gsum$`mean.case`,
    median_control= gsum$`median.control`,
    median_case   = gsum$`median.case`,
    delta_median  = med_diff,
    p_value     = p_raw
  )
})

# ---- McNemar for binary metrics ----
mcnemar_results <- map_df(binary_metrics, function(met) {
  tmp <- pairs_df %>%
    mutate(val = coerce_to_numeric(.data[[met]])) %>%
    filter(!is.na(val)) %>%
    arrange(SetID, Group)
  
  # Build 2x2 table on complete pairs
  wide <- tmp %>%
    group_by(SetID) %>%
    filter(n_distinct(Group) == 2) %>%
    summarise(
      control = as.integer(val[Group == "control"]),
      case    = as.integer(val[Group == "case"]),
      .groups = "drop"
    )
  
  if (nrow(wide) == 0) {
    return(tibble(
      metric = met, test = "mcnemar", n_pairs = 0,
      prop_control = NA_real_, prop_case = NA_real_,
      b_01 = NA_integer_, c_10 = NA_integer_,
      p_value = NA_real_
    ))
  }
  
  tbl <- table(wide$control, wide$case)
  # discordant cells
  b_01 <- ifelse(ncol(tbl) > 1, tbl[1,2], 0)
  c_10 <- ifelse(nrow(tbl) > 1, tbl[2,1], 0)
  
  p_raw <- tryCatch(mcnemar.test(tbl)$p.value, error = function(e) NA_real_)
  
  tibble(
    metric       = met,
    test         = "mcnemar",
    n_pairs      = nrow(wide),
    prop_control = mean(wide$control, na.rm = TRUE),
    prop_case    = mean(wide$case,    na.rm = TRUE),
    b_01         = as.integer(b_01),   # control=0 -> case=1
    c_10         = as.integer(c_10),   # control=1 -> case=0
    p_value      = p_raw
  )
})

# ---- Combine and add BH-adjusted p-values ----
stats_results <- bind_rows(wilcox_results, mcnemar_results) %>%
  arrange(p_value) %>%
  mutate(
    p_value_adj = p.adjust(p_value, method = "BH")
  )

# Optional: nicer ordering by test then p_adj
stats_results <- stats_results %>%
  arrange(test, p_value_adj)

# ---- Write Excel (two sheets) ----
write_xlsx(
  list(
    metrics_summary = metrics_summary,
    stats_results   = stats_results
  ),
  path = file.path(output_dir, "metrics_summary_plus.xlsx")
)

### samples with no clusters at all 
metrics %>%
  filter(coerce_to_numeric(total_clusters) == 0) %>%
  count(Group)

pairs_df %>%
  mutate(no_clusters = coerce_to_numeric(total_clusters) == 0) %>%
  arrange(SetID, Group) %>%
  group_by(SetID) %>%
  summarise(
    control = as.integer(no_clusters[Group == "control"]),
    case    = as.integer(no_clusters[Group == "case"]),
    .groups = "drop"
  ) %>%
  { mcnemar.test(table(.$control, .$case)) }

# 1) Add binary 'no_clusters' metric
metrics <- metrics %>%
  mutate(no_clusters = coerce_to_numeric(total_clusters) == 0)

# 2) Update metric lists
all_metric_cols <- setdiff(names(metrics), c("sample_id", "SetID", "Group"))
binary_metrics <- unique(c(binary_metrics, "no_clusters"))

# 3) Summary table
metrics_summary <- metrics %>%
  group_by(Group) %>%
  summarise(
    across(
      .cols = all_of(all_metric_cols),
      .fns = list(
        mean   = ~mean(coerce_to_numeric(.x), na.rm = TRUE),
        sd     = ~sd(coerce_to_numeric(.x), na.rm = TRUE),
        median = ~median(coerce_to_numeric(.x), na.rm = TRUE),
        IQR    = ~IQR(coerce_to_numeric(.x), na.rm = TRUE),
        n      = ~sum(!is.na(coerce_to_numeric(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# 4) Matched pairs
pairs_df <- metrics %>% group_by(SetID) %>% filter(n() == 2) %>% ungroup()

# Wilcoxon
cont_like <- setdiff(all_metric_cols, binary_metrics)
wilcox_results <- map_df(cont_like, ~{
  tmp <- pairs_df %>% mutate(..y.. = coerce_to_numeric(.data[[.x]]))
  tibble(
    metric = .x,
    test = "wilcoxon_paired",
    n_pairs = n_distinct(tmp$SetID[!is.na(tmp$..y..)]),
    mean_control = mean(tmp$..y..[tmp$Group == "control"], na.rm = TRUE),
    mean_case    = mean(tmp$..y..[tmp$Group == "case"], na.rm = TRUE),
    median_control = median(tmp$..y..[tmp$Group == "control"], na.rm = TRUE),
    median_case    = median(tmp$..y..[tmp$Group == "case"], na.rm = TRUE),
    delta_median   = median(tmp$..y..[tmp$Group == "case"], na.rm = TRUE) -
      median(tmp$..y..[tmp$Group == "control"], na.rm = TRUE),
    p_value = tryCatch(
      tmp %>% wilcox_test(..y.. ~ Group, paired = TRUE) %>% pull(p),
      error = function(e) NA_real_
    )
  )
})

# McNemar
mcnemar_results <- map_df(binary_metrics, ~{
  tmp <- pairs_df %>% mutate(..y.. = coerce_to_numeric(.data[[.x]]))
  dat <- tmp %>%
    arrange(SetID, Group) %>%
    group_by(SetID) %>%
    summarise(
      control = as.integer(..y..[Group == "control"]),
      case    = as.integer(..y..[Group == "case"]),
      .groups = "drop"
    )
  if (nrow(dat) == 0) return(tibble())
  tbl <- table(dat$control, dat$case)
  mc_p <- tryCatch(mcnemar.test(tbl)$p.value, error = function(e) NA_real_)
  tibble(
    metric = .x,
    test = "mcnemar",
    n_pairs = nrow(dat),
    prop_control = mean(dat$control, na.rm = TRUE),
    prop_case    = mean(dat$case, na.rm = TRUE),
    b_01 = tbl["0","1"],  # control=0, case=1
    c_10 = tbl["1","0"],  # control=1, case=0
    p_value = mc_p
  )
})

# 5) Merge and adjust
stats_results <- bind_rows(wilcox_results, mcnemar_results) %>%
  arrange(p_value) %>%
  mutate(p_value_adj = p.adjust(p_value, method = "fdr"))

# 6) Write to Excel
write_xlsx(
  list(metrics_summary = metrics_summary, stats_results = stats_results),
  path = file.path(output_dir, "metrics_summary_plus.xlsx")
)


## ADD DENSITY TO CLUSTERS 

# Clustered vs dispersed immune cell DENSITY (cells/mm²)
# - computes, per sample & cell type: density inside clusters vs outside clusters
# - keeps complete 1:1 pairs
# - paired Wilcoxon (Case vs Control) with BH correction
# - publication-style boxplots with p-value labels
# - Excel summary (means/medians + stats)


# --- SETTINGS ---
cell_types_of_interest <- c("Cytotoxic T-cell","T-cell","B-cell lineage","Macrophage")
min_pts  <- 70    # use your HDBSCAN setting
n_tiles  <- 6     # use your tiling setting
feat_col <- "Cell.Type"

# ---- Helpers ----
.area_mm2 <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  ((max(x) - min(x)) * (max(y) - min(y))) / (1000^2)
}

# Ensure Neighborhood exists on SPE; if missing, run tiled HDBSCAN and attach labels
.ensure_clusters <- function(spe,
                             cell_types_of_interest = c("Cytotoxic T-cell","T-cell","B-cell lineage","Macrophage"),
                             min_pts = 70,
                             feature_colname = "Cell.Type",
                             n_tiles = 6) {
  cd <- as.data.frame(colData(spe))
  if (!"Neighborhood" %in% colnames(cd) || length(cd$Neighborhood) != ncol(spe)) {
    colData(spe)$Neighborhood <- rep(NA_character_, ncol(spe))
  } else {
    # already present
    return(spe)
  }
  
  coords <- SpatialExperiment::spatialCoords(spe)
  stopifnot(feature_colname %in% names(cd))
  
  cells <- data.frame(
    Cell.ID = cd$Cell.ID,
    Cell.Type = cd[[feature_colname]],
    Cell.X.Position = coords[,1],
    Cell.Y.Position = coords[,2],
    stringsAsFactors = FALSE
  )
  
  coi <- cells[cells$Cell.Type %in% cell_types_of_interest, ]
  if (nrow(coi) == 0) return(spe)
  
  # Tile along X to keep memory bounded
  x_breaks <- stats::quantile(coi$Cell.X.Position, probs = seq(0, 1, length.out = n_tiles + 1), na.rm = TRUE)
  coi$tile <- cut(coi$Cell.X.Position, breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
  
  full_cluster_vec <- rep(NA_integer_, nrow(coi))
  cluster_id_counter <- 0L
  
  for (tile_num in seq_len(n_tiles)) {
    tile_cells <- coi[coi$tile == tile_num, ]
    if (nrow(tile_cells) < min_pts) next
    pts <- as.matrix(tile_cells[, c("Cell.X.Position","Cell.Y.Position")])
    res <- try(dbscan::hdbscan(pts, minPts = min_pts), silent = TRUE)
    if (inherits(res, "try-error")) next
    shifted <- ifelse(res$cluster == 0, NA_integer_, res$cluster + cluster_id_counter)
    full_cluster_vec[coi$tile == tile_num] <- shifted
    mx <- suppressWarnings(max(shifted, na.rm = TRUE))
    if (is.finite(mx)) cluster_id_counter <- max(cluster_id_counter, mx)
  }
  
  labels <- ifelse(is.na(full_cluster_vec), NA_character_, paste0("Cluster_", full_cluster_vec))
  idx <- match(coi$Cell.ID, cd$Cell.ID)
  colData(spe)$Neighborhood[idx] <- labels
  spe
}

# ---- 1) Build clustered vs dispersed density table per sample ----
spe_objects <- ls(pattern = "^SPE.*no$")
stopifnot(length(spe_objects) > 0)

dens_list <- lapply(spe_objects, function(nm) {
  spe <- get(nm)
  spe <- .ensure_clusters(
    spe,
    cell_types_of_interest = cell_types_of_interest,
    min_pts = min_pts,
    feature_colname = feat_col,
    n_tiles = n_tiles
  )
  
  coords <- SpatialExperiment::spatialCoords(spe)
  cd     <- as.data.frame(colData(spe))
  stopifnot(all(c("Cell.ID", feat_col, "Neighborhood") %in% names(cd)))
  
  # Whole-slide area from ALL immune cells of interest
  keep <- cd[[feat_col]] %in% cell_types_of_interest
  if (!any(keep)) return(NULL)
  
  x_all <- coords[keep, 1]; y_all <- coords[keep, 2]
  area_ws <- .area_mm2(x_all, y_all)
  
  # Per-cell table (immune-of-interest only)
  dat <- data.frame(
    Cell.Type = cd[[feat_col]],
    clustered = !is.na(cd$Neighborhood)
  )
  dat <- dat[dat$Cell.Type %in% cell_types_of_interest, , drop = FALSE]
  
  # Counts per Cell.Type × clustered
  counts <- dat |>
    dplyr::count(Cell.Type, clustered, name = "n_cells") |>
    tidyr::complete(
      Cell.Type = cell_types_of_interest,
      clustered = c(FALSE, TRUE),
      fill = list(n_cells = 0)
    ) |>
    as.data.frame()
  
  # --- robust density calculation (avoid any masking) ---
  area_ws <- as.numeric(area_ws)
  stopifnot(length(area_ws) == 1, is.finite(area_ws), area_ws > 0)
  
  counts$n_cells <- as.numeric(counts$n_cells)               # ensure numeric
  counts$Density_per_mm2 <- counts$n_cells / area_ws         # do division explicitly
  counts$Metric <- ifelse(counts$clustered, "Clustered_density", "Dispersed_density")
  
  # Add identifiers
  counts$Sample  <- nm
  counts$Pair_ID <- sub("^SPE_(\\d+)[P][0-9A-Za-z]*no$", "\\1", nm)
  counts$Group   <- ifelse(grepl("P0", nm), "Control", "Case")
  
  # Final order & return
  counts[, c("Sample","Pair_ID","Group","Cell.Type","n_cells","Density_per_mm2","Metric")]
})



dens_df <- dplyr::bind_rows(dens_list)
stopifnot(nrow(dens_df) > 0)

dens_df <- dens_df %>%
  dplyr::mutate(Group = factor(Group, levels = c("Control","Case")))

# Also compute TOTAL density per Cell.Type (clustered + dispersed)
total_df <- dens_df %>%
  dplyr::group_by(Sample, Pair_ID, Group, Cell.Type) %>%
  dplyr::summarise(
    Density_per_mm2 = sum(Density_per_mm2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Metric = "Total_density",
                n_cells = NA_integer_) %>%
  dplyr::select(Sample, Pair_ID, Group, Cell.Type, n_cells, Density_per_mm2, Metric)

dens_all <- dplyr::bind_rows(dens_df, total_df)

# ---- 2) Keep complete 1:1 pairs ----
valid_pairs <- dens_all %>%
  dplyr::distinct(Pair_ID, Group) %>%
  dplyr::count(Pair_ID, Group) %>%
  tidyr::pivot_wider(names_from = Group, values_from = n, values_fill = 0) %>%
  dplyr::filter(Control == 1, Case == 1) %>%
  dplyr::pull(Pair_ID)

dens_all <- dens_all %>% dplyr::filter(Pair_ID %in% valid_pairs)

# ---- 3) Paired Wilcoxon per (Cell.Type, Metric) + BH correction ----
stats <- dens_all %>%
  dplyr::group_by(Cell.Type, Metric) %>%
  dplyr::summarise(
    p = {
      wide <- tidyr::pivot_wider(cur_data_all(),
                                 id_cols = Pair_ID,
                                 names_from = Group,
                                 values_from = Density_per_mm2)
      if (all(c("Case","Control") %in% names(wide)) && nrow(wide) > 0) {
        stats::wilcox.test(wide$Case, wide$Control, paired = TRUE)$p.value
      } else NA_real_
    },
    n_pairs = dplyr::n_distinct(Pair_ID),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj = p.adjust(p, method = "BH"),
    label = paste0("p=", signif(p, 3), ifelse(!is.na(p_adj) & p_adj < 0.05, "*", ""))
  )

# ---- 4) Group summaries for Excel ----
summ <- dens_all %>%
  dplyr::group_by(Cell.Type, Metric, Group) %>%
  dplyr::summarise(
    mean_density   = mean(Density_per_mm2, na.rm = TRUE),
    median_density = median(Density_per_mm2, na.rm = TRUE),
    n_samples      = dplyr::n(),
    .groups = "drop"
  )

summary_table <- summ %>%
  dplyr::left_join(stats, by = c("Cell.Type","Metric")) %>%
  dplyr::arrange(Cell.Type, Metric, Group)

openxlsx::write.xlsx(summary_table, "clustered_vs_dispersed_density_stats.xlsx")

# ---- 5) Plot: boxplots with p-value labels (rows=Metric, cols=Cell.Type) ----
p_cdense <- ggplot(dens_all, aes(x = Group, y = Density_per_mm2, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5) +
  geom_text(
    data = stats,
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE, vjust = 1.3, size = 3
  ) +
  scale_fill_manual(values = c("Case" = "#E69F00", "Control" = "#56B4E9")) +
  facet_grid(Metric ~ Cell.Type, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = NULL, y = "Density (cells/mm²)") +
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

print(p_cdense)
ggsave("clustered_vs_dispersed_density_boxplots.png", p_cdense, width = 12, height = 8, dpi = 300)
ggsave("clustered_vs_dispersed_density_boxplots.pdf", p_cdense, width = 12, height = 8)

# ---- 6) Optional: Cytotoxic T-cell only plot ----
cyto_only <- dens_all %>% dplyr::filter(Cell.Type == "Cytotoxic T-cell")
stats_cyto <- stats %>% dplyr::filter(Cell.Type == "Cytotoxic T-cell")

p_cyto <- ggplot(cyto_only, aes(x = Group, y = Density_per_mm2, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5) +
  geom_text(
    data = stats_cyto,
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE, vjust = 1.3, size = 3
  ) +
  scale_fill_manual(values = c("Case" = "#E69F00", "Control" = "#56B4E9")) +
  facet_wrap(~ Metric, ncol = 3, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = NULL, y = "Cytotoxic T-cell density (cells/mm²)") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

print(p_cyto)
ggsave("cytotoxic_clustered_vs_dispersed_density.png", p_cyto, width = 8, height = 5, dpi = 300)



# ============================================================
# Supplementary Figure X: HDBSCAN cellular neighborhoods
# D-ESMEL study  —  version 4 (combined)
#
# Layout:
#   Row 1 (A–C): representative spatial examples
#   Row 2 (D–F): paired case–control results for 3 key metrics
#     D: total_clusters
#     E: frac_clustered_cells
#     F: frac_CD8_dominant_clusters
# ============================================================


# ---- 0) Settings ----
# output_dir is defined at the top of the script
fig_out_dir <- file.path(output_dir, "supp_figure_v4")
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

min_pts                <- 70
n_tiles                <- 6
feature_colname        <- "Cell.Type"
cell_types_of_interest <- c("T-cell", "B-cell lineage", "Cytotoxic T-cell", "Macrophage")

# ---- Colour scheme ----
cell_type_colours <- c(
  "Cytotoxic T-cell" = "#2a9d8f",
  "T-cell"           = "#4895ef",
  "B-cell lineage"   = "#e76f51",
  "Macrophage"       = "#9b5de5"
)
cell_type_labels <- c(
  "Cytotoxic T-cell" = "CD8\u207a T-cell",
  "T-cell"           = "CD3\u207a T-cell",
  "B-cell lineage"   = "B-cell",
  "Macrophage"       = "Macrophage"
)
group_colours <- c("Control" = "#457b9d", "Case" = "#e63946")

ALPHA_CLUSTERED <- 0.80
ALPHA_DISPERSED <- 0.12
SIZE_CLUSTERED  <- 0.55
SIZE_DISPERSED  <- 0.35

# ============================================================
# STEP 1: Load metrics
# ============================================================


metrics_files <- list.files(output_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(metrics_files) > 0)

metrics <- metrics_files %>%
  setNames(basename(.)) %>%
  purrr::map_dfr(
    ~{ obj <- readRDS(.x); if (is.data.frame(obj)) tibble::as_tibble(obj) else NULL },
    .id = "filename"
  ) %>%
  mutate(
    sample_id      = sub("\\.rds$", "", filename),
    SetID          = sub("P[01].*$", "", sample_id),
    Group          = factor(if_else(grepl("P0", sample_id), "Control", "Case"),
                            levels = c("Control", "Case")),
    total_clusters = as.integer(total_clusters)
  ) %>%
  select(-filename)

# Complete matched pairs only
pairs_df <- metrics %>%
  group_by(SetID) %>%
  filter(n_distinct(Group) == 2) %>%
  ungroup()

# ============================================================
# STEP 2: Select representative spatial samples
# ============================================================

immune_cap       <- quantile(metrics$total_immune_cells, 0.99, na.rm = TRUE)
metrics_filtered <- metrics %>%
  filter(total_immune_cells <= immune_cap, total_immune_cells > 0)

pick_rep <- function(df, target_clusters) {
  df %>%
    mutate(score = abs(total_clusters - target_clusters) +
             abs(total_immune_cells - median(df$total_immune_cells, na.rm = TRUE)) / 2000) %>%
    arrange(score) %>%
    slice_head(n = 1) %>%
    pull(sample_id)
}

q75 <- quantile(metrics_filtered$total_clusters, 0.75, na.rm = TRUE)
med <- median(metrics_filtered$total_clusters, na.rm = TRUE)

sample_high   <- pick_rep(filter(metrics_filtered, total_clusters >= q75), q75)
sample_medium <- pick_rep(filter(metrics_filtered,
                                 total_clusters >= floor(med),
                                 total_clusters < q75), med)

zero_candidates <- filter(metrics_filtered, total_clusters == 0)
if (nrow(zero_candidates) > 0) {
  sample_sparse <- zero_candidates %>%
    arrange(desc(total_immune_cells)) %>%
    slice_head(n = 1) %>%
    pull(sample_id)
  sparse_title <- "No clusters detected"
} else {
  sample_sparse <- pick_rep(filter(metrics_filtered, total_clusters <= 2), 0)
  sparse_title  <- "Low clustering"
}

tier_map <- list(
  list(id = sample_high,   title = "High clustering"),
  list(id = sample_medium, title = "Moderate clustering"),
  list(id = sample_sparse, title = sparse_title)
)

# ============================================================
# STEP 3: Extract spatial data
# ============================================================

extract_clusters_for_plot <- function(spe_object, cell_types_of_interest,
                                      min_pts = 70, feature_colname = "Cell.Type",
                                      n_tiles = 6) {
  cell_data <- as.data.frame(cbind(colData(spe_object), spatialCoords(spe_object)))
  all_immune <- cell_data %>%
    filter(.data[[feature_colname]] %in% cell_types_of_interest) %>%
    transmute(Cell.ID, Cell.Type = .data[[feature_colname]],
              x = Cell.X.Position, y = Cell.Y.Position,
              cluster_id = NA_character_)
  if (nrow(all_immune) == 0) return(all_immune)
  
  x_breaks        <- quantile(all_immune$x, probs = seq(0, 1, length.out = n_tiles + 1))
  all_immune$tile <- cut(all_immune$x, breaks = x_breaks,
                         include.lowest = TRUE, labels = FALSE)
  cluster_counter <- 0L
  
  for (tile_num in seq_len(n_tiles)) {
    idx   <- which(all_immune$tile == tile_num)
    cells <- all_immune[idx, ]
    if (nrow(cells) < min_pts) next
    pts    <- as.matrix(cells[, c("x", "y")])
    result <- try(dbscan::hdbscan(pts, minPts = min_pts), silent = TRUE)
    if (inherits(result, "try-error")) next
    shifted <- ifelse(result$cluster == 0L, NA_integer_,
                      result$cluster + cluster_counter)
    all_immune$cluster_id[idx] <- ifelse(is.na(shifted), NA_character_,
                                         paste0("Cluster_", shifted))
    mx <- suppressWarnings(max(shifted, na.rm = TRUE))
    if (is.finite(mx)) cluster_counter <- max(cluster_counter, mx)
  }
  all_immune
}

plot_data_list <- list()

for (tier in tier_map) {
  spe_name <- paste0("SPE_", tier$id)
  if (!exists(spe_name)) {
    # Uncomment to load from disk:
    # spe_path <- file.path("Z:/path/to/spe", paste0(spe_name, ".rds"))
    # if (file.exists(spe_path)) assign(spe_name, readRDS(spe_path)) else {
    #   warning("SPE not found: ", spe_name); next }
    warning("SPE not in environment: ", spe_name); next
  }
  
  spe <- get(spe_name)
  cat("Extracting:", tier$id, "\n")
  df <- extract_clusters_for_plot(spe, cell_types_of_interest,
                                  min_pts, feature_colname, n_tiles)
  rm(spe); gc()
  
  df <- mutate(df,
               x = x - min(x, na.rm = TRUE),
               y = y - min(y, na.rm = TRUE))
  
  n_cl      <- length(na.omit(unique(df$cluster_id)))
  pct_clust <- round(100 * mean(!is.na(df$cluster_id)), 0)
  
  df$panel_title    <- tier$title
  df$panel_subtitle <- paste0(
    n_cl, " cluster", if (n_cl != 1) "s" else "",
    "  \u2022  ", format(nrow(df), big.mark = ","), " immune cells",
    "  \u2022  ", pct_clust, "% clustered"
  )
  plot_data_list[[tier$title]] <- df
}

stopifnot(length(plot_data_list) > 0)

# ============================================================
# STEP 4: Spatial panel builder
# ============================================================

make_scalebar <- function(df, bar_um = 1000, gap_frac = 0.04) {
  xrange  <- diff(range(df$x, na.rm = TRUE))
  yrange  <- diff(range(df$y, na.rm = TRUE))
  x_end   <- max(df$x, na.rm = TRUE) - xrange * gap_frac
  x_start <- x_end - bar_um
  y_bar   <- min(df$y, na.rm = TRUE) + yrange * gap_frac * 1.5
  bar_label <- if (bar_um >= 1000) paste0(bar_um / 1000, " mm") else paste0(bar_um, " \u00b5m")
  list(
    segment = data.frame(x = x_start, xend = x_end, y = y_bar, yend = y_bar),
    label   = data.frame(x = (x_start + x_end) / 2,
                         y = y_bar + yrange * 0.028, label = bar_label)
  )
}

make_spatial_panel <- function(df) {
  df_dispersed <- filter(df, is.na(cluster_id))
  df_clustered <- filter(df, !is.na(cluster_id))
  bar_um <- if (diff(range(df$x, na.rm = TRUE)) >= 6000) 2000 else 1000
  sb     <- make_scalebar(df, bar_um = bar_um)
  
  ggplot() +
    geom_point(data = df_dispersed,
               aes(x = x, y = y, colour = Cell.Type),
               size = SIZE_DISPERSED, alpha = ALPHA_DISPERSED, shape = 16) +
    geom_point(data = df_clustered,
               aes(x = x, y = y, colour = Cell.Type),
               size = SIZE_CLUSTERED, alpha = ALPHA_CLUSTERED, shape = 16) +
    geom_segment(data = sb$segment,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 colour = "grey20", linewidth = 0.75,
                 lineend = "square", inherit.aes = FALSE) +
    geom_text(data = sb$label,
              aes(x = x, y = y, label = label),
              size = 2.4, colour = "grey20",
              hjust = 0.5, vjust = 0, inherit.aes = FALSE) +
    scale_colour_manual(
      values = cell_type_colours, labels = cell_type_labels,
      name = "Cell type",
      guide = guide_legend(
        override.aes = list(size = 2.8, alpha = 1, shape = 16), nrow = 1)
    ) +
    coord_fixed(clip = "off") +
    labs(title = unique(df$panel_title), subtitle = unique(df$panel_subtitle)) +
    theme_classic(base_size = 8) +
    theme(
      plot.title    = element_text(size = 8, face = "bold",
                                   colour = "grey10", margin = margin(b = 1)),
      plot.subtitle = element_text(size = 6.5, colour = "grey40",
                                   margin = margin(b = 4)),
      axis.title = element_blank(), axis.text  = element_blank(),
      axis.ticks = element_blank(), axis.line  = element_blank(),
      panel.border    = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
      legend.position = "none",
      plot.margin     = margin(4, 8, 4, 8)
    )
}

spatial_panels <- map(plot_data_list, make_spatial_panel)

# ============================================================
# STEP 5: Paired line plot builder
# ============================================================

# Load adjusted p-values from your existing Excel output
stats_file <- file.path(output_dir, "metrics_summary_plus.xlsx")
if (file.exists(stats_file)) {
  stats_results <- readxl::read_excel(stats_file, sheet = "stats_results")
} else {
  stop("Cannot find metrics_summary_plus.xlsx — run the main metrics script first.")
}

format_p <- function(p) {
  if (is.na(p))  return("p = NA")
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01)  return(paste0("p = ", sprintf("%.3f", p)))
  paste0("p = ", sprintf("%.2f", p))
}

get_p_adj <- function(metric_name) {
  row <- stats_results %>% filter(metric == metric_name)
  if (nrow(row) == 0 || is.null(row$p_value_adj)) return(NA_real_)
  row$p_value_adj[1]
}

make_paired_plot <- function(metric_col, y_label, title_label,
                             y_percent = FALSE,
                             log_scale  = FALSE) {
  df_plot <- pairs_df %>%
    transmute(SetID, Group,
              value = coerce_to_numeric(.data[[metric_col]])) %>%
    filter(!is.na(value)) %>%
    group_by(SetID) %>%
    filter(n_distinct(Group) == 2) %>%
    ungroup()
  
  p_adj   <- get_p_adj(metric_col)
  p_label <- paste0("adj. ", format_p(p_adj))
  
  # p-value bracket position: based on untransformed max
  y_max   <- max(df_plot$value, na.rm = TRUE)
  y_range <- diff(range(df_plot$value, na.rm = TRUE))
  
  # p-value bracket: annotate() uses original data scale even with trans,
  # so for log1p we supply raw values above the top axis break (200)
  if (log_scale) {
    y_ann      <- 260   # raw value; renders above 200 break on log1p axis
    y_ann_text <- 320
  } else {
    y_ann      <- y_max + y_range * 0.10
    y_ann_text <- y_max + y_range * 0.15
  }
  
  p <- ggplot(df_plot, aes(x = Group, y = value)) +
    # Paired connecting lines — slightly less dense
    geom_line(aes(group = SetID),
              colour = "grey75", linewidth = 0.22, alpha = 0.40) +
    # Individual points
    geom_point(aes(colour = Group),
               size = 0.9, alpha = 0.45, shape = 16) +
    # Median crossbar — black so it always reads clearly
    stat_summary(fun = median, geom = "crossbar",
                 colour = "grey15",
                 width = 0.38, linewidth = 0.9, fatten = 0) +
    # p-value bracket
    annotate("segment",
             x = 1, xend = 2, y = y_ann, yend = y_ann,
             colour = "grey50", linewidth = 0.35) +
    annotate("text",
             x = 1.5, y = y_ann_text,
             label = p_label, size = 2.3, colour = "grey30", hjust = 0.5) +
    scale_colour_manual(values = group_colours, guide = "none") +
    labs(title = title_label, x = NULL, y = y_label) +
    theme_classic(base_size = 8) +
    theme(
      plot.title   = element_text(size = 8, face = "bold",
                                  colour = "grey10", margin = margin(b = 2)),
      axis.title.y = element_text(size = 7, colour = "grey30"),
      axis.text    = element_text(size = 7, colour = "grey30"),
      axis.line    = element_line(colour = "grey70", linewidth = 0.3),
      axis.ticks   = element_line(colour = "grey70", linewidth = 0.3),
      panel.border = element_blank(),
      plot.margin  = margin(4, 14, 4, 8)
    )
  
  # Apply the appropriate y scale
  if (log_scale) {
    p <- p +
      scale_y_continuous(
        trans   = "log1p",
        breaks  = c(0, 1, 5, 10, 25, 50, 100, 200),
        expand  = expansion(mult = c(0.05, 0.25))
      )
  } else if (y_percent) {
    p <- p +
      scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        expand = expansion(mult = c(0.05, 0.25))
      )
  } else {
    p <- p +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
  }
  
  p
}

# Panel D: log1p scale to handle the high-clustering outlier
panel_D <- make_paired_plot(
  "total_clusters",
  y_label     = "Number of clusters (log scale)",
  title_label = "Total immune clusters",
  log_scale   = TRUE
)

panel_E <- make_paired_plot(
  "frac_clustered_cells",
  y_label     = "Fraction of immune cells",
  title_label = "Clustered immune cells",
  y_percent   = TRUE
)

# Panel F: threshold noted in y-axis label
panel_F <- make_paired_plot(
  "frac_CD8_dominant_clusters",
  y_label     = "Fraction of clusters (>60% CD8\u207a)",
  title_label = "CD8\u207a-dominant clusters",
  y_percent   = TRUE
)

results_row <- list(panel_D, panel_E, panel_F)

# ============================================================
# STEP 6: Assemble composite figure
# ============================================================

# Row 1: spatial (A–C) with shared cell type legend at bottom
row1 <- wrap_plots(spatial_panels, nrow = 1, guides = "collect") &
  theme(
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.text      = element_text(size = 7.5),
    legend.title     = element_text(size = 8, face = "bold"),
    legend.key.size  = unit(0.45, "cm"),
    legend.spacing.x = unit(0.4, "cm")
  )

# Row 2: paired plots (D–F)
row2 <- wrap_plots(results_row, nrow = 1)

# Stack; spatial row slightly taller
final_fig <- (row1 / row2) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag    = element_text(size = 9, face = "bold", colour = "grey10"),
    plot.margin = margin(4, 4, 4, 4)
  )

# ============================================================
# STEP 7: Save
# ============================================================

ggsave(
  filename = file.path(fig_out_dir, "Supplementary_Figure_HDBSCAN_v4.pdf"),
  plot     = final_fig,
  width    = 10.5, height = 8.0,
  device   = cairo_pdf
)

ggsave(
  filename = file.path(fig_out_dir, "Supplementary_Figure_HDBSCAN_v4.png"),
  plot     = final_fig,
  width    = 10.5, height = 8.0,
  dpi      = 300
)

cat("\nFigure saved to:\n  ", fig_out_dir, "\n")

