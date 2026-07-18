rm(list=ls())

# ===================== 1. Setup Environment =====================
setwd("XXXX")
options(repos = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

required_pkgs <- c("vegan", "ggplot2", "dplyr", "ggrepel", "caret")
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

theme_imeta <- function() {
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 9),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    plot.margin = margin(10, 10, 10, 10, "mm")
  )
}

# ===================== 2. Load Data Functions =====================
read_microbial_data <- function(file_path, min_samples = 3) {
  df <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  numeric_cols <- sapply(df, is.numeric)
  df <- df[, numeric_cols, drop = FALSE]
  
  present <- rowSums(df > 0) >= min_samples
  df <- df[present, , drop = FALSE]
  
  df_t <- as.data.frame(t(df))
  rownames(df_t) <- make.names(rownames(df_t), unique = TRUE)
  if (any(is.na(df_t))) df_t[is.na(df_t)] <- 0
  return(df_t)
}

read_soil_data <- function(file_path) {
  df <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  numeric_cols <- sapply(df, is.numeric)
  df <- df[, numeric_cols, drop = FALSE]
  rownames(df) <- make.names(rownames(df), unique = TRUE)
  
  if (any(is.na(df))) {
    df <- as.data.frame(lapply(df, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)))
  }
  return(df)
}

# ===================== 3. Procrustes Core Function =====================
perform_procrustes <- function(community_data, soil_data, title, 
                               method = "bray", transform = "hellinger", 
                               permutations = 999, save_plot = TRUE) {
  
  common_samples <- intersect(rownames(community_data), rownames(soil_data))
  community_filtered <- community_data[common_samples, , drop = FALSE]
  soil_filtered <- soil_data[common_samples, , drop = FALSE]
  
  if (ncol(soil_filtered) > 1) {
    cor_matrix <- cor(soil_filtered)
    high_cor <- caret::findCorrelation(cor_matrix, cutoff = 0.9)
    if (length(high_cor) > 0) soil_filtered <- soil_filtered[, -high_cor, drop = FALSE]
  }
  
  community_transformed <- if (transform == "hellinger") decostand(community_filtered, method = "hellinger") else community_filtered
  community_dist <- vegdist(community_transformed, method = "bray")
  soil_dist <- dist(soil_filtered, method = "euclidean")
  
  k <- min(nrow(community_filtered) - 1, 2)
  community_pcoa <- cmdscale(community_dist, k = k, eig = TRUE)
  soil_pcoa <- cmdscale(soil_dist, k = k, eig = TRUE)
  
  procrustes_result <- procrustes(community_pcoa$points, soil_pcoa$points)
  perm_test <- protest(community_pcoa$points, soil_pcoa$points, permutations = permutations)
  
  rotated_community <- as.data.frame(community_pcoa$points)
  rotated_soil <- as.data.frame(soil_pcoa$points)
  colnames(rotated_community) <- colnames(rotated_soil) <- c("PC1", "PC2")
  
  rotated_community$Sample <- rotated_soil$Sample <- rownames(rotated_community)
  rotated_community$Group <- rotated_soil$Group <- case_when(
    grepl("^PSS", rotated_community$Sample) ~ "PSS",
    grepl("^NPSS", rotated_community$Sample) ~ "NPSS",
    TRUE ~ "Unknown"
  )
  
  community_eig <- round(community_pcoa$eig[1:2]/sum(community_pcoa$eig)*100, 1)
  group_colors <- c("PSS" = "#d00000", "NPSS" = "#2f7eb0", "Unknown" = "gray50")
  
  p <- ggplot() +
    geom_point(data = rotated_community, aes(PC1, PC2, color = Group), size = 3, shape = 16, alpha = 0.8) +
    geom_point(data = rotated_soil, aes(PC1, PC2, color = Group), size = 3, shape = 17, alpha = 0.8) +
    geom_segment(data = data.frame(x=rotated_community$PC1, y=rotated_community$PC2, xend=rotated_soil$PC1, yend=rotated_soil$PC2, Sample=rotated_community$Sample),
                 aes(x, y, xend=xend, yend=yend), linetype="dashed", color="gray50", alpha=0.6) +
    geom_text_repel(data = rotated_community, aes(PC1, PC2, label=Sample, color=Group), size=3) +
    labs(title = title, subtitle = paste0("M² = ", round(procrustes_result$ss, 3), ", p = ", round(perm_test$signif, 3)),
         x = paste0("PC1 (", community_eig[1], "%)"), y = paste0("PC2 (", community_eig[2], "%)")) +
    scale_color_manual(values = group_colors) + theme_imeta() + theme(legend.position = "bottom")
  
  if (save_plot) {
    ggsave(paste0(title, "_procrustes_plot.png"), plot = p, width = 180, height = 150, units = "mm", dpi = 300)
    ggsave(paste0(title, "_procrustes_plot.pdf"), plot = p, width = 180, height = 150, units = "mm", device = "pdf")
  }
  
  return(list(plot = p, result = procrustes_result, perm_test = perm_test))
}

process_data <- function(microbial_file, soil_data, title) {
  microbial_data <- read_microbial_data(microbial_file)
  common <- intersect(rownames(microbial_data), rownames(soil_data))
  perform_procrustes(microbial_data[common, , drop=FALSE], soil_data[common, , drop=FALSE], title)
}

# ===================== 4. Execution =====================
soil <- read_soil_data("XXXX.csv")
archaea_proc <- process_data("XXXX.csv", soil, "Archaea")
bacteria_proc <- process_data("XXXX.csv", soil, "Bacteria")
protists_proc <- process_data("XXXX.csv", soil, "Protists")
fungi_proc <- process_data("XXXX.csv", soil, "Fungi")

save(archaea_proc, bacteria_proc, protists_proc, fungi_proc, file = "XXXX.RData")
