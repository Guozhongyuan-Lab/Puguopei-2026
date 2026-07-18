rm(list = ls())

# ===================== 1. Setup Environment =====================
required_pkgs <- c("vegan", "dplyr", "readr", "ggplot2")
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = FALSE, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

# ===================== 2. Core Parameters =====================
setwd("XXXX")

plastic_samples <- c("PSS1", "PSS2", "PSS3", "PSS4", "PSS5", "PSS6")
non_plastic_samples <- c("NPSS1", "NPSS2", "NPSS3", "NPSS4", "NPSS5", "NPSS6")
all_expected_samples <- c(plastic_samples, non_plastic_samples)
color_palette <- list(
  PSS = "#D00000",
  NPSS = "#2F7EB0"
)

# ===================== 3. Custom Function: Read Microbial Data =====================
read_csv_microbe <- function(file_path, microbe_type) {
  df <- read_csv(file_path, show_col_types = FALSE) %>% as.data.frame()
  if (ncol(df) < 2) stop(paste("File", file_path, "format error, at least 2 columns required"))
  
  asv_names <- df[, 1]
  sample_data <- df[, -1, drop = FALSE]
  colnames(sample_data) <- gsub("[^A-Za-z0-9_]", "", colnames(sample_data))
  
  t_data <- as.data.frame(t(sample_data))
  colnames(t_data) <- paste0(microbe_type, "_", asv_names)
  t_data$sample_id <- rownames(t_data)
  t_data <- t_data[t_data$sample_id %in% all_expected_samples, , drop = FALSE]
  return(t_data)
}

# ===================== 4. Load and Merge Microbial Data =====================
microbe_files <- list(
  Bacteria = "XXXX.csv",
  Archaea = "XXXX.csv",
  Eukaryote = "XXXX.csv",
  Fungus = "XXXX.csv"
)

microbe_dfs <- list()
for (type in names(microbe_files)) {
  df <- read_csv_microbe(microbe_files[[type]], type)
  microbe_dfs[[type]] <- df
}

microbe_step1 <- inner_join(microbe_dfs$Bacteria, microbe_dfs$Archaea, by = "sample_id")
microbe_step2 <- inner_join(microbe_dfs$Eukaryote, microbe_dfs$Fungus, by = "sample_id")
microbe_merged <- inner_join(microbe_step1, microbe_step2, by = "sample_id")

microbe_matrix <- microbe_merged[, !colnames(microbe_merged) %in% "sample_id", drop = FALSE]
microbe_matrix <- as.data.frame(apply(microbe_matrix, 2, as.numeric))
rownames(microbe_matrix) <- microbe_merged$sample_id
microbe_matrix <- microbe_matrix[, colMeans(microbe_matrix, na.rm = TRUE) > 0.0001, drop = FALSE]

# ===================== 5. Load and Process Metabolite Data =====================
metab_df <- read_csv("XXXX.csv", show_col_types = FALSE) %>% as.data.frame()
if ("...1" %in% colnames(metab_df)) {
  colnames(metab_df)[colnames(metab_df) == "...1"] <- "sample_id_raw"
} else if (!"sample_id_raw" %in% colnames(metab_df)) {
  metab_df$sample_id_raw <- rownames(metab_df)
}

rownames(metab_df) <- metab_df$sample_id_raw
metab_df <- metab_df[, !colnames(metab_df) %in% "sample_id_raw", drop = FALSE]
metab_df$sample_id <- gsub("[^A-Za-z0-9_]", "", rownames(metab_df))
metab_matrix <- metab_df[, !colnames(metab_df) %in% "sample_id", drop = FALSE]
metab_matrix <- as.data.frame(apply(metab_matrix, 2, as.numeric))
rownames(metab_matrix) <- metab_df$sample_id
metab_matrix <- metab_matrix[rownames(metab_matrix) %in% all_expected_samples, , drop = FALSE]

# ===================== 6. Sample Matching and Standardization =====================
common_samples <- intersect(rownames(microbe_matrix), rownames(metab_matrix))
microbe_filtered <- microbe_matrix[common_samples, , drop = FALSE]
metab_filtered <- metab_matrix[common_samples, , drop = FALSE]

plastic_idx <- common_samples %in% plastic_samples
non_plastic_idx <- common_samples %in% non_plastic_samples

microbe_hellinger <- decostand(microbe_filtered, method = "hellinger")
metab_scaled <- scale(metab_filtered)

# ===================== 7. Procrustes Analysis =====================
microbe_pca <- rda(microbe_hellinger)
microbe_scores <- scores(microbe_pca, display = "sites", scaling = 1)
metab_pca <- rda(metab_scaled)
metab_scores <- scores(metab_pca, display = "sites", scaling = 1)

pro_fit <- procrustes(X = microbe_scores[, 1:2], Y = metab_scores[, 1:2], symmetric = TRUE)
pro_test <- protest(X = microbe_scores[, 1:2], Y = metab_scores[, 1:2], permutations = 999)
r_value <- round(sqrt(pro_fit$ss / sum(pro_fit$Ytot)), 3)

metab_scores_rotated <- pro_fit$Yrot
rownames(metab_scores_rotated) <- rownames(microbe_scores)
sample_residuals <- sqrt(rowSums((microbe_scores[, 1:2] - metab_scores_rotated[, 1:2])^2))

# ===================== 8. Visualization =====================
pdf("XXXX.pdf", width = 8, height = 7)
par(mar = c(5, 5, 5, 3), cex.axis = 1.1, cex.lab = 1.3, cex.main = 1.4, font.lab = 2, font.main = 2, xpd = TRUE)

x_lim <- range(c(microbe_scores[, 1], metab_scores[, 1])) * 1.1
y_lim <- range(c(microbe_scores[, 2], metab_scores[, 2])) * 1.1

microbe_pca1_var <- round(eigenvals(microbe_pca)[1] / sum(eigenvals(microbe_pca)) * 100, 1)
microbe_pca2_var <- round(eigenvals(microbe_pca)[2] / sum(eigenvals(microbe_pca)) * 100, 1)

plot(x = x_lim, y = y_lim, type = "n",
     xlab = paste0("Microbial Community PCA1 (", microbe_pca1_var, "%)"),
     ylab = paste0("Microbial Community PCA2 (", microbe_pca2_var, "%)"),
     main = "Procrustes Analysis of Microbial Communities and C/N Metabolites",
     xaxt = "n", yaxt = "n")

x_ticks <- seq(round(x_lim[1], 1), round(x_lim[2], 1), by = round(diff(x_lim)/6, 1))
y_ticks <- seq(round(y_lim[1], 1), round(y_lim[2], 1), by = round(diff(y_lim)/6, 1))
axis(1, at = x_ticks, labels = x_ticks, cex.axis = 1.1)
axis(2, at = y_ticks, labels = y_ticks, cex.axis = 1.1)

points(microbe_scores[non_plastic_idx, 1], microbe_scores[non_plastic_idx, 2], pch = 16, col = color_palette$NPSS, cex = 1.8)
points(metab_scores[non_plastic_idx, 1], metab_scores[non_plastic_idx, 2], pch = 17, col = color_palette$NPSS, cex = 1.8)
points(microbe_scores[plastic_idx, 1], microbe_scores[plastic_idx, 2], pch = 16, col = color_palette$PSS, cex = 1.8)
points(metab_scores[plastic_idx, 1], metab_scores[plastic_idx, 2], pch = 17, col = color_palette$PSS, cex = 1.8)

for (s in common_samples) {
  lines(x = c(microbe_scores[s, 1], metab_scores[s, 1]),
        y = c(microbe_scores[s, 2], metab_scores[s, 2]),
        col = ifelse(s %in% plastic_samples, color_palette$PSS, color_palette$NPSS),
        lty = 2, lwd = 1.2)
}

legend("topright", legend = c("PSS-Microbiome", "PSS-Metabolites", "NPSS-Microbiome", "NPSS-Metabolites", "Matching Line"),
       pch = c(16, 17, 16, 17, NA), lty = c(NA, NA, NA, NA, 2),
       col = c(color_palette$PSS, color_palette$PSS, color_palette$NPSS, color_palette$NPSS, "gray50"),
       bty = "n", cex = 1.1, pt.cex = 1.3)

text(x = x_lim[1] + 0.1 * diff(x_lim), y = y_lim[2] - 0.1 * diff(y_lim), labels = paste0("Procrustes r = ", r_value), cex = 1.3, font = 2, adj = 0)
text(x = x_lim[1] + 0.1 * diff(x_lim), y = y_lim[2] - 0.15 * diff(y_lim), labels = paste0("P = ", pro_test$signif), cex = 1.3, font = 2, adj = 0)
dev.off()

# ===================== 9. Save Results =====================
result_table <- data.frame(
  Sample_ID = common_samples,
  Group = ifelse(plastic_idx, "Plastic (PSS)", "Non-plastic (NPSS)"),
  Microbe_PCA1 = round(microbe_scores[common_samples, 1], 4),
  Microbe_PCA2 = round(microbe_scores[common_samples, 2], 4),
  Metabolite_PCA1 = round(metab_scores[common_samples, 1], 4),
  Metabolite_PCA2 = round(metab_scores[common_samples, 2], 4),
  Residuals = round(sample_residuals[common_samples], 4)
)
write.csv(result_table, "XXXX.csv", row.names = FALSE)

# ===================== 10. Statistical Analysis =====================
residual_stats <- result_table %>% group_by(Group) %>%
  summarise(Mean = round(mean(Residuals), 4), SD = round(sd(Residuals), 4))
write.csv(residual_stats, "XXXX.csv", row.names = FALSE)

shapiro_test <- shapiro.test(result_table$Residuals)
if (shapiro_test$p.value > 0.05) {
  test_res <- t.test(Residuals ~ Group, data = result_table)
} else {
  test_res <- wilcox.test(Residuals ~ Group, data = result_table)
}
write.csv(as.data.frame(do.call(cbind, test_res[c("statistic", "p.value")])), "XXXX.csv", row.names = FALSE)

# ===================== 11. Statistical Report =====================
stat_report <- data.frame(
  Item = c("Microbial Variance", "Metabolite Variance", "Procrustes r", "P-value"),
  Value = c(
    paste0(round(sum(eigenvals(microbe_pca)[1:2]) / sum(eigenvals(microbe_pca)) * 100, 1), "%"),
    paste0(round(sum(eigenvals(metab_pca)[1:2]) / sum(eigenvals(metab_pca)) * 100, 1), "%"),
    as.character(r_value), as.character(pro_test$signif)
  )
)
write.csv(stat_report, "XXXX.csv", row.names = FALSE)
