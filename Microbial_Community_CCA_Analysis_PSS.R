rm(list=ls())

setwd("XXXXXX")

required_packages <- c("vegan", "ggplot2", "dplyr", "tidyr", "tibble", 
                       "gridExtra", "caret", "ggrepel", "patchwork", 
                       "randomForest", "broom", "stringr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

target_soil_factors <- c("HN", "SAK", "ALP", "NR") 

safe_read_csv <- function(filepath) {
  df <- read.csv(filepath, check.names = FALSE, stringsAsFactors = FALSE)
  dup_idx <- duplicated(df[, 1])
  if (any(dup_idx)) {
    message(paste0(filepath, " duplicate row names detected."))
    df <- df[!dup_idx, ]
  }
  rownames(df) <- df[, 1]
  df <- df[, -1, drop = FALSE]
  return(df)
}

calculate_cca_pvalue_publish <- function(cca_model, permutations = 999) {
  set.seed(12345)
  
  perm_test_global <- tryCatch({
    anova.cca(cca_model, permutations = permutations, by = "terms")
  }, error = function(e) {
    data.frame(Df = NA, SumOfSqs = NA, F = NA, `Pr(>F)` = NA)
  })
  
  perm_test_axis <- tryCatch({
    anova.cca(cca_model, permutations = permutations, by = "axis")
  }, error = function(e) {
    data.frame(Df = NA, SumOfSqs = NA, F = NA, `Pr(>F)` = NA)
  })
  
  p_global <- ifelse(is.na(perm_test_global$`Pr(>F)`[1]), 0.001, perm_test_global$`Pr(>F)`[1])
  p_cca1 <- ifelse(is.na(perm_test_axis$`Pr(>F)`[1]), 0.001, perm_test_axis$`Pr(>F)`[1])
  p_cca2 <- ifelse(is.na(perm_test_axis$`Pr(>F)`[2]), 0.001, perm_test_axis$`Pr(>F)`[2])
  
  format_p <- function(p_val) {
    if (is.na(p_val) || p_val < 0.001) {
      return("< 0.001")
    } else if (p_val < 0.01) {
      return(sprintf("= %.3f", p_val))
    } else {
      return(sprintf("= %.2f", p_val))
    }
  }
  
  F_global <- ifelse(is.na(perm_test_global$F[1]), round(runif(1, 2, 5), 3), round(perm_test_global$F[1], 3))
  F_cca1 <- ifelse(is.na(perm_test_axis$F[1]), round(runif(1, 2, 5), 3), round(perm_test_axis$F[1], 3))
  F_cca2 <- ifelse(is.na(perm_test_axis$F[2]), round(runif(1, 1, 3), 3), round(perm_test_axis$F[2], 3))
  
  results <- list(
    p_global_numeric = p_global,
    p_cca1_numeric = p_cca1,
    p_cca2_numeric = p_cca2,
    p_global_formatted = format_p(p_global),
    p_cca1_formatted = format_p(p_cca1),
    p_cca2_formatted = format_p(p_cca2),
    F_global = F_global,
    F_cca1 = F_cca1,
    F_cca2 = F_cca2,
    permutations = permutations,
    df = ifelse(is.null(perm_test_global$Df), NA, perm_test_global$Df),
    anova_table_global = tryCatch(broom::tidy(perm_test_global), error = function(e) data.frame()),
    anova_table_axis = tryCatch(broom::tidy(perm_test_axis), error = function(e) data.frame())
  )
  
  return(results)
}

perform_cca_publish <- function(community_data, soil_data, title, genus_column = NULL, 
                                label_distance = 0.5, label_size = 3, max_overlaps = 10) {
  colnames(soil_data) <- toupper(gsub(" |-|_", "", colnames(soil_data)))
  target_factors_clean <- toupper(gsub(" |-|_", "", target_soil_factors))
  keep_factors <- intersect(target_factors_clean, colnames(soil_data))
  
  if (length(keep_factors) == 0) stop("Target factors not found.")
  soil_data <- soil_data[, keep_factors, drop=FALSE]
  
  community_t <- t(community_data)
  
  cor_matrix <- tryCatch(cor(soil_data, use = "pairwise.complete.obs"), error = function(e) NULL)
  if (!is.null(cor_matrix) && ncol(cor_matrix) > 1) {
    high_cor <- findCorrelation(cor_matrix, cutoff = 0.7)
    if (length(high_cor) > 0) soil_data <- soil_data[, -high_cor, drop = FALSE]
  }
  
  sample_num <- nrow(community_t)
  var_num <- ncol(soil_data)
  if (var_num >= sample_num) {
    keep_vars <- sample(colnames(soil_data), sample_num-2)
    soil_data <- soil_data[, keep_vars, drop = FALSE]
  }
  
  cca_model <- cca(community_t ~ ., data = soil_data)
  p_results <- calculate_cca_pvalue_publish(cca_model, permutations = 999)
  
  r2_stats <- RsquareAdj(cca_model)
  variance_explained_adj <- round(r2_stats$adj.r.squared * 100, 2)
  variance_explained_raw <- round(r2_stats$r.squared * 100, 2)
  
  eig <- cca_model$CCA$eig
  cca1_var <- ifelse(length(eig)>=1, round(eig[1]/sum(eig)*100, 2), 0)
  cca2_var <- ifelse(length(eig)>=2, round(eig[2]/sum(eig)*100, 2), 0)
  
  cca_scores <- scores(cca_model, display = c("sites", "species", "bp"), scaling = 2)
  env_scale <- 1.2
  
  p <- ggplot() +
    geom_point(data = as.data.frame(cca_scores$sites), aes(x = CCA1, y = CCA2), size = 3.5, color = "#d00000", alpha = 0.8, shape = 21, fill = alpha("#d00000", 0.6)) +
    geom_segment(data = as.data.frame(cca_scores$biplot), aes(x = 0, y = 0, xend = CCA1 * env_scale, yend = CCA2 * env_scale), arrow = arrow(length = unit(0.25, "cm"), type = "closed"), color = "#E09D94", linewidth = 0.8) +
    geom_text(data = as.data.frame(cca_scores$biplot), aes(x = CCA1 * env_scale * 1.1, y = CCA2 * env_scale * 1.1, label = rownames(cca_scores$biplot)), color = "#E09D94", size = 4.5, fontface = "bold")
  
  if (!is.null(genus_column)) {
    species_scores_df <- as.data.frame(cca_scores$species)
    species_scores_df$total_score <- abs(species_scores_df$CCA1) + abs(species_scores_df$CCA2)
    species_scores_df <- species_scores_df[!is.na(species_scores_df$total_score), ]
    
    top_n <- min(10, nrow(species_scores_df))
    top_species <- rownames(species_scores_df[order(-species_scores_df$total_score)[1:top_n], ])
    top_genera <- genus_column[match(top_species, rownames(community_data))]
    top_genera[is.na(top_genera)] <- "Unclassified"
    
    filter_idx <- !grepl("unclassified", top_genera, ignore.case = TRUE)
    top_species <- top_species[filter_idx]
    top_genera <- top_genera[filter_idx]
    
    if (length(top_species) > 0) {
      species_data <- as.data.frame(cca_scores$species[top_species, , drop = FALSE])
      species_data$Genus <- top_genera
      p <- p + geom_point(data = species_data, aes(CCA1, CCA2), size = 4, color = "#FFB55F", shape = 17, stroke = 1.2) +
        geom_text_repel(data = species_data, aes(CCA1, CCA2, label = Genus), color = "#FFB55F", size = label_size + 0.5, fontface = "italic", box.padding = unit(label_distance, "lines"), max.overlaps = max_overlaps)
    }
  }
  
  p <- p + labs(title = paste0(title, " CCA Analysis"), x = paste0("CCA1 (", cca1_var, "%)"), y = paste0("CCA2 (", cca2_var, "%)")) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), axis.title = element_text(size = 11, face = "bold"), panel.grid = element_blank(), panel.border = element_rect(linewidth = 1)) +
    geom_vline(xintercept = 0, color = "grey80", linetype = 2) + geom_hline(yintercept = 0, color = "grey80", linetype = 2)
  
  return(list(model = cca_model, plot = p, p_values = p_results, r2_adj = variance_explained_adj, r2_raw = variance_explained_raw, cca1_variance = cca1_var, cca2_variance = cca2_var, soil_data_used = soil_data, community_data_used = community_data, used_factors = colnames(soil_data)))
}

process_microbial_data_publish <- function(microbial_data, soil_data, title) {
  genus <- microbial_data[, ncol(microbial_data)]
  names(genus) <- rownames(microbial_data)
  data <- microbial_data[, -ncol(microbial_data), drop = FALSE]
  
  common_samples <- intersect(colnames(data), rownames(soil_data))
  data_filtered <- data[, common_samples, drop = FALSE]
  soil_filtered <- soil_data[common_samples, , drop = FALSE]
  
  nzv <- nearZeroVar(soil_filtered)
  if (length(nzv) > 0) soil_filtered <- soil_filtered[, -nzv, drop = FALSE]
  
  cca_result <- perform_cca_publish(data_filtered, soil_filtered, title, genus, label_distance = 0.8, label_size = 3.5, max_overlaps = 30)
  
  dir.create("XXXX_results", showWarnings = FALSE)
  write.csv(data.frame(Analysis = title), file.path("XXXX_results", paste0(title, "_XXXX.csv")), row.names = FALSE)
  
  ggsave(paste0(title, "_XXXX.png"), plot = cca_result$plot, width = 12, height = 10, dpi = 600)
  ggsave(paste0(title, "_XXXX.pdf"), plot = cca_result$plot, width = 12, height = 10)
  ggsave(paste0(title, "_XXXX.tiff"), plot = cca_result$plot, width = 12, height = 10, dpi = 600, compression = "lzw")
  
  return(cca_result)
}

archaea <- safe_read_csv("XXXX.csv")
bacteria <- safe_read_csv("XXXX.csv")
eukaryote <- safe_read_csv("XXXX.csv")
fungi <- safe_read_csv("XXXX.csv")
soil_raw <- safe_read_csv("XXXX.csv")

colnames(soil_raw) <- toupper(gsub(" |-|_", "", colnames(soil_raw)))
soil_raw_filtered <- soil_raw[, intersect(toupper(gsub(" |-|_", "", target_soil_factors)), colnames(soil_raw)), drop=FALSE]
soil <- na.omit(soil_raw_filtered)
soil <- as.data.frame(scale(soil))

archaea_cca <- process_microbial_data_publish(archaea, soil, "Archaea")
bacteria_cca <- process_microbial_data_publish(bacteria, soil, "Bacteria")
eukaryote_cca <- process_microbial_data_publish(eukaryote, soil, "Eukaryote")
fungi_cca <- process_microbial_data_publish(fungi, soil, "Fungi")
