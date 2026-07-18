rm(list = ls())

# ===================== 1. Setup Environment =====================
packages <- c("psych", "dplyr", "tidyr", "tibble", "readr", "pheatmap", "RColorBrewer", "stringi")
installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages], dependencies = TRUE, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}
invisible(lapply(packages, library, character.only = TRUE))

work_dir <- "XXXX"
setwd(work_dir)

if (!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
}

clean_text <- function(text) {
  text <- stri_trans_general(text, "latin-ascii")
  text <- gsub("[^A-Za-z0-9 \\-_\\(\\)\\.]", " ", text)
  text <- gsub("\\s+", " ", text)
  text <- trimws(text)
  ifelse(text == "", "unknown", text)
}

# ===================== 2. Load and Process Microbial Data =====================
microbe_files <- list(
  Bacteria = "XXXX.csv",
  Archaea = "XXXX.csv",
  Protists = "XXXX.csv",
  Fungi = "XXXX.csv"
)

process_microbe_data <- function(file_path, microbe_type) {
  data <- read_csv(file_path, show_col_types = FALSE)
  colnames(data) <- sapply(colnames(data), clean_text)
  
  data_t <- as.data.frame(t(data))
  new_colnames <- sapply(as.character(data_t[1, ]), clean_text)
  data_t <- data_t[-1, , drop = FALSE]
  colnames(data_t) <- new_colnames
  
  raw_rownames <- make.unique(toupper(gsub("[^A-Za-z0-9]", "", rownames(data_t))), sep = "_")
  idx_empty <- which(raw_rownames == "")
  if(length(idx_empty) > 0) raw_rownames[idx_empty] <- paste0("Sample_", seq_along(idx_empty))
  rownames(data_t) <- raw_rownames
  
  data_t[] <- lapply(data_t, function(x) as.numeric(as.character(x)))
  colnames(data_t) <- paste0(microbe_type, "_", colnames(data_t))
  return(data_t)
}

microbe_list <- list()
for(type in names(microbe_files)){
  fp <- microbe_files[[type]]
  if(file.exists(fp)) microbe_list[[type]] <- process_microbe_data(fp, type)
}

merged_microbe <- bind_cols(microbe_list)

# ===================== 3. Load and Process Metabolite Data =====================
metab_matrix <- read.csv("XXXX.csv", row.names = 1, stringsAsFactors = FALSE)
colnames(metab_matrix) <- make.unique(sapply(colnames(metab_matrix), clean_text), sep = "_")
rownames(metab_matrix) <- make.unique(toupper(gsub("[^A-Za-z0-9]", "", rownames(metab_matrix))), sep = "_")

# ===================== 4. Match Samples and Filtering =====================
common_samples <- intersect(rownames(merged_microbe), rownames(metab_matrix))
merged_microbe <- merged_microbe[common_samples, , drop = FALSE]
metab_matrix <- metab_matrix[common_samples, , drop = FALSE]

mean_ab <- apply(merged_microbe, 2, mean, na.rm = TRUE)
prop_present <- apply(merged_microbe > 0, 2, function(x) sum(x, na.rm = TRUE)/length(x))
filtered_microbes <- merged_microbe[, which(mean_ab > 0.01 & prop_present > 0.5), drop = FALSE]

metab_stats <- as.data.frame(t(apply(metab_matrix, 2, function(col) {
  non_zero <- col[col > 0]
  if (length(non_zero) == 0) return(c(cv = 0, mean_ab = 0))
  return(c(cv = sd(non_zero) / mean(non_zero), mean_ab = mean(non_zero)))
})))

keep_metab <- which(colSums(metab_matrix > 0) >= 5 & metab_stats$cv >= 0.2 & metab_stats$mean_ab >= 1000)
filtered_metabs <- metab_matrix[, keep_metab, drop = FALSE]

if(ncol(filtered_metabs) > 15) {
  filtered_metabs_top15 <- filtered_metabs[, order(colMeans(filtered_metabs), decreasing = TRUE)[1:15], drop = FALSE]
} else {
  filtered_metabs_top15 <- filtered_metabs
}
write.csv(filtered_metabs_top15, file.path(work_dir, "XXXX.csv"), row.names = TRUE)

# ===================== 5. Correlation Analysis =====================
cor_result <- corr.test(x = filtered_microbes, y = filtered_metabs, method = "spearman", adjust = "fdr")
cor_matrix <- cor_result$r
p_matrix <- cor_result$p

sig_indices <- which(abs(cor_matrix) >= 0.6 & p_matrix <= 0.05, arr.ind = TRUE)
sig_cor_df <- data.frame(
  Microbe = rownames(cor_matrix)[sig_indices[, "row"]],
  Metabolite = colnames(cor_matrix)[sig_indices[, "col"]],
  Correlation = cor_matrix[sig_indices],
  P_value = p_matrix[sig_indices]
) %>% distinct(Microbe, Metabolite, .keep_all = TRUE)

# ===================== 6. Visualization and Export =====================
if (nrow(sig_cor_df) > 0) {
  sig_microbes <- unique(sig_cor_df$Microbe)
  sig_cor_submatrix <- cor_matrix[sig_microbes, unique(sig_cor_df$Metabolite), drop = FALSE]
  
  microbe_annotation <- data.frame(Microbe_Type = factor(gsub("_.+", "", rownames(sig_cor_submatrix))))
  rownames(microbe_annotation) <- rownames(sig_cor_submatrix)
  
  pheatmap(sig_cor_submatrix, color = colorRampPalette(c("#0318FB", "white", "#FC322E"))(200),
           breaks = seq(-1, 1, length.out = 200), annotation_row = microbe_annotation,
           fontsize = 6, angle_col = 45, filename = file.path(work_dir, "XXXX.pdf"))
  
  sig_cor_df <- sig_cor_df %>% mutate(Microbe_Type = gsub("_.+", "", Microbe))
  write.csv(sig_cor_df, file.path(work_dir, "XXXX.csv"), row.names = FALSE)
  write.csv(filtered_microbes[, sig_microbes, drop = FALSE], file.path(work_dir, "XXXX.csv"), row.names = TRUE)
  write.csv(filtered_metabs, file.path(work_dir, "XXXX.csv"), row.names = TRUE)
} else {
  write.csv(filtered_microbes, file.path(work_dir, "XXXX.csv"), row.names = TRUE)
  write.csv(filtered_metabs, file.path(work_dir, "XXXX.csv"), row.names = TRUE)
}
