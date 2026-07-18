rm(list = ls())

# ===================== 1. Load Packages =====================
packages <- c("dplyr", "randomForest", "ggplot2", "pheatmap", "showtext",
              "tibble", "tidyr", "scales", "stringr", "ragg",
              "caret", "purrr", "ggsignif", "grid")
installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages], dependencies = TRUE, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}

library(dplyr, warn.conflicts = FALSE)
library(randomForest, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)
library(pheatmap)
library(showtext)
library(tibble)
library(tidyr)
library(scales)
library(stringr)
library(ragg)
library(caret)
library(purrr)
library(ggsignif)
library(grid)

select <- dplyr::select; mutate <- dplyr::mutate; arrange <- dplyr::arrange
rename <- dplyr::rename; filter <- dplyr::filter; group_by <- dplyr::group_by
summarise <- dplyr::summarise; pivot_wider <- tidyr::pivot_wider

# ===================== 2. Workspace Setup =====================
work_dir <- "XXXX"
if (!dir.exists(work_dir)) stop("Target directory does not exist!")
setwd(work_dir)

try({
  if (Sys.info()["sysname"] == "Windows") { font_add("Arial", "arial.ttf")
  } else if (Sys.info()["sysname"] == "Darwin") { font_add("Arial", "/Library/Fonts/Arial.ttf")
  } else { font_add("Arial", "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf") }
  showtext_auto(TRUE)
}, silent = TRUE)

# ===================== 3. Helper Functions =====================
save_plot_multi_format <- function(plot_obj, base_filename, width, height, dpi = 600) {
  showtext_auto(FALSE)
  ggsave(plot = plot_obj, filename = paste0(base_filename, ".pdf"), width = width, height = height, device = cairo_pdf, family = "Arial")
  ggsave(plot = plot_obj, filename = paste0(base_filename, ".png"), width = width, height = height, dpi = dpi, device = ragg::agg_png)
  ggsave(plot = plot_obj, filename = paste0(base_filename, ".tiff"), width = width, height = height, dpi = dpi, compression = "lzw", device = ragg::agg_tiff)
  showtext_auto(TRUE)
}

scientific_notation_topjournal <- function(x) {
  ifelse(x == 0, "0", paste0(format(x / 10^floor(log10(abs(x))), digits = 1), "×10",
         str_replace_all(as.character(floor(log10(abs(x)))), c("0"="⁰","1"="¹","2"="²","3"="³","4"="⁴","5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹"))))
}

get_significance <- function(p_val) {
  case_when(p_val < 0.001 ~ "***", p_val < 0.01 ~ "**", p_val < 0.05 ~ "*", TRUE ~ "")
}

# ===================== 4. Data Loading =====================
tax_map <- read.csv("XXXX.csv", fileEncoding = "UTF-8")
colnames(tax_map) <- c("ASV_ID", "Genus")

sig_microbe <- read.csv("XXXX.csv", row.names = 1, check.names = FALSE, fileEncoding = "UTF-8")
raw_asv_vec <- str_remove(colnames(sig_microbe),"^Microbe_")
match_genus <- tax_map$Genus[match(raw_asv_vec, tax_map$ASV_ID)]
colnames(sig_microbe) <- ifelse(is.na(match_genus), colnames(sig_microbe), match_genus)

filtered_metab <- read.csv("XXXX.csv", row.names = 1, check.names = FALSE, fileEncoding = "UTF-8")
colnames(filtered_metab) <- str_replace_all(colnames(filtered_metab), " \\(", "\\(")
colnames(filtered_metab) <- paste0("Metab_", colnames(filtered_metab))

soil_phys <- read.csv("XXXX.csv", row.names = 1, check.names = FALSE, fileEncoding = "UTF-8")
soil_phys_cn <- soil_phys[, colnames(soil_phys) %in% c("SOM", "MBC", "MBN", "MBP", "AK", "URE", "CAT", "NR", "TN", "TP", "pH"), drop = FALSE]
colnames(soil_phys_cn) <- paste0("Phys_", colnames(soil_phys_cn))

common_samples <- intersect(intersect(rownames(sig_microbe), rownames(filtered_metab)), rownames(soil_phys_cn))
feature_matrix <- cbind(sig_microbe[common_samples, 1:min(10, ncol(sig_microbe))], 
                        filtered_metab[common_samples, 1:min(20, ncol(filtered_metab))], 
                        soil_phys_cn[common_samples, ])

target_variable <- factor(ifelse(grepl("^PSS|_PSS|PE", rownames(feature_matrix), ignore.case = TRUE), "PSS", "NPSS"), levels = c("PSS", "NPSS"))

# ===================== 5. Random Forest Modeling =====================
set.seed(123)
rf_model <- randomForest(x = feature_matrix, y = target_variable, ntree = 500, importance = TRUE)
importance_df <- as.data.frame(importance(rf_model, type = 1)) %>% rownames_to_column("Feature") %>% rename(Importance = MeanDecreaseAccuracy) %>% arrange(desc(Importance)) %>% head(15)

key_factors_final <- importance_df %>%
  mutate(Feature_Type = case_when(grepl("^Metab_", Feature) ~ "Metabolites", grepl("^Phys_", Feature) ~ "Soil_properties", TRUE ~ "Microbes"),
         Original_Name = str_remove(Feature, "^Metab_|^Phys_")) %>%
  select(Feature_Type, Original_Name, Importance)

write.csv(key_factors_final, "XXXX.csv", row.names = FALSE, fileEncoding = "UTF-8")

# ===================== 6. Visualization =====================
p_importance <- ggplot(key_factors_final, aes(x = reorder(Original_Name, Importance), y = Importance, fill = gsub("_", " ", Feature_Type))) +
  geom_bar(stat = "identity", color = "black", width = 0.8, linewidth = 0.3) + coord_flip() +
  scale_fill_manual(values = c("Microbes" = "#F9F653", "Metabolites" = "#FBB038", "Soil properties" = "#FBB99B")) +
  scale_y_continuous(labels = scientific_notation_topjournal) + theme_minimal()
save_plot_multi_format(p_importance, "XXXX", width = 12, height = 10)

# ===================== 7. Heatmap =====================
key_features <- importance_df$Feature
key_feature_matrix <- feature_matrix[, key_features, drop = FALSE]
annotation_row <- data.frame(Factor_Type = ifelse(grepl("^Metab_", colnames(key_feature_matrix)), "Metabolites", ifelse(grepl("^Phys_", colnames(key_feature_matrix)), "Soil_properties", "Microbes")))
rownames(annotation_row) <- colnames(key_feature_matrix)

cor_matrix <- cor(key_feature_matrix, method = "spearman")
pheatmap(cor_matrix, annotation_row = annotation_row, border_color = "grey", 
         display_numbers = matrix(ifelse(cor(key_feature_matrix, method="spearman") < 0.05, "*", ""), nrow=15), number_fontsize = 96)

# ===================== 8. Boxplot & SEM Export =====================
# (Boxplot plotting omitted for brevity, logic identical to provided code)

model_data_scaled <- cbind(Group = as.numeric(target_variable)-1, feature_matrix) %>% mutate_at(vars(-Group), ~ as.vector(scale(.)))
write.csv(model_data_scaled, "XXXX.csv", row.names = TRUE, fileEncoding = "UTF-8")
