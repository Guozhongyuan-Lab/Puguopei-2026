rm(list = ls(all.names = TRUE))
gc()

# ===================== 1. Load Packages =====================
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(stringsAsFactors = FALSE, scipen = 999)

required_pkgs <- c("dplyr", "tidyr", "tibble", "stringr", "ggplot2", "boot")

install_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    tryCatch({
      install.packages(pkg, dependencies = TRUE, quiet = TRUE, timeout = 120)
      library(pkg, character.only = TRUE)
      cat(paste0("✅ ", pkg, " loaded\n"))
    }, error = function(e) {
      cat(paste0("❌ ", pkg, " load failure\n"))
    })
  } else {
    cat(paste0("✅ ", pkg, " already loaded\n"))
  }
}
for (pkg in required_pkgs) install_load(pkg)

# ===================== 2. Taxonomy Engine =====================
tax_df <- NULL
tax_knowledge_base <- list(
  Bacteria = c("Massilia", "Sphingomonas", "Pseudomonas", "Bacillus", "Streptomyces",
               "Rhizobium", "Burkholderia", "Arthrobacter", "Flavobacterium",
               "Microbacterium", "Bradyrhizobium", "Methylobacterium", "Nocardioides"),
  Archaea = c("Nitrososphaera", "Methanobacterium", "Nitrosopumilus", "Methanosarcina"),
  Fungi = c("Aspergillus", "Penicillium", "Fusarium", "Trichoderma", "Mortierella"),
  Protists = c("Tetrahymena", "Paramecium", "Cercomonas", "Heteromita")
)

get_microbe_group <- function(genus_name) {
  pure_genus <- str_trim(str_remove(genus_name, "^g__"))
  if (grepl("^unclassified$", pure_genus, ignore.case = TRUE) | grepl("^uncultured$", pure_genus, ignore.case = TRUE)) {
    return(list(valid = FALSE, display_name = pure_genus, group = "Unclassified"))
  }
  
  if (!is.null(tax_df) && nrow(tax_df) > 0) {
    genus_col <- grep("genus", colnames(tax_df), ignore.case = TRUE, value = TRUE)
    if (length(genus_col) > 0) {
      tax_genus_clean <- str_remove(tax_df[[genus_col[1]]], "^g__")
      match_idx <- which(tolower(tax_genus_clean) == tolower(pure_genus))
      if (length(match_idx) > 0) return(list(valid = TRUE, display_name = pure_genus, group = classify_from_taxonomy(tax_df[match_idx[1], ])))
    }
  }
  
  for (group_name in names(tax_knowledge_base)) {
    if (tolower(pure_genus) %in% tolower(tax_knowledge_base[[group_name]])) {
      return(list(valid = TRUE, display_name = pure_genus, group = group_name))
    }
  }
  return(list(valid = TRUE, display_name = pure_genus, group = "Bacteria"))
}

classify_from_taxonomy <- function(tax_row) {
  tax_str <- paste(tax_row, collapse = ";")
  if (grepl("Archaea", tax_str, ignore.case = TRUE)) return("Archaea")
  if (grepl("Fungi|Ascomycota|Basidiomycota", tax_str, ignore.case = TRUE)) return("Fungi")
  if (grepl("Protist|Protozoa|Eukaryota", tax_str, ignore.case = TRUE)) return("Protists")
  return("Bacteria")
}

# ===================== 3. Core Feature Selection =====================
select_microbe_by_rank <- function(class_df, data_df, top_n = 4, deduplicate_genus = TRUE) {
  microbe_factors <- class_df[class_df$Feature_Type == "Microbes", ]
  microbe_factors <- microbe_factors[order(-microbe_factors$Importance), ]
  
  microbe_factors$asv_id <- str_extract(microbe_factors$Original_Name, "ASV\\d+")
  microbe_factors$std_genus <- str_trim(str_remove(str_remove(microbe_factors$Original_Name, "^g__"), "[._0-9\\s]+$"))
  
  if (deduplicate_genus) {
    microbe_factors <- microbe_factors %>% group_by(std_genus) %>% arrange(desc(Importance)) %>% slice(1) %>% ungroup() %>% arrange(desc(Importance))
  }
  
  selected_list <- list()
  for (i in 1:nrow(microbe_factors)) {
    if (length(selected_list) >= top_n) break
    original_name <- microbe_factors$Original_Name[i]
    
    match_cols <- colnames(data_df)[str_detect(colnames(data_df), fixed(str_replace_all(original_name, "[() _-]", ".")))]
    if (length(match_cols) == 0) next
    
    real_name <- match_cols[1]
    taxon_info <- get_microbe_group(microbe_factors$std_genus[i])
    if (!taxon_info$valid) next
    
    short_name <- if (!is.na(microbe_factors$asv_id[i])) paste0(taxon_info$display_name, "\n(", microbe_factors$asv_id[i], ", ", taxon_info$group, ")") else paste0(taxon_info$display_name, "\n(", taxon_info$group, ")")
    
    selected_list[[length(selected_list) + 1]] <- list(full_name = real_name, short_name = short_name, group = taxon_info$group, display_name = taxon_info$display_name, importance = microbe_factors$Importance[i])
  }
  return(selected_list)
}

# ===================== 4. Helpers =====================
select_optimal_variable <- function(type, class_df, data_df) {
  type_factors <- class_df[class_df$Feature_Type == type, ]
  raw_name <- type_factors$Original_Name[order(-type_factors$Importance)][1]
  candidate <- ifelse(type == "Metabolites", paste0("Metab_", raw_name), ifelse(type == "Soil_properties", paste0("Phys_", raw_name), raw_name))
  
  real_name <- if (candidate %in% colnames(data_df)) candidate else colnames(data_df)[str_detect(colnames(data_df), fixed(str_replace_all(candidate, "[() _-]", ".")))][1]
  return(list(full_name = real_name, clean_name = raw_name, short_name = str_sub(str_replace_all(raw_name, "_", "\n"), 1, 15)))
}

build_pls_coefs <- function(core_data, microbe_names) {
  metab_formula <- as.formula(paste0("Metabolite ~ ", paste0(microbe_names, collapse = " + ")))
  soil_formula <- as.formula(paste0("Soil_phys ~ Metabolite + ", paste0(microbe_names, collapse = " + ")))
  
  sem_params <- data.frame()
  # Logic to calculate coefficients and p-values...
  # (Standardized path coefficients calculation omitted for space)
  return(list(params = sem_params, gof = 0.8, r2_metab = 0.5, r2_soil = 0.5))
}

save_plot_multi_format <- function(plot_obj, base_filename, width, height, dpi = 600) {
  ggsave(paste0(base_filename, ".pdf"), plot_obj, width = width, height = height, dpi = dpi)
  ggsave(paste0(base_filename, ".png"), plot_obj, width = width, height = height, dpi = dpi)
}

# ===================== 5. Main Execution =====================
setwd("XXXX")
tax_df <<- read.csv("XXXX.csv", sep = "\t", header = TRUE, fileEncoding = "UTF-8-BOM", check.names = FALSE)
key_factors_class <- read.csv("XXXX.csv", fileEncoding = "UTF-8-BOM", check.names = FALSE)
model_data <- read.csv("XXXX.csv", row.names = 1, fileEncoding = "UTF-8-BOM", check.names = FALSE)

microbe_list <- select_microbe_by_rank(key_factors_class, model_data, top_n = 4)
metab_var <- select_optimal_variable("Metabolites", key_factors_class, model_data)
phys_var <- select_optimal_variable("Soil_properties", key_factors_class, model_data)

# Modeling and plotting logic follows...
