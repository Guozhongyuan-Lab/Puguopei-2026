rm(list = ls())
gc()

# ===================== 1. Setup Environment =====================
pkg_list <- c("dplyr", "tidyr", "tibble", "stringr")
need_install <- pkg_list[!pkg_list %in% rownames(installed.packages())]
if (length(need_install) > 0) {
  install.packages(need_install, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}
invisible(lapply(pkg_list, library, character.only = TRUE))

work_path <- "XXXX"
setwd(work_path)

# ===================== 2. Load Metabolite Data =====================
neg_raw <- tryCatch({
  read.csv("XXXX.csv", header = TRUE, stringsAsFactors = FALSE, na.strings = c("", "NA", "NA."), fileEncoding = "UTF-8")
}, error = function(e) stop("Missing file: XXXX.csv"))

pos_raw <- tryCatch({
  read.csv("XXXX.csv", header = TRUE, stringsAsFactors = FALSE, na.strings = c("", "NA", "NA."), fileEncoding = "UTF-8")
}, error = function(e) stop("Missing file: XXXX.csv"))

abnormal_metabs <- c("Pseudo-Lactose", "Unknown", "Blank")

neg_clean <- neg_raw %>% filter(!Name %in% abnormal_metabs & !str_detect(Name, "Pseudo-|Unknown|Blank"))
pos_clean <- pos_raw %>% filter(!Name %in% abnormal_metabs & !str_detect(Name, "Pseudo-|Unknown|Blank"))

# ===================== 3. Merging and Filtering =====================
metab_merge <- bind_rows(neg_clean, pos_clean) %>% distinct(Name, .keep_all = TRUE)

all_col <- colnames(metab_merge)
sample_col_vec <- grep("PSS|NPSS", all_col, ignore.case = TRUE, value = TRUE)
if (length(sample_col_vec) == 0) stop("No PSS/NPSS sample columns detected!")

metab_sub <- metab_merge %>% select(Name, all_of(sample_col_vec))
metab_filter <- metab_sub %>% filter(if_any(all_of(sample_col_vec), ~ !is.na(.) & . > 0))

# ===================== 4. Long Format and CV Calculation =====================
metab_long <- metab_filter %>%
  pivot_longer(cols = all_of(sample_col_vec), names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Sample = str_extract(SampleID, "PSS\\d*|NPSS\\d*"))

metab_agg <- metab_long %>%
  group_by(Sample, Name) %>%
  summarise(
    Intensity_mean = mean(Intensity, na.rm = TRUE),
    sd_val = sd(Intensity, na.rm = TRUE),
    CV = ifelse(mean(Intensity, na.rm = TRUE) == 0, 999, sd_val / mean(Intensity, na.rm = TRUE)),
    .groups = "drop"
  )

metab_cv_filter <- metab_agg %>% filter(CV < 3 | is.na(CV)) %>% select(-CV, -sd_val)

# ===================== 5. Matrix Construction and QC =====================
metab_final <- metab_cv_filter %>%
  pivot_wider(names_from = Name, values_from = Intensity_mean) %>%
  column_to_rownames("Sample")

metab_final[is.na(metab_final)] <- apply(metab_final, 2, median, na.rm = TRUE)

write.csv(metab_final, "XXXX.csv", row.names = TRUE, fileEncoding = "UTF-8")

qc_summary <- data.frame(
  Neg_Metab_Count = nrow(neg_clean),
  Pos_Metab_Count = nrow(pos_clean),
  Merged_Metab_Count = nrow(metab_merge),
  Filtered_Metab_Count = ncol(metab_final),
  Sample_Count = nrow(metab_final)
)
write.csv(qc_summary, "XXXX.csv", row.names = FALSE, fileEncoding = "UTF-8")
