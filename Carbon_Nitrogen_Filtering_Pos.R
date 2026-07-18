rm(list = ls())
gc()

# ===================== 1. Setup Environment =====================
if (!require("dplyr")) {
  install.packages("dplyr", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(dplyr)
}

work_path <- "XXXX"
setwd(work_path)

meta_pos <- tryCatch({
  read.csv("XXXX.csv", header = TRUE, stringsAsFactors = FALSE, 
           na.strings = c("", "NA", "NA."), fileEncoding = "UTF-8")
}, error = function(e) {
  stop("File reading failed: XXXX.csv")
})

# ===================== 2. Define Keywords =====================
carbon_keywords <- c("carbon", "sugar", "glucose", "fructose", "lactose", "starch", "cellulose", 
                     "citrate", "malate", "succinate", "acetate", "pyruvate", "lactate", 
                     "glycolate", "oxalate", "formate", "methane", "ethanol", "acetaldehyde",
                     "tca", "krebs", "glycolysis", "pentose", "xylose", "ribose", "galactose")

nitrogen_keywords <- c("nitrogen", "amino", "ammonia", "ammonium", "nitrate", "nitrite", 
                       "urea", "glutamate", "glutamine", "aspartate", "asparagine", 
                       "alanine", "arginine", "lysine", "methionine", "cysteine", "serine", 
                       "threonine", "tryptophan", "tyrosine", "valine", "leucine", "isoleucine", 
                       "proline", "histidine", "phenylalanine", "amide", "amine", "nitrile", 
                       "purine", "pyrimidine")

# ===================== 3. Filtering Process =====================
meta_pos_clean <- meta_pos %>%
  mutate(
    Name_lower = tolower(Name),
    Class_I_lower = tolower(replace(Class_I, is.na(Class_I), "")),
    Class_II_lower = tolower(replace(Class_II, is.na(Class_II), "")),
    Class_III_lower = tolower(replace(Class_III, is.na(Class_III), "")),
    Class_IV_lower = tolower(replace(Class_IV, is.na(Class_IV), ""))
  )

carbon_metabolites <- meta_pos_clean %>%
  filter(
    grepl(paste(carbon_keywords, collapse = "|"), Name_lower) |
    grepl(paste(carbon_keywords, collapse = "|"), Class_I_lower) |
    grepl(paste(carbon_keywords, collapse = "|"), Class_II_lower) |
    grepl(paste(carbon_keywords, collapse = "|"), Class_III_lower) |
    grepl(paste(carbon_keywords, collapse = "|"), Class_IV_lower)
  )

nitrogen_metabolites <- meta_pos_clean %>%
  filter(
    grepl(paste(nitrogen_keywords, collapse = "|"), Name_lower) |
    grepl(paste(nitrogen_keywords, collapse = "|"), Class_I_lower) |
    grepl(paste(nitrogen_keywords, collapse = "|"), Class_II_lower) |
    grepl(paste(nitrogen_keywords, collapse = "|"), Class_III_lower) |
    grepl(paste(nitrogen_keywords, collapse = "|"), Class_IV_lower)
  )

cn_pos_all <- bind_rows(carbon_metabolites, nitrogen_metabolites) %>%
  distinct(Name, .keep_all = TRUE)

# ===================== 4. Export =====================
cat("Carbon Metabolites Count:", nrow(carbon_metabolites), "\n")
cat("Nitrogen Metabolites Count:", nrow(nitrogen_metabolites), "\n")
cat("Total Unique C/N Metabolites:", nrow(cn_pos_all), "\n")

write.csv(cn_pos_all, "XXXX.csv", row.names = FALSE, fileEncoding = "UTF-8")
