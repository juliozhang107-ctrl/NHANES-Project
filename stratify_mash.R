# ==============================================================================
# NHANES MASH (FIB-4 > 2.67) Cohort Definition & Comorbidity Stratification
# ==============================================================================
# This script defines the MASH cohort (MASLD participants with high risk of 
# advanced liver fibrosis: FIB-4 > 2.67) from 1999 to 2023, and stratifies 
# them across comorbidities and lifestyle factors.
# ==============================================================================

# 1. Load Required Libraries
cat("Loading required libraries...\n")
library(nhanesA)
library(dplyr)
library(survey)
library(readr)

# Set survey option to adjust for strata with lonely PSUs
options(survey.lonely.psu = "adjust")

# Configurations
CACHE_DIR <- "nhanes_cache"

# 2. Define Mappings for CBC tables (containing platelet count LBXPLTSI)
cbc_cycles <- list(
  "1999-2000" = "LAB25",
  "2001-2002" = "L25_B",
  "2003-2004" = "L25_C",
  "2005-2006" = "CBC_D",
  "2007-2008" = "CBC_E",
  "2009-2010" = "CBC_F",
  "2011-2012" = "CBC_G",
  "2013-2014" = "CBC_H",
  "2015-2016" = "CBC_I",
  "2017-2020 (Pre-Pandemic)" = "P_CBC",
  "2021-2023" = "CBC_L"
)

# Biochemistry mappings (to load cached AST and ALT)
biopro_cycles <- list(
  "1999-2000" = "LAB18",
  "2001-2002" = "L40_B",
  "2003-2004" = "L40_C",
  "2005-2006" = "BIOPRO_D",
  "2007-2008" = "BIOPRO_E",
  "2009-2010" = "BIOPRO_F",
  "2011-2012" = "BIOPRO_G",
  "2013-2014" = "BIOPRO_H",
  "2015-2016" = "BIOPRO_I",
  "2017-2020 (Pre-Pandemic)" = "P_BIOPRO",
  "2021-2023" = "BIOPRO_L"
)

# Helper Function: Load and Cache NHANES Datasets
load_nhanes_cached <- function(table_name, cache_dir = CACHE_DIR) {
  cache_file <- file.path(cache_dir, paste0(table_name, ".rds"))
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  } else {
    cat(paste0("  Downloading table ", table_name, "...\n"))
    df <- tryCatch({
      nhanes(table_name)
    }, error = function(e) {
      cat(paste0("  Warning: failed to download ", table_name, ": ", e$message, "\n"))
      return(NULL)
    })
    if (!is.null(df)) {
      saveRDS(df, cache_file)
    }
    return(df)
  }
}

# 3. Load Base Comorbidities Cohort Dataset
cat("Loading comorbidity dataset...\n")
base_data <- read.csv("usfli_nhanes_comorbidities.csv")
cat("Starting with ", nrow(base_data), " MASLD cases (expected: 8673)\n")

# 4. Fetch and Process Platelets, AST, and ALT for each Cycle
lab_list <- list()

for (cyc_name in names(cbc_cycles)) {
  cat("\nProcessing cycle:", cyc_name, "...\n")
  cyc_base <- base_data %>% filter(Cycle == cyc_name)
  if (nrow(cyc_base) == 0) next
  
  # Load CBC table (for platelets) and Biochemistry (for AST/ALT)
  cbc_table_name <- cbc_cycles[[cyc_name]]
  biopro_table_name <- biopro_cycles[[cyc_name]]
  
  cbc <- load_nhanes_cached(cbc_table_name)
  biopro <- load_nhanes_cached(biopro_table_name)
  
  # Clean tables
  cbc_clean <- if (!is.null(cbc) && "LBXPLTSI" %in% colnames(cbc)) cbc %>% select(SEQN, Platelets = LBXPLTSI) else data.frame(SEQN=integer(), Platelets=numeric())
  
  biopro_clean <- if (!is.null(biopro) && "LBXSASSI" %in% colnames(biopro) && "LBXSATSI" %in% colnames(biopro)) {
    biopro %>% select(SEQN, AST = LBXSASSI, ALT = LBXSATSI)
  } else {
    data.frame(SEQN=integer(), AST=numeric(), ALT=numeric())
  }
  
  # Merge with base data
  cyc_merged <- cyc_base %>%
    left_join(cbc_clean, by = "SEQN") %>%
    left_join(biopro_clean, by = "SEQN")
  
  lab_list[[cyc_name]] <- cyc_merged
}

# Combine all cycles
masld_lab <- bind_rows(lab_list)

# 5. Calculate FIB-4 and Risk Strata
cat("\nCalculating FIB-4 index...\n")

# FIB-4 = (Age * AST) / (Platelets * sqrt(ALT))
masld_lab <- masld_lab %>%
  mutate(
    # Compute FIB-4
    FIB4 = (RIDAGEYR * AST) / (Platelets * sqrt(ALT)),
    
    # Classify Risk Strata
    FIB4_Group = case_when(
      FIB4 < 1.30 ~ "Low Risk (<1.30)",
      FIB4 >= 1.30 & FIB4 <= 2.67 ~ "Intermediate Risk (1.30-2.67)",
      FIB4 > 2.67 ~ "High Risk (>2.67)",
      TRUE ~ NA_character_
    )
  )

# Write patient-level dataset with FIB-4 details
write.csv(masld_lab, "usfli_nhanes_mash.csv", row.names = FALSE)
cat("Saved detailed patient dataset with FIB-4 to 'usfli_nhanes_mash.csv'\n")

# Exclude participants with missing FIB-4
masld_valid_fib4 <- masld_lab %>% filter(!is.na(FIB4))
cat("\nMASLD Cases with complete FIB-4 data: ", nrow(masld_valid_fib4), " / ", nrow(masld_lab), "\n")

# Set up complex survey design
nhanes_design <- svydesign(
  id = ~SDMVPSU, 
  strata = ~SDMVSTRA, 
  weights = ~FastingWeight, 
  data = masld_valid_fib4, 
  nest = TRUE
)

# 6. Estimate FIB-4 Strata Distribution
cat("\n============================================================\n")
cat("FIB-4 RISK STRATA DISTRIBUTION IN MASLD COHORT\n")
cat("============================================================\n")

# Raw counts
raw_fib4_table <- table(masld_valid_fib4$FIB4_Group, useNA = "no")
raw_fib4_total <- sum(raw_fib4_table)

# Factor levels coercion
masld_valid_fib4$FIB4_Group <- factor(masld_valid_fib4$FIB4_Group, levels = names(raw_fib4_table))
nhanes_design_updated <- svydesign(
  id = ~SDMVPSU, 
  strata = ~SDMVSTRA, 
  weights = ~FastingWeight, 
  data = masld_valid_fib4, 
  nest = TRUE
)

# Survey weighted prevalence
weighted_fib4 <- svymean(~FIB4_Group, nhanes_design_updated, na.rm = TRUE)
weighted_fib4_pct <- as.numeric(weighted_fib4) * 100
weighted_fib4_se <- as.numeric(SE(weighted_fib4)) * 100

fib4_summary <- data.frame(
  Stratum = names(raw_fib4_table),
  Raw_N = as.numeric(raw_fib4_table),
  Raw_Percent = (as.numeric(raw_fib4_table) / raw_fib4_total) * 100,
  Weighted_Percent = weighted_fib4_pct,
  Weighted_SE = weighted_fib4_se,
  stringsAsFactors = FALSE
)
print(fib4_summary)

# 7. Perform Comorbidity Stratification on the High-Risk (FIB-4 > 2.67) MASH Cohort
cat("\n============================================================\n")
cat("STRATIFICATION OF MASH COHORT (FIB-4 > 2.67)\n")
cat("============================================================\n")

# Subset design to High Risk group
design_mash <- subset(nhanes_design_updated, FIB4_Group == "High Risk (>2.67)")
cat("MASH Cohort Size (High Risk Raw N): ", sum(masld_valid_fib4$FIB4_Group == "High Risk (>2.67)"), "\n")

# Re-use our stratification helper function
get_stratification_stats <- function(design_obj, var_name, category_title) {
  f <- as.formula(paste0("~", var_name))
  df_data <- design_obj$variables
  
  # Filter out missing for raw counts
  raw_table <- table(df_data[[var_name]], useNA = "no")
  raw_total <- sum(raw_table)
  
  if (raw_total == 0) {
    return(data.frame(
      Category = category_title,
      Stratum = "No Cases",
      Raw_N = 0,
      Raw_Percent = 0,
      Weighted_Percent = NA,
      Weighted_SE = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  design_obj$variables[[var_name]] <- factor(df_data[[var_name]], levels = names(raw_table))
  
  weighted_prev <- svymean(f, design_obj, na.rm = TRUE)
  weighted_vals <- as.numeric(weighted_prev) * 100
  weighted_se <- as.numeric(SE(weighted_prev)) * 100
  
  strata_names <- names(raw_table)
  svy_names <- gsub(paste0("^", var_name), "", names(coef(weighted_prev)))
  
  result_rows <- list()
  for (i in seq_along(strata_names)) {
    sname <- strata_names[i]
    raw_n <- as.numeric(raw_table[sname])
    raw_pct <- (raw_n / raw_total) * 100
    
    svy_idx <- which(svy_names == sname)
    w_pct <- if (length(svy_idx) > 0) weighted_vals[svy_idx] else NA
    w_se <- if (length(svy_idx) > 0) weighted_se[svy_idx] else NA
    
    result_rows[[i]] <- data.frame(
      Category = category_title,
      Stratum = sname,
      Raw_N = raw_n,
      Raw_Percent = raw_pct,
      Weighted_Percent = w_pct,
      Weighted_SE = w_se,
      stringsAsFactors = FALSE
    )
  }
  
  # Add a row for missing values
  missing_n <- sum(is.na(df_data[[var_name]]))
  result_rows[[length(result_rows) + 1]] <- data.frame(
    Category = category_title,
    Stratum = "Missing / Unknown",
    Raw_N = missing_n,
    Raw_Percent = (missing_n / nrow(df_data)) * 100,
    Weighted_Percent = NA,
    Weighted_SE = NA,
    stringsAsFactors = FALSE
  )
  
  return(do.call(rbind, result_rows))
}

# Run comorbidity stratification on the MASH cohort
obesity_stats <- get_stratification_stats(design_mash, "Obesity_Cat", "Obesity/BMI")
diabetes_stats <- get_stratification_stats(design_mash, "Diabetes_Cat", "Type 2 DM")
hypertension_stats <- get_stratification_stats(design_mash, "Hypertension_Cat", "Hypertension")
dyslipidemia_stats <- get_stratification_stats(design_mash, "Dyslipidemia_Cat", "Dyslipidemia")
diet_stats <- get_stratification_stats(design_mash, "Diet_Cat", "Diet (Caloric Intake)")
pa_stats <- get_stratification_stats(design_mash, "PA_Cat", "Physical Activity")
smoking_stats <- get_stratification_stats(design_mash, "Smoking_Cat", "Smoking Status")

# Mortality linkage (1999-2018)
design_mortality <- subset(design_mash, Cycle != "2021-2023" & !is.na(MORTSTAT))
mortality_stats <- get_stratification_stats(design_mortality, "Mortality_Cat", "Mortality (1999-2018 cohorts)")

# Combine MASH results
mash_stratification <- bind_rows(
  obesity_stats,
  diabetes_stats,
  hypertension_stats,
  dyslipidemia_stats,
  diet_stats,
  pa_stats,
  smoking_stats,
  mortality_stats
)
print(mash_stratification)

# Save summary results
write.csv(mash_stratification, "mash_comorbidities_stratification.csv", row.names = FALSE)
cat("\nSaved MASH stratification summary to 'mash_comorbidities_stratification.csv'\n")
cat("MASH analysis completed successfully!\n")
