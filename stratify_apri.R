# ==============================================================================
# NHANES Liver Fibrosis (APRI > 0.5, > 1.0, and >= 2.0) Cohort Definition & Stratification
# ==============================================================================
# This script defines liver fibrosis and cirrhosis cohorts using the AST to 
# Platelet Ratio Index (APRI) from 1999 to 2023, and stratifies them across 
# comorbidities.
# ==============================================================================

# 1. Load Required Libraries
cat("Loading required libraries...\n")
library(dplyr)
library(survey)

# Set survey option to adjust for lonely PSUs
options(survey.lonely.psu = "adjust")

# 2. Load Patient-Level Dataset (which already contains AST and Platelets)
cat("Loading patient-level dataset...\n")
base_data <- read.csv("usfli_nhanes_mash.csv")
cat("Loaded ", nrow(base_data), " MASLD cases.\n")

# 3. Calculate APRI and Risk Strata
# APRI = ((AST / 40) * 100) / Platelets
cat("\nCalculating APRI index...\n")
base_data <- base_data %>%
  mutate(
    # Compute APRI
    APRI = (AST / 40) * 100 / Platelets,
    
    # Classify Risk Strata
    APRI_Group = case_when(
      APRI <= 0.5 ~ "Low Risk (<=0.5)",
      APRI > 0.5 & APRI <= 1.0 ~ "Intermediate Risk (0.5-1.0)",
      APRI > 1.0 & APRI < 2.0 ~ "High Risk (1.0-2.0)",
      APRI >= 2.0 ~ "Cirrhosis (>=2.0)",
      TRUE ~ NA_character_
    ),
    
    # Binary indicators for sub-cohort analysis
    APRI_GT_0_5 = ifelse(APRI > 0.5, 1, 0),
    APRI_GT_1_0 = ifelse(APRI > 1.0, 1, 0),
    APRI_GE_2_0 = ifelse(APRI >= 2.0, 1, 0)
  )

# Write finalized patient-level dataset with both FIB-4 and APRI details
write.csv(base_data, "usfli_nhanes_fibrosis.csv", row.names = FALSE)
cat("Saved detailed patient dataset with APRI to 'usfli_nhanes_fibrosis.csv'\n")

# Filter to complete APRI cases
apri_valid <- base_data %>% filter(!is.na(APRI))
cat("\nMASLD Cases with complete APRI data: ", nrow(apri_valid), " / ", nrow(base_data), "\n")

# 4. Estimate APRI Risk Strata Distribution
cat("\n============================================================\n")
cat("APRI RISK STRATA DISTRIBUTION IN MASLD COHORT\n")
cat("============================================================\n")

# Set up complex survey design
nhanes_design <- svydesign(
  id = ~SDMVPSU, 
  strata = ~SDMVSTRA, 
  weights = ~FastingWeight, 
  data = apri_valid, 
  nest = TRUE
)

# Raw counts
raw_apri_table <- table(apri_valid$APRI_Group, useNA = "no")
raw_apri_total <- sum(raw_apri_table)

# Factor levels coercion
apri_valid$APRI_Group <- factor(apri_valid$APRI_Group, levels = names(raw_apri_table))
nhanes_design_updated <- svydesign(
  id = ~SDMVPSU, 
  strata = ~SDMVSTRA, 
  weights = ~FastingWeight, 
  data = apri_valid, 
  nest = TRUE
)

# Survey weighted prevalence
weighted_apri <- svymean(~APRI_Group, nhanes_design_updated, na.rm = TRUE)
weighted_apri_pct <- as.numeric(weighted_apri) * 100
weighted_apri_se <- as.numeric(SE(weighted_apri)) * 100

apri_summary <- data.frame(
  Stratum = names(raw_apri_table),
  Raw_N = as.numeric(raw_apri_table),
  Raw_Percent = (as.numeric(raw_apri_table) / raw_apri_total) * 100,
  Weighted_Percent = weighted_apri_pct,
  Weighted_SE = weighted_apri_se,
  stringsAsFactors = FALSE
)
print(apri_summary)

# 5. Helper Function for Stratification
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

# 6. Perform Comorbidity Stratification on the Significant Fibrosis Cohort (APRI > 0.5)
cat("\n============================================================\n")
cat("STRATIFICATION OF SIGNIFICANT FIBROSIS COHORT (APRI > 0.5)\n")
cat("============================================================\n")

design_apri_0_5 <- subset(nhanes_design_updated, APRI_GT_0_5 == 1)
cat("Significant Fibrosis Cohort Size (APRI > 0.5 Raw N): ", sum(apri_valid$APRI_GT_0_5 == 1), "\n")

obesity_05 <- get_stratification_stats(design_apri_0_5, "Obesity_Cat", "Obesity/BMI")
diabetes_05 <- get_stratification_stats(design_apri_0_5, "Diabetes_Cat", "Type 2 DM")
hypertension_05 <- get_stratification_stats(design_apri_0_5, "Hypertension_Cat", "Hypertension")
dyslipidemia_05 <- get_stratification_stats(design_apri_0_5, "Dyslipidemia_Cat", "Dyslipidemia")
diet_05 <- get_stratification_stats(design_apri_0_5, "Diet_Cat", "Diet (Caloric Intake)")
pa_05 <- get_stratification_stats(design_apri_0_5, "PA_Cat", "Physical Activity")
smoking_05 <- get_stratification_stats(design_apri_0_5, "Smoking_Cat", "Smoking Status")

design_mortality_05 <- subset(design_apri_0_5, Cycle != "2021-2023" & !is.na(MORTSTAT))
mortality_05 <- get_stratification_stats(design_mortality_05, "Mortality_Cat", "Mortality (1999-2018 cohorts)")

apri_05_strat <- bind_rows(obesity_05, diabetes_05, hypertension_05, dyslipidemia_05, diet_05, pa_05, smoking_05, mortality_05) %>%
  mutate(Cohort = "Significant Fibrosis (APRI > 0.5)")

# 7. Perform Comorbidity Stratification on the Advanced Fibrosis/Cirrhosis Cohort (APRI > 1.0)
cat("\n============================================================\n")
cat("STRATIFICATION OF ADVANCED FIBROSIS/CIRRHOSIS COHORT (APRI > 1.0)\n")
cat("============================================================\n")

design_apri_1_0 <- subset(nhanes_design_updated, APRI_GT_1_0 == 1)
cat("Advanced Fibrosis Cohort Size (APRI > 1.0 Raw N): ", sum(apri_valid$APRI_GT_1_0 == 1), "\n")

obesity_10 <- get_stratification_stats(design_apri_1_0, "Obesity_Cat", "Obesity/BMI")
diabetes_10 <- get_stratification_stats(design_apri_1_0, "Diabetes_Cat", "Type 2 DM")
hypertension_10 <- get_stratification_stats(design_apri_1_0, "Hypertension_Cat", "Hypertension")
dyslipidemia_10 <- get_stratification_stats(design_apri_1_0, "Dyslipidemia_Cat", "Dyslipidemia")
diet_10 <- get_stratification_stats(design_apri_1_0, "Diet_Cat", "Diet (Caloric Intake)")
pa_10 <- get_stratification_stats(design_apri_1_0, "PA_Cat", "Physical Activity")
smoking_10 <- get_stratification_stats(design_apri_1_0, "Smoking_Cat", "Smoking Status")

design_mortality_10 <- subset(design_apri_1_0, Cycle != "2021-2023" & !is.na(MORTSTAT))
mortality_10 <- get_stratification_stats(design_mortality_10, "Mortality_Cat", "Mortality (1999-2018 cohorts)")

apri_10_strat <- bind_rows(obesity_10, diabetes_10, hypertension_10, dyslipidemia_10, diet_10, pa_10, smoking_10, mortality_10) %>%
  mutate(Cohort = "Advanced Fibrosis (APRI > 1.0)")

# 8. Perform Comorbidity Stratification on the Cirrhosis Cohort (APRI >= 2.0)
cat("\n============================================================\n")
cat("STRATIFICATION OF CIRRHOSIS COHORT (APRI >= 2.0)\n")
cat("============================================================\n")

design_apri_2_0 <- subset(nhanes_design_updated, APRI_GE_2_0 == 1)
cat("Cirrhosis Cohort Size (APRI >= 2.0 Raw N): ", sum(apri_valid$APRI_GE_2_0 == 1), "\n")

obesity_20 <- get_stratification_stats(design_apri_2_0, "Obesity_Cat", "Obesity/BMI")
diabetes_20 <- get_stratification_stats(design_apri_2_0, "Diabetes_Cat", "Type 2 DM")
hypertension_20 <- get_stratification_stats(design_apri_2_0, "Hypertension_Cat", "Hypertension")
dyslipidemia_20 <- get_stratification_stats(design_apri_2_0, "Dyslipidemia_Cat", "Dyslipidemia")
diet_20 <- get_stratification_stats(design_apri_2_0, "Diet_Cat", "Diet (Caloric Intake)")
pa_20 <- get_stratification_stats(design_apri_2_0, "PA_Cat", "Physical Activity")
smoking_20 <- get_stratification_stats(design_apri_2_0, "Smoking_Cat", "Smoking Status")

design_mortality_20 <- subset(design_apri_2_0, Cycle != "2021-2023" & !is.na(MORTSTAT))
mortality_20 <- get_stratification_stats(design_mortality_20, "Mortality_Cat", "Mortality (1999-2018 cohorts)")

apri_20_strat <- bind_rows(obesity_20, diabetes_20, hypertension_20, dyslipidemia_20, diet_20, pa_20, smoking_20, mortality_20) %>%
  mutate(Cohort = "Cirrhosis (APRI >= 2.0)")

# Combine and save results
final_apri_strat <- bind_rows(apri_05_strat, apri_10_strat, apri_20_strat)
write.csv(final_apri_strat, "apri_comorbidities_stratification.csv", row.names = FALSE)
cat("\nSaved APRI stratification summary to 'apri_comorbidities_stratification.csv'\n")

# Print combined outputs
print(final_apri_strat)
cat("APRI analysis completed successfully!\n")
