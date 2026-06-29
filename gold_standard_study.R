# ==============================================================================
# NHANES Gold Standard Liver Study (2017-2023)
# ==============================================================================
# This script defines MASLD, MASH, and Cirrhosis using transient elastography 
# (FibroScan) gold standard definitions for the 2017-2023 cohorts, combines 
# survey weights, and stratifies each cohort across comorbidities and lifestyle factors.
# ==============================================================================

# 1. Load Required Libraries
cat("Loading required libraries...\n")
library(nhanesA)
library(dplyr)
library(survey)
library(readr)

# Set survey option to adjust for lonely PSUs
options(survey.lonely.psu = "adjust")

# Configurations
CACHE_DIR <- "nhanes_cache"

# Helper Function: Check "Yes" or "1" in a rowwise fashion for survey columns
check_yes_row <- function(df_subset) {
  if (ncol(df_subset) == 0) return(rep(FALSE, nrow(df_subset)))
  logical_matrix <- sapply(df_subset, function(col) {
    val_char <- as.character(col)
    val_char == "Yes" | val_char == "1"
  })
  if (is.vector(logical_matrix)) {
    logical_matrix <- matrix(logical_matrix, ncol = 1)
  }
  rowSums(logical_matrix, na.rm = TRUE) > 0
}

# Helper Function: Load and Cache NHANES Datasets
load_nhanes_cached <- function(table_name, cache_dir = CACHE_DIR) {
  cache_file = file.path(cache_dir, paste0(table_name, ".rds"))
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

# Helper Function: Load and Cache NCHS Mortality Files
load_mortality_cached <- function(cycle_name, cache_dir = CACHE_DIR) {
  filename_map <- list(
    "2017-2020 (Pre-Pandemic)" = "NHANES_2017_2018_MORT_2019_PUBLIC.dat" # 2017-18 mortality for 2017-20 cohort
  )
  if (!(cycle_name %in% names(filename_map))) return(NULL)
  fname <- filename_map[[cycle_name]]
  dest_file <- file.path(cache_dir, fname)
  if (!file.exists(dest_file)) {
    url <- paste0("https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/", fname)
    cat("  Downloading mortality file:", fname, "...\n")
    tryCatch({
      download.file(url, dest_file, mode = "wb", quiet = TRUE)
    }, error = function(e) {
      cat("  Warning: failed to download mortality file:", e$message, "\n")
      return(NULL)
    })
  }
  if (!file.exists(dest_file)) return(NULL)
  df <- read_fwf(
    file = dest_file,
    col_types = cols(
      SEQN = col_integer(),
      ELIGSTAT = col_integer(),
      MORTSTAT = col_integer(),
      UCOD_LEADING = col_integer(),
      DIABETES = col_integer(),
      HYPERTEN = col_integer(),
      PERMTH_INT = col_integer(),
      PERMTH_EXM = col_integer()
    ),
    fwf_cols(
      SEQN         = c(1, 14),
      ELIGSTAT     = c(15, 15),
      MORTSTAT     = c(16, 16),
      UCOD_LEADING = c(17, 19),
      DIABETES     = c(20, 20),
      HYPERTEN     = c(21, 21),
      PERMTH_INT   = c(43, 45),
      PERMTH_EXM   = c(46, 48)
    ),
    na = c("", ".")
  )
  return(df)
}

# 3. Define Recent Cycles (2017-2023)
recent_cycles <- list(
  "2017-2020 (Pre-Pandemic)" = list(
    demo="P_DEMO", bmx="P_BMX", biopro="P_BIOPRO", glu="P_GLU", ins="P_INS",
    diq="P_DIQ", ghb="P_GHB", bpq="P_BPQ", bpx="P_BPX", chol="P_TCHOL",
    hdl="P_HDL", diet="P_DR1TOT", pa="P_PAQ", smq="P_SMQ", lux="P_LUX",
    weight_var="WTSAFPRP"
  ),
  "2021-2023" = list(
    demo="DEMO_L", bmx="BMX_L", biopro="BIOPRO_L", glu="GLU_L", ins="INS_L",
    diq="DIQ_L", ghb="GHB_L", bpq="BPQ_L", bpx="BPX_L", chol="TCHOL_L",
    hdl="HDL_L", diet="DR1TOT_L", pa="PAQ_L", smq="SMQ_L", lux="LUX_L",
    weight_var="WTSAF2YR"
  )
)

# 4. Fetch and Clean Datasets Cycle-by-Cycle
clean_list <- list()

for (cyc_name in names(recent_cycles)) {
  cat("\nProcessing cycle:", cyc_name, "...\n")
  info <- recent_cycles[[cyc_name]]
  
  # Load core tables
  demo <- load_nhanes_cached(info$demo)
  bmx <- load_nhanes_cached(info$bmx)
  biopro <- load_nhanes_cached(info$biopro)
  glu <- load_nhanes_cached(info$glu)
  ins <- load_nhanes_cached(info$ins)
  diq <- load_nhanes_cached(info$diq)
  ghb <- load_nhanes_cached(info$ghb)
  bpq <- load_nhanes_cached(info$bpq)
  bpx <- load_nhanes_cached(info$bpx)
  chol <- load_nhanes_cached(info$chol)
  hdl <- load_nhanes_cached(info$hdl)
  diet <- load_nhanes_cached(info$diet)
  pa  <- load_nhanes_cached(info$pa)
  smq <- load_nhanes_cached(info$smq)
  lux <- load_nhanes_cached(info$lux)
  mort <- load_mortality_cached(cyc_name)
  
  if (is.null(demo) || is.null(lux) || is.null(biopro)) next
  
  # Select variables
  # Demographics: Age >= 20
  demo_clean <- demo %>%
    filter(RIDAGEYR >= 20) %>%
    select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH1, SDMVPSU, SDMVSTRA)
  
  # Elastography (LUX): LUXCAPM (CAP) and LUXSMED (LSM)
  lux_clean <- lux %>%
    select(SEQN, LUXCAPM, LUXSMED) %>%
    filter(!is.na(LUXCAPM), !is.na(LUXSMED))
  
  # AST / ALT
  biopro_clean <- biopro %>%
    select(SEQN, AST = LBXSASSI, ALT = LBXSATSI)
  
  # BMI
  bmx_clean <- bmx %>%
    select(SEQN, BMXBMI)
  
  # Diabetes and Fasting Glucose
  diq_clean <- if (!is.null(diq) && "DIQ010" %in% colnames(diq)) diq %>% select(SEQN, DIQ010) else data.frame(SEQN=integer(), DIQ010=numeric())
  ghb_clean <- if (!is.null(ghb) && "LBXGH" %in% colnames(ghb)) ghb %>% select(SEQN, LBXGH) else data.frame(SEQN=integer(), LBXGH=numeric())
  
  # Fasting weight is in the fasting glucose file
  if (!is.null(glu) && "LBXGLU" %in% colnames(glu)) {
    glu_cols <- intersect(colnames(glu), c("LBXGLU", info$weight_var))
    glu_clean <- glu %>% 
      select(SEQN, Glucose = LBXGLU, FastingWeight = !!sym(info$weight_var)) %>%
      filter(!is.na(FastingWeight))
  } else {
    glu_clean <- data.frame(SEQN=integer(), Glucose=numeric(), FastingWeight=numeric())
  }
  
  # Blood pressure
  bpq_cols <- intersect(colnames(bpq), c("BPQ020", "BPQ040A", "BPQ050A", "BPQ080", "BPQ090D"))
  bpq_clean <- if (!is.null(bpq)) bpq %>% select(SEQN, one_of(bpq_cols)) else data.frame(SEQN=integer())
  
  if (!is.null(bpx)) {
    sy_cols <- intersect(colnames(bpx), c("BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4"))
    di_cols <- intersect(colnames(bpx), c("BPXDI1", "BPXDI2", "BPXDI3", "BPXDI4"))
    bpx_processed <- bpx %>%
      mutate(
        MeanSystolic = rowMeans(select(., one_of(sy_cols)), na.rm = TRUE),
        MeanDiastolic = rowMeans(select(., one_of(di_cols)), na.rm = TRUE)
      ) %>%
      mutate(
        MeanSystolic = ifelse(is.nan(MeanSystolic), NA, MeanSystolic),
        MeanDiastolic = ifelse(is.nan(MeanDiastolic), NA, MeanDiastolic)
      ) %>%
      select(SEQN, MeanSystolic, MeanDiastolic)
  } else {
    bpx_processed <- data.frame(SEQN=integer(), MeanSystolic=numeric(), MeanDiastolic=numeric())
  }
  
  # Lipids
  chol_clean <- if (!is.null(chol) && "LBXTC" %in% colnames(chol)) chol %>% select(SEQN, LBXTC) else data.frame(SEQN=integer(), LBXTC=numeric())
  hdl_col <- intersect(colnames(hdl), c("LBDHDD", "LBXHDD", "LBDHDL"))
  hdl_clean <- if (!is.null(hdl) && length(hdl_col) > 0) hdl %>% select(SEQN, HDL = !!sym(hdl_col[1])) else data.frame(SEQN=integer(), HDL=numeric())
  
  # Diet (Calorie)
  diet_col <- intersect(colnames(diet), c("DR1TKCAL", "DRXTKCAL"))
  diet_clean <- if (!is.null(diet) && length(diet_col) > 0) diet %>% select(SEQN, Kcal = !!sym(diet_col[1])) else data.frame(SEQN=integer(), Kcal=numeric())
  
  # Physical Activity
  if (!is.null(pa)) {
    pa_cols_late <- intersect(colnames(pa), c("PAQ605", "PAQ620", "PAQ650", "PAQ665"))
    pa_processed <- pa %>%
      mutate(
        PA_Active = as.numeric(check_yes_row(select(., one_of(pa_cols_late))))
      ) %>%
      select(SEQN, PA_Active)
  } else {
    pa_processed <- data.frame(SEQN=integer(), PA_Active=numeric())
  }
  
  # Smoking
  smq_clean <- if (!is.null(smq)) {
    smq_cols <- intersect(colnames(smq), c("SMQ020", "SMQ040"))
    smq %>% select(SEQN, one_of(smq_cols))
  } else {
    data.frame(SEQN=integer())
  }
  
  # Mortality
  mort_clean <- if (!is.null(mort)) mort %>% select(SEQN, ELIGSTAT, MORTSTAT) else data.frame(SEQN=integer(), ELIGSTAT=numeric(), MORTSTAT=numeric())
  
  # Merge - must be inner_join with glu_clean to keep only fasting subsample
  cyc_merged <- demo_clean %>%
    inner_join(glu_clean, by = "SEQN") %>%
    inner_join(lux_clean, by = "SEQN") %>%
    left_join(biopro_clean, by = "SEQN") %>%
    left_join(bmx_clean, by = "SEQN") %>%
    left_join(diq_clean, by = "SEQN") %>%
    left_join(ghb_clean, by = "SEQN") %>%
    left_join(bpq_clean, by = "SEQN") %>%
    left_join(bpx_processed, by = "SEQN") %>%
    left_join(chol_clean, by = "SEQN") %>%
    left_join(hdl_clean, by = "SEQN") %>%
    left_join(diet_clean, by = "SEQN") %>%
    left_join(pa_processed, by = "SEQN") %>%
    left_join(smq_clean, by = "SEQN") %>%
    left_join(mort_clean, by = "SEQN") %>%
    mutate(Cycle = cyc_name)
  
  clean_list[[cyc_name]] <- cyc_merged
}

# Combine all recent cohorts (2017-2023)
gold_raw <- bind_rows(clean_list)
cat("\nTotal fasting participants with valid demographics and elastography: ", nrow(gold_raw), "\n")

# 5. Classify Comorbidities
cat("\nCategorizing comorbidity and lifestyle factors...\n")

# Diet tertiles (calculated within this recent cohort)
kcal_tertiles <- quantile(gold_raw$Kcal, probs = c(1/3, 2/3), na.rm = TRUE)

gold_raw <- gold_raw %>%
  mutate(
    # Obesity Category
    Obesity_Cat = case_when(
      BMXBMI < 18.5 ~ "Underweight (<18.5)",
      BMXBMI >= 18.5 & BMXBMI < 25 ~ "Normal (18.5-24.9)",
      BMXBMI >= 25 & BMXBMI < 30 ~ "Overweight (25.0-29.9)",
      BMXBMI >= 30 ~ "Obese (>=30.0)",
      TRUE ~ NA_character_
    ),
    
    # Type 2 DM Category
    Diabetes_Cat = case_when(
      as.character(DIQ010) == "Yes" | as.character(DIQ010) == "1" | Glucose >= 126 | LBXGH >= 6.5 ~ "Diabetes",
      as.character(DIQ010) == "No" | as.character(DIQ010) == "2" | Glucose < 126 | LBXGH < 6.5 ~ "No Diabetes",
      TRUE ~ NA_character_
    ),
    
    # Hypertension Category
    Hypertension_Cat = case_when(
      as.character(BPQ020) == "Yes" | as.character(BPQ020) == "1" | 
        as.character(BPQ040A) == "Yes" | as.character(BPQ040A) == "1" | 
        as.character(BPQ050A) == "Yes" | as.character(BPQ050A) == "1" | 
        MeanSystolic >= 130 | MeanDiastolic >= 85 ~ "Hypertension",
      as.character(BPQ020) == "No" | as.character(BPQ020) == "2" | 
        MeanSystolic < 130 | MeanDiastolic < 85 ~ "No Hypertension",
      TRUE ~ NA_character_
    ),
    
    # Dyslipidemia Category
    gender_male = (RIAGENDR == "Male" | RIAGENDR == "1" | as.numeric(as.factor(RIAGENDR)) == 1),
    Dyslipidemia_Cat = case_when(
      as.character(BPQ080) == "Yes" | as.character(BPQ080) == "1" | 
        as.character(BPQ090D) == "Yes" | as.character(BPQ090D) == "1" | 
        LBXTC >= 200 | (gender_male & HDL < 40) | (!gender_male & HDL < 50) ~ "Dyslipidemia",
      as.character(BPQ080) == "No" | as.character(BPQ080) == "2" | 
        LBXTC < 200 | (gender_male & HDL >= 40) | (!gender_male & HDL >= 50) ~ "No Dyslipidemia",
      TRUE ~ NA_character_
    ),
    
    # Diet Category
    Diet_Cat = case_when(
      Kcal <= kcal_tertiles[1] ~ "Low Calorie (Tertile 1)",
      Kcal > kcal_tertiles[1] & Kcal <= kcal_tertiles[2] ~ "Medium Calorie (Tertile 2)",
      Kcal > kcal_tertiles[2] ~ "High Calorie (Tertile 3)",
      TRUE ~ NA_character_
    ),
    
    # Physical Activity Category
    PA_Cat = case_when(
      PA_Active == 1 ~ "Active",
      PA_Active == 0 ~ "Inactive",
      TRUE ~ NA_character_
    ),
    
    # Smoking Category
    Smoking_Cat = case_when(
      as.character(SMQ020) == "No" | as.character(SMQ020) == "2" ~ "Never Smoker",
      (as.character(SMQ020) == "Yes" | as.character(SMQ020) == "1") & (as.character(SMQ040) == "Not at all" | as.character(SMQ040) == "3") ~ "Former Smoker",
      (as.character(SMQ020) == "Yes" | as.character(SMQ020) == "1") & (as.character(SMQ040) == "Every day" | as.character(SMQ040) == "Some days" | as.character(SMQ040) == "1" | as.character(SMQ040) == "2") ~ "Current Smoker",
      TRUE ~ NA_character_
    ),
    
    # Mortality Category (only for 2017-2020 cycle, 2021-2023 set to NA)
    Mortality_Cat = case_when(
      MORTSTAT == 1 ~ "Deceased",
      MORTSTAT == 0 ~ "Alive",
      TRUE ~ NA_character_
    )
  )

# 6. Apply Gold Standard Definitions
# MASLD: CAP >= 285 + at least one Cardiometabolic Risk Factor (CMR)
# MASH: FAST score >= 0.35 (optimal screening cutoff for active MASH: NAS >= 4, Fibrosis >= F2)
# Cirrhosis: Liver Stiffness Measurement (LSM) > 13.6 kPa
cat("\nApplying Gold Standard Liver Disease Definitions...\n")

gold_raw <- gold_raw %>%
  mutate(
    # Cardiometabolic Risk (CMR)
    CMR = ifelse(BMXBMI >= 25 | Diabetes_Cat == "Diabetes" | Glucose >= 100 | LBXGH >= 5.7 | 
                   Hypertension_Cat == "Hypertension" | Dyslipidemia_Cat == "Dyslipidemia", 1, 0),
    
    # Gold MASLD
    MASLD_Gold = ifelse(LUXCAPM >= 285 & CMR == 1, 1, 0),
    
    # FAST score formula: x = -1.65 + 1.07*ln(LSM) + (2.66e-8)*CAP^3 - 63.3/AST
    # FAST = exp(x)/(1+exp(x))
    L = -1.65 + 1.07 * log(LUXSMED) + (2.66e-8) * (LUXCAPM^3) - 63.3 / AST,
    FAST = exp(L) / (1 + exp(L)),
    
    # Gold MASH (Active MASH)
    MASH_Gold = ifelse(MASLD_Gold == 1 & FAST >= 0.35, 1, 0),
    
    # Gold Cirrhosis
    Cirrhosis_Gold = ifelse(MASLD_Gold == 1 & LUXSMED > 13.6, 1, 0)
  )

# Write finalized gold standard patient-level dataset
write.csv(gold_raw, "usfli_nhanes_gold.csv", row.names = FALSE)
cat("Saved detailed gold standard dataset to 'usfli_nhanes_gold.csv'\n")

# Combine survey weights for recent 6-year period (2017-2023)
# Scaled weight: 4-year cycle (2017-2020) scaled by 4/6 = 2/3
# Scaled weight: 2-year cycle (2021-2023) scaled by 2/6 = 1/3
gold_raw <- gold_raw %>%
  mutate(
    ScaledWeight = ifelse(Cycle == "2017-2020 (Pre-Pandemic)", FastingWeight * (2/3), FastingWeight * (1/3))
  )

# 7. Print Cohort N Summary
cat("\n============================================================\n")
cat("GOLD STANDARD COHORT POPULATION (N = ", nrow(gold_raw), " with valid elastography)\n")
cat("============================================================\n")
cat("MASLD (CAP >= 285 + CMR):       Raw N = ", sum(gold_raw$MASLD_Gold == 1, na.rm=TRUE), "\n")
cat("MASH (MASLD + FAST >= 0.35):    Raw N = ", sum(gold_raw$MASH_Gold == 1, na.rm=TRUE), "\n")
cat("Cirrhosis (MASLD + LSM > 13.6): Raw N = ", sum(gold_raw$Cirrhosis_Gold == 1, na.rm=TRUE), "\n")

# Establish survey design object
gold_design <- svydesign(
  id = ~SDMVPSU, 
  strata = ~SDMVSTRA, 
  weights = ~ScaledWeight, 
  data = gold_raw, 
  nest = TRUE
)

# 8. Helper Function for Stratification
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
      Stratum = "bullet", # placeholder, will resolve to sname later
      ActualStratum = sname,
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
    Stratum = "bullet",
    ActualStratum = "Missing / Unknown",
    Raw_N = missing_n,
    Raw_Percent = (missing_n / nrow(df_data)) * 100,
    Weighted_Percent = NA,
    Weighted_SE = NA,
    stringsAsFactors = FALSE
  )
  
  # Adjust dataframe to return Stratum as ActualStratum (dropping extra columns)
  res_df <- do.call(rbind, result_rows)
  res_df$Stratum <- res_df$ActualStratum
  res_df$ActualStratum <- NULL
  return(res_df)
}

# 9. Perform Comorbidity Stratification on the Three Gold Standard Cohorts

# A. MASLD (Gold)
cat("\nRunning stratification for Gold MASLD...\n")
design_masld <- subset(gold_design, MASLD_Gold == 1)
obesity_masld <- get_stratification_stats(design_masld, "Obesity_Cat", "Obesity/BMI")
diabetes_masld <- get_stratification_stats(design_masld, "Diabetes_Cat", "Type 2 DM")
hypertension_masld <- get_stratification_stats(design_masld, "Hypertension_Cat", "Hypertension")
dyslipidemia_masld <- get_stratification_stats(design_masld, "Dyslipidemia_Cat", "Dyslipidemia")
diet_masld <- get_stratification_stats(design_masld, "Diet_Cat", "Diet (Caloric Intake)")
pa_masld <- get_stratification_stats(design_masld, "PA_Cat", "Physical Activity")
smoking_masld <- get_stratification_stats(design_masld, "Smoking_Cat", "Smoking Status")

design_mortality_masld <- subset(design_masld, Cycle == "2017-2020 (Pre-Pandemic)" & !is.na(MORTSTAT))
mortality_masld <- get_stratification_stats(design_mortality_masld, "Mortality_Cat", "Mortality (2017-2020 cohort only)")

strat_masld <- bind_rows(obesity_masld, diabetes_masld, hypertension_masld, dyslipidemia_masld, diet_masld, pa_masld, smoking_masld, mortality_masld) %>%
  mutate(Cohort = "Gold MASLD (CAP >= 285 + CMR)")

# B. MASH (Gold)
cat("Running stratification for Gold MASH...\n")
design_mash <- subset(gold_design, MASH_Gold == 1)
obesity_mash <- get_stratification_stats(design_mash, "Obesity_Cat", "Obesity/BMI")
diabetes_mash <- get_stratification_stats(design_mash, "Diabetes_Cat", "Type 2 DM")
hypertension_mash <- get_stratification_stats(design_mash, "Hypertension_Cat", "Hypertension")
dyslipidemia_mash <- get_stratification_stats(design_mash, "Dyslipidemia_Cat", "Dyslipidemia")
diet_mash <- get_stratification_stats(design_mash, "Diet_Cat", "Diet (Caloric Intake)")
pa_mash <- get_stratification_stats(design_mash, "PA_Cat", "Physical Activity")
smoking_mash <- get_stratification_stats(design_mash, "Smoking_Cat", "Smoking Status")

design_mortality_mash <- subset(design_mash, Cycle == "2017-2020 (Pre-Pandemic)" & !is.na(MORTSTAT))
mortality_mash <- get_stratification_stats(design_mortality_mash, "Mortality_Cat", "Mortality (2017-2020 cohort only)")

strat_mash <- bind_rows(obesity_mash, diabetes_mash, hypertension_mash, dyslipidemia_mash, diet_mash, pa_mash, smoking_mash, mortality_mash) %>%
  mutate(Cohort = "Gold MASH (FAST >= 0.35)")

# C. Cirrhosis (Gold)
cat("Running stratification for Gold Cirrhosis...\n")
design_cirrhosis <- subset(gold_design, Cirrhosis_Gold == 1)
obesity_cirr <- get_stratification_stats(design_cirrhosis, "Obesity_Cat", "Obesity/BMI")
diabetes_cirr <- get_stratification_stats(design_cirrhosis, "Diabetes_Cat", "Type 2 DM")
hypertension_cirr <- get_stratification_stats(design_cirrhosis, "Hypertension_Cat", "Hypertension")
dyslipidemia_cirr <- get_stratification_stats(design_cirrhosis, "Dyslipidemia_Cat", "Dyslipidemia")
diet_cirr <- get_stratification_stats(design_cirrhosis, "Diet_Cat", "Diet (Caloric Intake)")
pa_cirr <- get_stratification_stats(design_cirrhosis, "PA_Cat", "Physical Activity")
smoking_cirr <- get_stratification_stats(design_cirrhosis, "Smoking_Cat", "Smoking Status")

design_mortality_cirr <- subset(design_cirrhosis, Cycle == "2017-2020 (Pre-Pandemic)" & !is.na(MORTSTAT))
mortality_cirr <- get_stratification_stats(design_mortality_cirr, "Mortality_Cat", "Mortality (2017-2020 cohort only)")

strat_cirr <- bind_rows(obesity_cirr, diabetes_cirr, hypertension_cirr, dyslipidemia_cirr, diet_cirr, pa_cirr, smoking_cirr, mortality_cirr) %>%
  mutate(Cohort = "Gold Cirrhosis (LSM > 13.6)")

# Combine and save results
gold_stratification <- bind_rows(strat_masld, strat_mash, strat_cirr)
write.csv(gold_stratification, "gold_standard_stratification.csv", row.names = FALSE)
cat("\nSaved gold standard stratification summary to 'gold_standard_stratification.csv'\n")

# Print output summary
print(gold_stratification)
cat("Gold standard liver analysis completed successfully!\n")
