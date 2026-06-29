# ==============================================================================
# NHANES MASLD Cohort Comorbidity and Lifestyle Stratification
# ==============================================================================
# This script loads the cohort of 8,673 MASLD cases, downloads additional 
# questionnaire and laboratory datasets for comorbidities (obesity, diabetes, 
# hypertension, dyslipidemia) and lifestyle factors (diet, physical activity, 
# smoking, mortality), and calculates the stratified sample sizes and weighted 
# prevalence rates.
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
if (!dir.exists(CACHE_DIR)) {
  dir.create(CACHE_DIR)
}

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

# 2. Define Mappings for New Tables across the 11 Cohort Cycles
comorb_cycles <- list(
  "1999-2000" = list(diq="DIQ",     ghb="GHB",     bpq="BPQ",     bpx="BPX",     chol="LAB13",   hdl="LAB13",   diet="DRXTOT",   pa="PAQ",     smq="SMQ"),
  "2001-2002" = list(diq="DIQ_B",   ghb="GHB_B",   bpq="BPQ_B",   bpx="BPX_B",   chol="L13_B",   hdl="L13_B",   diet="DRXTOT_B", pa="PAQ_B",   smq="SMQ_B"),
  "2003-2004" = list(diq="DIQ_C",   ghb="GHB_C",   bpq="BPQ_C",   bpx="BPX_C",   chol="L13_C",   hdl="L13_C",   diet="DR1TOT_C", pa="PAQ_C",   smq="SMQ_C"),
  "2005-2006" = list(diq="DIQ_D",   ghb="GHB_D",   bpq="BPQ_D",   bpx="BMX_D",   chol="TCHOL_D", hdl="HDL_D",   diet="DR1TOT_D", pa="PAQ_D",   smq="SMQ_D"), # Note: bpx table is named BMX_D by CDC, wait: is BP exam BMX_D? No! The blood pressure exam is BPX_D! Wait, let's keep bpx="BPX_D".
  "2007-2008" = list(diq="DIQ_E",   ghb="GHB_E",   bpq="BPQ_E",   bpx="BPX_E",   chol="TCHOL_E", hdl="HDL_E",   diet="DR1TOT_E", pa="PAQ_E",   smq="SMQ_E"),
  "2009-2010" = list(diq="DIQ_F",   ghb="GHB_F",   bpq="BPQ_F",   bpx="BPX_F",   chol="TCHOL_F", hdl="HDL_F",   diet="DR1TOT_F", pa="PAQ_F",   smq="SMQ_F"),
  "2011-2012" = list(diq="DIQ_G",   ghb="GHB_G",   bpq="BPQ_G",   bpx="BPX_G",   chol="TCHOL_G", hdl="HDL_G",   diet="DR1TOT_G", pa="PAQ_G",   smq="SMQ_G"),
  "2013-2014" = list(diq="DIQ_H",   ghb="GHB_H",   bpq="BPQ_H",   bpx="BPX_H",   chol="TCHOL_H", hdl="HDL_H",   diet="DR1TOT_H", pa="PAQ_H",   smq="SMQ_H"),
  "2015-2016" = list(diq="DIQ_I",   ghb="GHB_I",   bpq="BPQ_I",   bpx="BPX_I",   chol="TCHOL_I", hdl="HDL_I",   diet="DR1TOT_I", pa="PAQ_I",   smq="SMQ_I"),
  "2017-2020 (Pre-Pandemic)" = list(diq="P_DIQ",   ghb="P_GHB",   bpq="P_BPQ",   bpx="P_BPX",   chol="P_TCHOL", hdl="P_HDL",   diet="P_DR1TOT", pa="P_PAQ",   smq="P_SMQ"),
  "2021-2023" = list(diq="DIQ_L",   ghb="GHB_L",   bpq="BPQ_L",   bpx="BPX_L",   chol="TCHOL_L", hdl="HDL_L",   diet="DR1TOT_L", pa="PAQ_L",   smq="SMQ_L")
)

# Fix: Ensure bpx is BPX_D in 2005-2006
comorb_cycles[["2005-2006"]]$bpx <- "BPX_D"

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

# Helper Function: Load and Cache Fixed-Width NCHS Linked Mortality Files
load_mortality_cached <- function(cycle_name, cache_dir = CACHE_DIR) {
  filename_map <- list(
    "1999-2000" = "NHANES_1999_2000_MORT_2019_PUBLIC.dat",
    "2001-2002" = "NHANES_2001_2002_MORT_2019_PUBLIC.dat",
    "2003-2004" = "NHANES_2003_2004_MORT_2019_PUBLIC.dat",
    "2005-2006" = "NHANES_2005_2006_MORT_2019_PUBLIC.dat",
    "2007-2008" = "NHANES_2007_2008_MORT_2019_PUBLIC.dat",
    "2009-2010" = "NHANES_2009_2010_MORT_2019_PUBLIC.dat",
    "2011-2012" = "NHANES_2011_2012_MORT_2019_PUBLIC.dat",
    "2013-2014" = "NHANES_2013_2014_MORT_2019_PUBLIC.dat",
    "2015-2016" = "NHANES_2015_2016_MORT_2019_PUBLIC.dat",
    "2017-2018" = "NHANES_2017_2018_MORT_2019_PUBLIC.dat",
    "2017-2020 (Pre-Pandemic)" = "NHANES_2017_2018_MORT_2019_PUBLIC.dat" # Merge 2017-18 NDI linkage for 2017-20 cohort
  )
  
  if (!(cycle_name %in% names(filename_map))) {
    return(NULL) # For 2021-2023, no public mortality file is available
  }
  
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

# 3. Load MASLD Cohort File
cat("Loading MASLD cohort from 'usfli_nhanes_merged.csv'...\n")
merged_base <- read.csv("usfli_nhanes_merged.csv")

# Filter to the cohort of 8,673 MASLD cases: USFLI > 30 and excluding the 2017-2018 standalone cycle
masld_cohort <- merged_base %>% 
  filter(Cycle != "2017-2018", USFLI_GT_30 == 1)

cat("Target MASLD Cohort Size: ", nrow(masld_cohort), " rows (expected: 8673)\n")

# 4. Fetch and Process Comorbidity and Lifestyle Data for Each Cycle
comorb_list <- list()

for (cyc_name in names(comorb_cycles)) {
  cat("\nProcessing cycle:", cyc_name, "...\n")
  info <- comorb_cycles[[cyc_name]]
  
  # Filter cohort to this cycle
  cyc_masld <- masld_cohort %>% filter(Cycle == cyc_name)
  if (nrow(cyc_masld) == 0) next
  
  # Load tables from cache
  diq <- load_nhanes_cached(info$diq)
  ghb <- load_nhanes_cached(info$ghb)
  bpq <- load_nhanes_cached(info$bpq)
  bpx <- load_nhanes_cached(info$bpx)
  chol <- load_nhanes_cached(info$chol)
  hdl <- load_nhanes_cached(info$hdl)
  diet <- load_nhanes_cached(info$diet)
  pa  <- load_nhanes_cached(info$pa)
  smq <- load_nhanes_cached(info$smq)
  mort <- load_mortality_cached(cyc_name)
  
  # Load BMX to extract BMXBMI (Obesity)
  bmx_tbl_name <- if (cyc_name == "2017-2020 (Pre-Pandemic)") "P_BMX" else if (cyc_name == "2021-2023") "BMX_L" else paste0("BMX", sub(".*-20", "_", cyc_name))
  # Standardize suffix mapping
  if (cyc_name == "1999-2000") bmx_tbl_name <- "BMX"
  if (cyc_name == "2001-2002") bmx_tbl_name <- "BMX_B"
  if (cyc_name == "2003-2004") bmx_tbl_name <- "BMX_C"
  if (cyc_name == "2005-2006") bmx_tbl_name <- "BMX_D"
  if (cyc_name == "2007-2008") bmx_tbl_name <- "BMX_E"
  if (cyc_name == "2009-2010") bmx_tbl_name <- "BMX_F"
  if (cyc_name == "2011-2012") bmx_tbl_name <- "BMX_G"
  if (cyc_name == "2013-2014") bmx_tbl_name <- "BMX_H"
  if (cyc_name == "2015-2016") bmx_tbl_name <- "BMX_I"
  
  bmx <- load_nhanes_cached(bmx_tbl_name)
  
  # Join and select specific columns
  # BMX: BMXBMI
  bmx_clean <- if (!is.null(bmx) && "BMXBMI" %in% colnames(bmx)) bmx %>% select(SEQN, BMXBMI) else data.frame(SEQN=integer(), BMXBMI=numeric())
  
  # DIQ: DIQ010 (ever told you have diabetes)
  diq_clean <- if (!is.null(diq) && "DIQ010" %in% colnames(diq)) diq %>% select(SEQN, DIQ010) else data.frame(SEQN=integer(), DIQ010=numeric())
  
  # GHB: LBXGH (glycohemoglobin HbA1c)
  ghb_clean <- if (!is.null(ghb) && "LBXGH" %in% colnames(ghb)) ghb %>% select(SEQN, LBXGH) else data.frame(SEQN=integer(), LBXGH=numeric())
  
  # BPQ: BPQ020 (told high BP), BPQ040A (taking BP med), BPQ050A (now taking BP med), BPQ080 (told high chol), BPQ090D (now taking chol med)
  bpq_cols <- intersect(colnames(bpq), c("BPQ020", "BPQ030", "BPQ040A", "BPQ050A", "BPQ080", "BPQ090D"))
  bpq_clean <- if (!is.null(bpq)) bpq %>% select(SEQN, one_of(bpq_cols)) else data.frame(SEQN=integer())
  
  # BPX: measured systolic and diastolic BP
  # Average the non-missing systolic (BPXSY1-4) and diastolic (BPXDI1-4) values
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
  
  # Total Cholesterol: LBXTC
  chol_clean <- if (!is.null(chol) && "LBXTC" %in% colnames(chol)) chol %>% select(SEQN, LBXTC) else data.frame(SEQN=integer(), LBXTC=numeric())
  
  # HDL Cholesterol: direct HDL (LBDHDL or LBXHDD or LBDHDD)
  hdl_col <- intersect(colnames(hdl), c("LBDHDD", "LBXHDD", "LBDHDL"))
  hdl_clean <- if (!is.null(hdl) && length(hdl_col) > 0) hdl %>% select(SEQN, HDL = !!sym(hdl_col[1])) else data.frame(SEQN=integer(), HDL=numeric())
  
  # Diet (Caloric intake): DR1TKCAL or DRXTKCAL
  diet_col <- intersect(colnames(diet), c("DR1TKCAL", "DRXTKCAL"))
  diet_clean <- if (!is.null(diet) && length(diet_col) > 0) diet %>% select(SEQN, Kcal = !!sym(diet_col[1])) else data.frame(SEQN=integer(), Kcal=numeric())
  
  # Physical Activity
  # Harmonize active physical activity indicator using check_yes_row
  if (!is.null(pa)) {
    pa_cols_early <- intersect(colnames(pa), c("PAD020", "PAD080"))
    pa_cols_late <- intersect(colnames(pa), c("PAQ605", "PAQ620", "PAQ650", "PAQ665"))
    
    if (length(pa_cols_early) > 0) {
      pa_processed <- pa %>%
        mutate(
          PA_Active = as.numeric(check_yes_row(select(., one_of(pa_cols_early))))
        ) %>%
        select(SEQN, PA_Active)
    } else if (length(pa_cols_late) > 0) {
      pa_processed <- pa %>%
        mutate(
          PA_Active = as.numeric(check_yes_row(select(., one_of(pa_cols_late))))
        ) %>%
        select(SEQN, PA_Active)
    } else {
      pa_processed <- data.frame(SEQN=integer(), PA_Active=numeric())
    }
  } else {
    pa_processed <- data.frame(SEQN=integer(), PA_Active=numeric())
  }
  
  # Smoking: SMQ020 (smoked 100 cigs), SMQ040 (smoke now)
  smq_clean <- if (!is.null(smq)) {
    smq_cols <- intersect(colnames(smq), c("SMQ020", "SMQ040"))
    smq %>% select(SEQN, one_of(smq_cols))
  } else {
    data.frame(SEQN=integer())
  }
  
  # Mortality variables
  mort_clean <- if (!is.null(mort)) mort %>% select(SEQN, ELIGSTAT, MORTSTAT, PERMTH_INT, PERMTH_EXM) else data.frame(SEQN=integer(), ELIGSTAT=numeric(), MORTSTAT=numeric())
  
  # Merge comorbidity/lifestyle files into the cycle MASLD cohort
  cyc_masld_merged <- cyc_masld %>%
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
    left_join(mort_clean, by = "SEQN")
  
  comorb_list[[cyc_name]] <- cyc_masld_merged
}

# Combine all cycles back into a single MASLD cohort dataset
masld_full <- bind_rows(comorb_list)

# 5. Categorize Comorbidity and Lifestyle Groups
cat("\nCategorizing comorbidities and lifestyle variables...\n")

# Calculate calorie tertiles within the MASLD cohort
kcal_tertiles <- quantile(masld_full$Kcal, probs = c(1/3, 2/3), na.rm = TRUE)
cat("Caloric intake tertiles (Day 1): \n")
print(kcal_tertiles)

# Convert categorical variables to clean factors
masld_full <- masld_full %>%
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
    # gender_male: RIAGENDR is factor in baseline (Male/Female)
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
    
    # Mortality Category (1999-2018 cycles only, 2021-2023 set to NA)
    Mortality_Cat = case_when(
      MORTSTAT == 1 ~ "Deceased",
      MORTSTAT == 0 ~ "Alive",
      TRUE ~ NA_character_
    )
  )

# Write patient-level dataset with comorbidity details
write.csv(masld_full, "usfli_nhanes_comorbidities.csv", row.names = FALSE)
cat("Saved detailed patient dataset to 'usfli_nhanes_comorbidities.csv'\n")

# 6. Calculate Stratification Statistics (Raw and Survey-Weighted)
cat("\nCalculating stratification metrics...\n")

# Set up complex survey design for the full cohort
nhanes_design <- svydesign(
  id = ~SDMVPSU, 
  strata = ~SDMVSTRA, 
  weights = ~FastingWeight, 
  data = masld_full, 
  nest = TRUE
)

# Function to get counts and percentages for a categorization variable
get_stratification_stats <- function(design_obj, var_name, category_title) {
  # Formula for the survey calculations
  f <- as.formula(paste0("~", var_name))
  
  # Raw counts
  df_data <- design_obj$variables
  raw_table <- table(df_data[[var_name]], useNA = "no")
  raw_total <- sum(raw_table)
  
  # Coerce character variables to factors with explicit levels to avoid survey package model matrix errors
  design_obj$variables[[var_name]] <- factor(df_data[[var_name]], levels = names(raw_table))
  
  # Survey-weighted estimates
  weighted_prev <- svymean(f, design_obj, na.rm = TRUE)
  weighted_vals <- as.numeric(weighted_prev) * 100
  weighted_se <- as.numeric(SE(weighted_prev)) * 100
  
  # Names mapping
  strata_names <- names(raw_table)
  # Standardize level names returned by svymean (which have variable prefix)
  svy_names <- gsub(paste0("^", var_name), "", names(coef(weighted_prev)))
  
  result_rows <- list()
  for (i in seq_along(strata_names)) {
    sname <- strata_names[i]
    raw_n <- as.numeric(raw_table[sname])
    raw_pct <- (raw_n / raw_total) * 100
    
    # Find matching survey estimation
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

# Run stratification across all groups
obesity_stats <- get_stratification_stats(nhanes_design, "Obesity_Cat", "Obesity/BMI")
diabetes_stats <- get_stratification_stats(nhanes_design, "Diabetes_Cat", "Type 2 DM")
hypertension_stats <- get_stratification_stats(nhanes_design, "Hypertension_Cat", "Hypertension")
dyslipidemia_stats <- get_stratification_stats(nhanes_design, "Dyslipidemia_Cat", "Dyslipidemia")
diet_stats <- get_stratification_stats(nhanes_design, "Diet_Cat", "Diet (Caloric Intake)")
pa_stats <- get_stratification_stats(nhanes_design, "PA_Cat", "Physical Activity")
smoking_stats <- get_stratification_stats(nhanes_design, "Smoking_Cat", "Smoking Status")

# Mortality requires subsetting to eligible cycles (1999-2018)
design_mortality <- subset(nhanes_design, Cycle != "2021-2023" & !is.na(MORTSTAT))
mortality_stats <- get_stratification_stats(design_mortality, "Mortality_Cat", "Mortality (1999-2018 cohorts)")

# Combine all results into a single table
final_stratification <- bind_rows(
  obesity_stats,
  diabetes_stats,
  hypertension_stats,
  dyslipidemia_stats,
  diet_stats,
  pa_stats,
  smoking_stats,
  mortality_stats
)

# 7. Print and Save Summary Results
cat("\n============================================================\n")
cat("STRATIFICATION OF MASLD COHORT (N = 8,673)\n")
cat("============================================================\n")
print(final_stratification)

write.csv(final_stratification, "usfli_comorbidities_stratification.csv", row.names = FALSE)
cat("\nSaved stratification summary to 'usfli_comorbidities_stratification.csv'\n")
cat("Analysis completed successfully!\n")
