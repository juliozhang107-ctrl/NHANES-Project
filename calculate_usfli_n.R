# ==============================================================================
# NHANES MASLD / USFLI > 30 Sample Size and Prevalence Calculator
# ==============================================================================
# This script programmatically downloads NHANES data cycles from 1999-2000 to the
# present (August 2021 - August 2023), merges variables, calculates the U.S.
# Fatty Liver Index (USFLI), and estimates both the raw sample sizes (N) and 
# survey-weighted population sizes for participants with USFLI > 30.
# ==============================================================================

# 1. Load Required Libraries
cat("Loading required libraries...\n")
library(nhanesA)
library(dplyr)
library(survey)
library(ggplot2)

# Set survey option to adjust for strata with lonely PSUs (single PSU)
options(survey.lonely.psu = "adjust")

# 2. Configurations
ADULT_AGE_CUTOFF <- 20  # NHANES standard adult age cutoff (years)
CACHE_DIR <- "nhanes_cache"

# Create cache directory if it doesn't exist to speed up subsequent runs
if (!dir.exists(CACHE_DIR)) {
  dir.create(CACHE_DIR)
}

# 3. Define NHANES Cycles and Table Names
# Note: For 2017-2020, we use the CDC combined pre-pandemic cycle (P_ tables) 
# which is recommended by NCHS to maintain national representativeness.
cycles_info <- list(
  "1999-2000" = list(demo="DEMO",     bmx="BMX",     biopro="LAB18",     glu="LAB10AM", ins="LAB10AM", wt="WTSAF2YR"),
  "2001-2002" = list(demo="DEMO_B",   bmx="BMX_B",   biopro="L40_B",     glu="L10AM_B", ins="L10AM_B", wt="WTSAF2YR"),
  "2003-2004" = list(demo="DEMO_C",   bmx="BMX_C",   biopro="L40_C",     glu="L10AM_C", ins="L10AM_C", wt="WTSAF2YR"),
  "2005-2006" = list(demo="DEMO_D",   bmx="BMX_D",   biopro="BIOPRO_D",  glu="GLU_D",   ins="GLU_D",   wt="WTSAF2YR"),
  "2007-2008" = list(demo="DEMO_E",   bmx="BMX_E",   biopro="BIOPRO_E",  glu="GLU_E",   ins="GLU_E",   wt="WTSAF2YR"),
  "2009-2010" = list(demo="DEMO_F",   bmx="BMX_F",   biopro="BIOPRO_F",  glu="GLU_F",   ins="GLU_F",   wt="WTSAF2YR"),
  "2011-2012" = list(demo="DEMO_G",   bmx="BMX_G",   biopro="BIOPRO_G",  glu="GLU_G",   ins="GLU_G",   wt="WTSAF2YR"),
  "2013-2014" = list(demo="DEMO_H",   bmx="BMX_H",   biopro="BIOPRO_H",  glu="GLU_H",   ins="INS_H",   wt="WTSAF2YR"),
  "2015-2016" = list(demo="DEMO_I",   bmx="BMX_I",   biopro="BIOPRO_I",  glu="GLU_I",   ins="INS_I",   wt="WTSAF2YR"),
  "2017-2018" = list(demo="DEMO_J",   bmx="BMX_J",   biopro="BIOPRO_J",  glu="GLU_J",   ins="INS_J",   wt="WTSAF2YR"),
  "2017-2020 (Pre-Pandemic)" = list(demo="P_DEMO", bmx="P_BMX", biopro="P_BIOPRO", glu="P_GLU",   ins="P_INS",   wt="WTSAFPRP"),
  "2021-2023" = list(demo="DEMO_L",   bmx="BMX_L",   biopro="BIOPRO_L",  glu="GLU_L",   ins="INS_L",   wt="WTSAF2YR")
)

# 4. Helper Function: Load and Cache NHANES Datasets
load_nhanes_cached <- function(table_name, cache_dir = CACHE_DIR) {
  cache_file <- file.path(cache_dir, paste0(table_name, ".rds"))
  
  if (file.exists(cache_file)) {
    cat(paste0("  Loading table ", table_name, " from local cache...\n"))
    return(readRDS(cache_file))
  } else {
    cat(paste0("  Downloading table ", table_name, " from CDC NHANES...\n"))
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

# 5. Process Cycles and Calculate USFLI
summary_list <- list()
all_participants <- data.frame()

cat("\nProcessing NHANES cycles...\n")

for (cycle_name in names(cycles_info)) {
  cat("\n------------------------------------------------------------\n")
  cat("Cycle:", cycle_name, "\n")
  cat("------------------------------------------------------------\n")
  info <- cycles_info[[cycle_name]]
  
  # Load dataframes
  demo <- load_nhanes_cached(info$demo)
  bmx  <- load_nhanes_cached(info$bmx)
  biopro <- load_nhanes_cached(info$biopro)
  glu  <- load_nhanes_cached(info$glu)
  ins  <- load_nhanes_cached(info$ins)
  
  # Skip cycle if any key dataset could not be downloaded
  if (is.null(demo) || is.null(bmx) || is.null(biopro) || is.null(glu) || is.null(ins)) {
    cat("  Error: Missing required datasets for cycle", cycle_name, ". Skipping.\n")
    next
  }
  
  # Select only the columns of interest to simplify merging and prevent duplicate naming conflicts
  demo_clean <- demo %>% select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH1, SDMVPSU, SDMVSTRA)
  bmx_clean <- bmx %>% select(SEQN, BMXWAIST)
  
  # Identify GGT variable name (should be LBXSGTSI)
  ggt_col <- if("LBXSGTSI" %in% colnames(biopro)) "LBXSGTSI" else if("LBXSGT" %in% colnames(biopro)) "LBXSGT" else NULL
  if (is.null(ggt_col)) {
    cat("  Error: GGT variable not found in biochemistry profile. Skipping cycle.\n")
    next
  }
  biopro_clean <- biopro %>% select(SEQN, GGT = !!sym(ggt_col))
  
  # Fasting Glucose & Subsample Weight
  wt_col <- info$wt
  if (!(wt_col %in% colnames(glu))) {
    # If the weight column has a slightly different spelling, search for it
    wt_match <- grep("WTSAF", colnames(glu), value = TRUE)
    if (length(wt_match) > 0) {
      wt_col <- wt_match[1]
    } else {
      cat("  Error: Fasting weight column not found. Skipping cycle.\n")
      next
    }
  }
  glu_clean <- glu %>% select(SEQN, Glucose = LBXGLU, FastingWeight = !!sym(wt_col))
  
  # Fasting Insulin (convert LBXIN to pmol/L by multiplying by 6.0)
  ins_clean <- ins %>% select(SEQN, Insulin_uIU = LBXIN)
  
  # Merge all datasets
  cycle_merged <- demo_clean %>%
    inner_join(bmx_clean, by = "SEQN") %>%
    inner_join(biopro_clean, by = "SEQN") %>%
    inner_join(glu_clean, by = "SEQN") %>%
    inner_join(ins_clean, by = "SEQN")
  
  cat(paste0("  Merged sample size (raw rows): ", nrow(cycle_merged), "\n"))
  
  # Apply adult filter and handle missing values
  # USFLI formula is validated for adults and requires complete data for variables:
  # Age, Race, GGT, Waist Circumference, Insulin, and Glucose.
  cycle_cleaned <- cycle_merged %>%
    filter(
      RIDAGEYR >= ADULT_AGE_CUTOFF,
      !is.na(RIDRETH1),
      !is.na(GGT),
      !is.na(BMXWAIST),
      !is.na(Glucose),
      !is.na(Insulin_uIU),
      !is.na(FastingWeight),
      FastingWeight > 0
    )
  
  raw_total_n <- nrow(cycle_cleaned)
  cat(paste0("  Cleaned adult sample size (with complete data & positive weights): ", raw_total_n, "\n"))
  
  if (raw_total_n == 0) {
    cat("  Error: No valid observations left after cleaning. Skipping cycle.\n")
    next
  }
  
  # Calculate USFLI
  # Race coding: RIDRETH1 (1 = Mexican American, 4 = Non-Hispanic Black)
  cycle_cleaned <- cycle_cleaned %>%
    mutate(
      eth_str = as.character(RIDRETH1),
      NH_Black = ifelse(eth_str == "Non-Hispanic Black" | eth_str == "4", 1, 0),
      Mex_Amer = ifelse(eth_str == "Mexican American" | eth_str == "1", 1, 0),
      Insulin_pmol = Insulin_uIU * 6.0,  # Convert µIU/mL to pmol/L
      
      # Natural logs
      ln_GGT = log(GGT),
      ln_Insulin = log(Insulin_pmol),
      ln_Glucose = log(Glucose),
      
      # USFLI Linear Predictor
      LP = -0.8073 * NH_Black + 
            0.3458 * Mex_Amer + 
            0.0093 * RIDAGEYR + 
            0.6151 * ln_GGT + 
            0.0249 * BMXWAIST + 
            1.1792 * ln_Insulin + 
            0.8242 * ln_Glucose - 
            14.7812,
      
      # USFLI Score
      USFLI = (exp(LP) / (1 + exp(LP))) * 100,
      
      # USFLI > 30 Definition for MASLD
      USFLI_GT_30 = ifelse(USFLI > 30, 1, 0),
      
      # Save cycle name
      Cycle = cycle_name
    )
  
  raw_usfli_gt_30_n <- sum(cycle_cleaned$USFLI_GT_30)
  raw_prev <- (raw_usfli_gt_30_n / raw_total_n) * 100
  
  # Set up complex survey design
  cat("  Calculating survey-weighted statistics...\n")
  nhanes_design <- svydesign(
    id = ~SDMVPSU, 
    strata = ~SDMVSTRA, 
    weights = ~FastingWeight, 
    data = cycle_cleaned, 
    nest = TRUE
  )
  
  # Weighted prevalence estimate
  prev_est <- svymean(~USFLI_GT_30, nhanes_design, na.rm = TRUE)
  weighted_prev <- as.numeric(prev_est) * 100
  weighted_prev_se <- as.numeric(SE(prev_est)) * 100
  
  # Weighted population size estimate
  tot_est <- svytotal(~USFLI_GT_30, nhanes_design, na.rm = TRUE)
  weighted_n <- as.numeric(tot_est)
  weighted_n_se <- as.numeric(SE(tot_est))
  
  # Total population represented
  pop_represented_est <- svytotal(~rep(1, nrow(cycle_cleaned)), nhanes_design)
  total_pop_represented <- as.numeric(pop_represented_est)
  
  cat(paste0("  Raw USFLI > 30 N: ", raw_usfli_gt_30_n, " / ", raw_total_n, " (", round(raw_prev, 2), "%)\n"))
  cat(paste0("  Weighted Prevalence: ", round(weighted_prev, 2), "% (SE: ", round(weighted_prev_se, 2), "%)\n"))
  cat(paste0("  Weighted Population Size N: ", round(weighted_n), " (SE: ", round(weighted_n_se), ")\n"))
  
  # Append to summary list
  summary_list[[cycle_name]] <- data.frame(
    Cycle = cycle_name,
    Raw_Fasting_Total_N = raw_total_n,
    Raw_USFLI_GT_30_N = raw_usfli_gt_30_n,
    Raw_Prevalence_Pct = raw_prev,
    Weighted_Prevalence_Pct = weighted_prev,
    Weighted_Prevalence_SE_Pct = weighted_prev_se,
    Weighted_Population_Size_N = weighted_n,
    Weighted_Population_Size_SE_N = weighted_n_se,
    Total_Population_Represented_N = total_pop_represented,
    stringsAsFactors = FALSE
  )
  
  # Append to participant-level dataframe
  all_participants <- bind_rows(all_participants, cycle_cleaned)
}

# 6. Combine and Save Results
summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

cat("\n============================================================\n")
cat("SUMMARY RESULTS\n")
cat("============================================================\n")
print(summary_df %>% select(Cycle, Raw_Fasting_Total_N, Raw_USFLI_GT_30_N, Raw_Prevalence_Pct, Weighted_Prevalence_Pct))

# Save summary results
write.csv(summary_df, "usfli_summary_results.csv", row.names = FALSE)
cat("\nSaved cycle summary to 'usfli_summary_results.csv'\n")

# Save merged participant-level dataset
write.csv(all_participants %>% select(SEQN, Cycle, RIDAGEYR, RIAGENDR, RIDRETH1, GGT, BMXWAIST, Glucose, Insulin_uIU, Insulin_pmol, USFLI, USFLI_GT_30, FastingWeight, SDMVPSU, SDMVSTRA), 
          "usfli_nhanes_merged.csv", row.names = FALSE)
cat("Saved merged participant-level data to 'usfli_nhanes_merged.csv'\n")

# 7. Generate Trend Plot
cat("\nGenerating prevalence trend plot...\n")
plot_data <- summary_df %>%
  # Filter out the 2017-2020 Pre-Pandemic cycle for plotting if it overlaps with 2017-2018
  # (so we have a clean chronological sequence of 2-year cycles)
  filter(Cycle != "2017-2020 (Pre-Pandemic)")

# Calculate 95% Confidence Intervals for weighted prevalence
plot_data <- plot_data %>%
  mutate(
    CI_Lower = pmax(0, Weighted_Prevalence_Pct - 1.96 * Weighted_Prevalence_SE_Pct),
    CI_Upper = pmin(100, Weighted_Prevalence_Pct + 1.96 * Weighted_Prevalence_SE_Pct),
    # Create an ordered factor for cycle ordering on x-axis
    Cycle_Factor = factor(Cycle, levels = plot_data$Cycle)
  )

p <- ggplot(plot_data, aes(x = Cycle_Factor, y = Weighted_Prevalence_Pct, group = 1)) +
  geom_line(color = "#1a73e8", size = 1.2) +
  geom_point(color = "#1a73e8", size = 3) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "#5f6368", alpha = 0.8) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#e8eaed")
  ) +
  labs(
    title = "Weighted Prevalence of USFLI > 30 (MASLD surrogate) in NHANES Adults (Age >= 20)",
    x = "NHANES Survey Cycle",
    y = "Survey-Weighted Prevalence (%)",
    caption = "Error bars represent 95% confidence intervals. Fasting weights utilized. 2017-2020 combined cycle omitted from plot to prevent double counting."
  )

ggsave("usfli_prevalence_trend.png", plot = p, width = 10, height = 6, dpi = 300)
cat("Saved trend plot to 'usfli_prevalence_trend.png'\n")
cat("\nExecution completed successfully!\n")
