# ==============================================================================
# NHANES MASLD, MASH, and Fibrosis Prevalence Trends (1999-2023)
# ==============================================================================
# This script calculates the longitudinal prevalence trends of MASLD, MASH,
# and Significant Fibrosis from 1999 to 2023, and plots them on a single chart.
# ==============================================================================

# 1. Load Required Libraries
cat("Loading required libraries...\n")
library(dplyr)
library(survey)
library(ggplot2)

# Set survey option to adjust for lonely PSUs
options(survey.lonely.psu = "adjust")

# Configurations
CACHE_DIR <- "nhanes_cache"

# Helper Function: Load and Cache NHANES Datasets
load_nhanes_cached <- function(table_name, cache_dir = CACHE_DIR) {
  cache_file <- file.path(cache_dir, paste0(table_name, ".rds"))
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  } else {
    cat(paste0("  Warning: table ", table_name, " not found in cache.\n"))
    return(NULL)
  }
}

# 2. Define Mappings for CBC and Biochemistry tables
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

# 3. Load General Fasting Cohort Dataset
cat("Loading general fasting cohort dataset...\n")
base_data <- read.csv("usfli_nhanes_merged.csv")
cat("Loaded ", nrow(base_data), " fasting participants.\n")

# 4. Fetch and Process Platelets, AST, and ALT for each Cycle
lab_list <- list()

for (cyc_name in names(cbc_cycles)) {
  cat("Processing lab data for cycle:", cyc_name, "...\n")
  cyc_base <- base_data %>% filter(Cycle == cyc_name)
  if (nrow(cyc_base) == 0) next
  
  cbc_table_name <- cbc_cycles[[cyc_name]]
  biopro_table_name <- biopro_cycles[[cyc_name]]
  
  cbc <- load_nhanes_cached(cbc_table_name)
  biopro <- load_nhanes_cached(biopro_table_name)
  
  cbc_clean <- if (!is.null(cbc) && "LBXPLTSI" %in% colnames(cbc)) cbc %>% select(SEQN, Platelets = LBXPLTSI) else data.frame(SEQN=integer(), Platelets=numeric())
  
  biopro_clean <- if (!is.null(biopro) && "LBXSASSI" %in% colnames(biopro) && "LBXSATSI" %in% colnames(biopro)) {
    biopro %>% select(SEQN, AST = LBXSASSI, ALT = LBXSATSI)
  } else {
    data.frame(SEQN=integer(), AST=numeric(), ALT=numeric())
  }
  
  cyc_merged <- cyc_base %>%
    left_join(cbc_clean, by = "SEQN") %>%
    left_join(biopro_clean, by = "SEQN")
  
  lab_list[[cyc_name]] <- cyc_merged
}

all_data <- bind_rows(lab_list)

# 5. Calculate FIB-4 and APRI
cat("\nCalculating FIB-4 and APRI indexes...\n")
all_data <- all_data %>%
  mutate(
    FIB4 = (RIDAGEYR * AST) / (Platelets * sqrt(ALT)),
    APRI = (AST / 40) * 100 / Platelets,
    
    # Binary indicators (with fallback to 0 if NA or missing)
    MASLD = ifelse(USFLI_GT_30 == 1, 1, 0),
    MASH = ifelse(USFLI_GT_30 == 1 & !is.na(FIB4) & FIB4 > 2.67, 1, 0),
    Fibrosis = ifelse(USFLI_GT_30 == 1 & !is.na(APRI) & APRI > 0.5, 1, 0)
  )

# 6. Calculate Survey-Weighted Prevalence Trend Cycle-by-Cycle
summary_list <- list()

for (cyc_name in unique(all_data$Cycle)) {
  cat("Calculating survey prevalence for:", cyc_name, "...\n")
  cyc_data <- all_data %>% filter(Cycle == cyc_name)
  
  # Set up survey design
  design <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~FastingWeight,
    data = cyc_data,
    nest = TRUE
  )
  
  # MASLD prevalence
  prev_masld <- svymean(~MASLD, design, na.rm = TRUE)
  # MASH prevalence
  prev_mash <- svymean(~MASH, design, na.rm = TRUE)
  # Fibrosis prevalence
  prev_fib <- svymean(~Fibrosis, design, na.rm = TRUE)
  
  summary_list[[cyc_name]] <- data.frame(
    Cycle = cyc_name,
    MASLD_Prev = as.numeric(prev_masld) * 100,
    MASLD_SE = as.numeric(SE(prev_masld)) * 100,
    MASH_Prev = as.numeric(prev_mash) * 100,
    MASH_SE = as.numeric(SE(prev_mash)) * 100,
    Fibrosis_Prev = as.numeric(prev_fib) * 100,
    Fibrosis_SE = as.numeric(SE(prev_fib)) * 100,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

# Save trend summary CSV
write.csv(summary_df, "usfli_prevalence_trends_all.csv", row.names = FALSE)
cat("Saved detailed prevalence trends to 'usfli_prevalence_trends_all.csv'\n")

# 7. Generate Multi-Line Trend Plot
cat("\nGenerating multi-line prevalence trend plot...\n")
plot_data <- summary_df %>%
  filter(Cycle != "2017-2020 (Pre-Pandemic)")

# Reorganize data into long format for ggplot
plot_long <- bind_rows(
  plot_data %>% select(Cycle, Prev = MASLD_Prev, SE = MASLD_SE) %>% mutate(Condition = "MASLD (USFLI > 30)"),
  plot_data %>% select(Cycle, Prev = MASH_Prev, SE = MASH_SE) %>% mutate(Condition = "MASH (FIB-4 > 2.67)"),
  plot_data %>% select(Cycle, Prev = Fibrosis_Prev, SE = Fibrosis_SE) %>% mutate(Condition = "Sig. Fibrosis (APRI > 0.5)")
)

# Calculate 95% CI
plot_long <- plot_long %>%
  mutate(
    CI_Lower = pmax(0, Prev - 1.96 * SE),
    CI_Upper = pmin(100, Prev + 1.96 * SE),
    Cycle_Factor = factor(Cycle, levels = plot_data$Cycle)
  )

# Plot
p <- ggplot(plot_long, aes(x = Cycle_Factor, y = Prev, group = Condition, color = Condition)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, alpha = 0.7, size = 0.8) +
  scale_color_manual(values = c("MASLD (USFLI > 30)" = "#1a73e8", "MASH (FIB-4 > 2.67)" = "#d93025", "Sig. Fibrosis (APRI > 0.5)" = "#188038")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#e8eaed")
  ) +
  labs(
    title = "Weighted Prevalence Trends of MASLD, MASH, & Fibrosis in NHANES (1999-2023)",
    x = "NHANES Survey Cycle",
    y = "Survey-Weighted Prevalence (%) in U.S. Adults",
    color = "Cohort Definition",
    caption = "Error bars represent 95% confidence intervals. Fasting weights utilized. 2017-2020 combined cycle omitted from plot to prevent double counting."
  )

ggsave("usfli_prevalence_trend.png", plot = p, width = 11, height = 6.5, dpi = 300)
cat("Saved updated trend plot to 'usfli_prevalence_trend.png'\n")
cat("\nExecution completed successfully!\n")
