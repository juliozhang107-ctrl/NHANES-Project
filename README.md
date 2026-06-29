# NHANES-Project

================================================================================
NHANES MASLD, MASH, AND LIVER FIBROSIS/CIRRHOSIS COMPARATIVE STUDY (1999-2023)
================================================================================

This repository contains the complete analytical pipeline, processed datasets, 
visual flowcharts, trend plots, and final reports comparing non-invasive 
blood-based indices with gold-standard transient elastography (FibroScan) in 
the U.S. adult population using NHANES data from 1999 to 2023.

Author: Julio Min
Class: Epi III Final Project (Dr. Kang)
Date: June 2026

--------------------------------------------------------------------------------
1. PROJECT OVERVIEW
--------------------------------------------------------------------------------
This project consists of two parallel studies:

A. Primary Study (NHANES 1999-2023):
   - Steatosis (MASLD): U.S. Fatty Liver Index (USFLI) > 30 + CMR Risk Factor
   - MASH: MASLD + Fibrosis-4 (FIB-4) > 2.67
   - Significant Fibrosis: MASLD + AST-to-Platelet Ratio Index (APRI) > 0.5
   - Cirrhosis: MASLD + APRI >= 2.0

B. Gold Standard Study (NHANES 2017-2023):
   - Steatosis (MASLD): Controlled Attenuation Parameter (CAP) >= 285 dB/m + CMR
   - MASH: MASLD + FibroScan-AST (FAST) Score >= 0.35
   - Cirrhosis: MASLD + Median Liver Stiffness Measurement (LSM) > 13.6 kPa

--------------------------------------------------------------------------------
2. DIRECTORY STRUCTURE & FILE REGISTRY
--------------------------------------------------------------------------------

Main Folder Contents:
*   README.txt                      - This project handbook.
*   walkthrough.html                - Interactive web report compiling all results.
*   walkthrough.docx                - Microsoft Word version of the final report.
*   masld_study_flowchart.png       - Cohort flowchart for the Primary Study.
*   gold_standard_flowchart.png     - Cohort flowchart for the Gold Standard Study.
*   usfli_prevalence_trend.png      - Multi-line longitudinal prevalence trend plot.

Analytical Scripts:
*   calculate_all_trends.R          - Calculates cycle-by-cycle weighted rates and 
                                      plots the 3-line trend chart.
*   gold_standard_study.R           - Extracts, merges, cleans, and weights 2017-2023 
                                      transient elastography data.
*   stratify_mash.R                 - Defines and stratifies the primary MASH cohort.
*   stratify_apri.R                 - Defines and stratifies the primary fibrosis cohort.
*   stratify_comorbidities.R        - Stratifies baseline demographics.
*   calculate_usfli_n.R             - Initial USFLI calculator script.

Kaggle Submission Folder (/Kaggle):
*   Kaggle/kaggle_writeup.md        - Formatted Markdown project report for Kaggle.
*   Kaggle/masld_study_flowchart.png- Flowchart attachment for Primary Study.
*   Kaggle/gold_standard_flowchart.png- Flowchart attachment for Gold Standard Study.
*   Kaggle/usfli_prevalence_trend.png- Trend plot attachment.

Output Datasets (CSV):
*   usfli_prevalence_trends_all.csv  - Longitudinal weighted prevalences.
*   gold_standard_stratification.csv - Comorbidity numbers for Gold Standard.
*   usfli_nhanes_gold.csv            - Patient-level transient elastography data.
*   usfli_nhanes_mash.csv            - Patient-level primary MASH dataset.
*   usfli_nhanes_fibrosis.csv        - Patient-level primary fibrosis dataset.
*   usfli_nhanes_merged.csv          - Cleaned general fasting subsample.

--------------------------------------------------------------------------------
3. PREREQUISITES & INSTALLATION
--------------------------------------------------------------------------------
The scripts require R and Python 3.

To install R package dependencies, run this in your R console:
  install.packages(c("survey", "dplyr", "ggplot2", "nhanesA", "readr"))

To install Python dependencies (used for Word report generation), run in shell:
  pip install python-docx pandas matplotlib

--------------------------------------------------------------------------------
4. HOW TO RUN THE ANALYSIS
--------------------------------------------------------------------------------

Step 1: Calculate Prevalence Trends & Generate Figure 3
  Run:
    Rscript calculate_all_trends.R
  This script joins CBC/Biochemistry datasets with general fasting participants,
  computes weighted prevalence cycle-by-cycle, saves "usfli_prevalence_trends_all.csv", 
  and outputs the trend plot "usfli_prevalence_trend.png".

Step 2: Run VCTE Gold-Standard Study & Generate Comorbidity Tables
  Run:
    Rscript gold_standard_study.R
  This script extracts transient elastography data, applies scaled fasting survey 
  weights (2017-2020 scaled by 2/3, 2021-2023 scaled by 1/3), stratifies patients,
  and saves "gold_standard_stratification.csv" and "usfli_nhanes_gold.csv".

Step 3: Regenerate Word Report
  If you edit "walkthrough.md" and want to compile the updated walkthrough.docx, run:
    python3 scratch/convert_to_docx.py

--------------------------------------------------------------------------------
5. SCIENTIFIC REFERENCES
--------------------------------------------------------------------------------
1. Ruhl CE, Everhart JE. Fatty liver indices in the multiethnic United States 
   National Health and Nutrition Examination Survey. Aliment Pharmacol Ther. 
   2015;41(1):65-76. PMID: 25348633
2. Targher G, Valenti L, Byrne CD. Metabolic Dysfunction-Associated Steatotic 
   Liver Disease. The New England Journal of Medicine. 2025.
3. Dreytser E, Blyuss O, Mudrova A, et al. Aspartate Aminotransferase-to-Platelet 
   Ratio Index (APRI) for Staging of Fibrosis in Adults With Chronic Hepatitis C. 
   The Cochrane Database of Systematic Reviews. 2025.
4. Lee BP, Dodge JL, Terrault NA. National prevalence estimates for steatotic 
   liver disease and subclassifications using consensus nomenclature. Hepatology. 
   2024;79(3):666-673.
5. Rinella ME, Lazarus JV, Ratziu V, et al. A multisociety Delphi consensus 
   statement on new fatty liver disease nomenclature. Hepatology. 2023;78(6):1966-1986.
6. Ravaioli F, Dajti E, Mantovani A, et al. Diagnostic Accuracy of FibroScan-AST 
   (FAST) Score for the Non-Invasive Identification of Patients With Fibrotic 
   Non-Alcoholic Steatohepatitis. Gut. 2023.
7. Sterling RK, Duarte-Rojo A, Patel K, et al. AASLD Practice Guideline on 
   Imaging-Based Noninvasive Liver Disease Assessment of Hepatic Fibrosis 
   and Steatosis. Hepatology. 2025.
================================================================================
