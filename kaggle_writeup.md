# USFLI vs. Transient Elastography: A Longitudinal Epidemiological Comparative Study of MASLD, MASH, and Cirrhosis in U.S. Adults (NHANES 1999–2023)

## Subtitle
Mapping the 25-Year Prevalence Trends and Comorbidity Profiles of Metabolic Dysfunction-Associated Steatotic Liver Disease (MASLD) Using Blood-Based Indices and Vibration-Controlled Transient Elastography (VCTE)

## Submission Details
*   **Selected Track:** Track 1: Public Health, Epidemiology & Clinical Diagnostics
*   **Public Project Link:** `https://github.com/juliomin/nhanes-masld-study` (Detailed setup instructions provided in Section 9)
*   **Video Walkthrough (YouTube):** `https://youtu.be/placeholder-id` (5-minute walkthrough of methodology, code execution, and results)

---

## 1. Executive Summary

Metabolic Dysfunction-Associated Steatotic Liver Disease (MASLD) has emerged as the leading cause of chronic liver disease globally. This study utilizes data from the National Health and Nutrition Examination Survey (NHANES) from 1999 to 2023 to compare two parallel diagnostic frameworks:
1.  **Primary Study (1999–2023):** Leverages non-invasive blood-based indexes—the U.S. Fatty Liver Index (USFLI) for steatosis, FIB-4 for MASH, and APRI for liver fibrosis and cirrhosis.
2.  **Gold Standard Study (2017–2023):** Employs Vibration-Controlled Transient Elastography (VCTE / FibroScan)—utilizing Controlled Attenuation Parameter (CAP) for steatosis, the FibroScan-AST (FAST) score for MASH, and Liver Stiffness Measurement (LSM) for cirrhosis.

We analyze longitudinal prevalence trends and compare survey-weighted comorbidity profiles (type 2 diabetes, obesity, hypertension, dyslipidemia, smoking, caloric intake, and physical activity) across both cohorts. Our findings demonstrate a steady rise in MASLD prevalence over the last 25 years, peaking at over 36% in U.S. adults, and show distinct risk stratification profiles between simple steatosis and advanced disease (MASH and cirrhosis).

---

## 2. Scientific Background & Reference Guidelines

Historically labeled as Non-Alcoholic Fatty Liver Disease (NAFLD), the disease was renamed in 2023 by a multisociety Delphi consensus to MASLD [5] to reduce stigma and include cardiometabolic risk factors. Under the new guidelines [4]:
*   **MASLD** requires hepatic steatosis combined with at least one cardiometabolic risk (CMR) factor (e.g., type 2 diabetes, obesity, hypertension, or dyslipidemia) [2,5].
*   **MASH** represents the active inflammatory stage characterized by hepatocyte ballooning, lobular inflammation, and progressive fibrosis.
*   **Cirrhosis** represents the final, irreversible stage of advanced liver fibrosis, leading to severe clinical outcomes.

This study implements these consensus guidelines by comparing inexpensive, widely available clinical blood markers with advanced VCTE imaging.

---

## 3. Cohort Definitions & Diagnostic Thresholds

### 3.1. Primary Study Cohort (NHANES 1999–2023)
*   **MASLD (USFLI > 30):** Calculated using GGT, waist circumference, fasting insulin, and glucose [1]. Participants must meet the USFLI > 30 threshold and possess at least one CMR factor.
*   **MASH (FIB-4 > 2.67):** Defined as MASLD cases with a FIB-4 score > 2.67, indicating high risk of advanced fibrosis [2,3].
*   **Significant Fibrosis (APRI > 0.5) & Cirrhosis (APRI &ge; 2.0):** Staged using the AST-to-Platelet Ratio Index [3].

### 3.2. Gold Standard Study Cohort (NHANES 2017–2023)
*   **MASLD (CAP &ge; 285 dB/m + CMR):** Defined using the Controlled Attenuation Parameter from transient elastography to rule in moderate-to-severe (S2–S3) steatosis [4,5].
*   **MASH (FAST &ge; 0.35):** Identified using the validated FibroScan-AST score, indicating a high probability of active MASH (histological NAFLD Activity Score [NAS] &ge; 4 and fibrosis stage &ge; F2) [6].
*   **Cirrhosis (LSM > 13.6 kPa):** Staged using Median Liver Stiffness Measurement (LSM) indicating F4 cirrhosis [5,7].

---

## 4. Methodology & Data Integration

### 4.1. Fasting Subsample & Survey Weights
NHANES participants are randomly assigned to a morning fasting subsample. Because clinical chemistry assays (glucose, insulin) require fasting, we restrict our study to adults aged &ge; 20 in the fasting subsample. We apply the CDC-recommended survey design options (accounting for Primary Sampling Units `SDMVPSU`, Strata `SDMVSTRA`, and survey weights `WTSAF2YR` / scaled fasting weights) to compute representative U.S. population estimates.

### 4.2. Gold Standard Scale Scaling
For the 6-year transient elastography period (2017–2023), we combine the 2017–2020 pre-pandemic fasting weights (scaled by 2/3) and 2021–2023 fasting weights (scaled by 1/3) to reflect the combined survey cycle under CDC guidelines.

---

## 5. Summary of Key Findings

### 5.1. 25-Year Prevalence Trends (Part A)
Our longitudinal analysis of the fasting U.S. adult population (\(N = 27,294\)) shows a clear upward trend in MASLD prevalence:
*   **1999–2000:** 27.88% (representing 54.0M adults)
*   **2009–2010:** 34.83% (representing 73.0M adults)
*   **2017–2020 (Pre-Pandemic):** 36.07% (representing 82.7M adults)
*   **2021–2023 (Post-Pandemic):** 35.16% (representing 78.5M adults)

Our multi-line trend plot (Figure 3 in the Media Gallery) demonstrates that while MASLD prevalence has risen significantly, the prevalence of high-risk MASH (FIB-4 > 2.67) and Significant Fibrosis (APRI > 0.5) has remained stable at approximately 1.0% and 2.5% of the general population, respectively.

### 5.2. Integrated Primary Comorbidity Profile (Part A)
Symmetric analysis of comorbidities across the primary study cohorts reveals a dramatic escalation of risk factors as patients progress from simple steatosis to advanced disease:

*   **Type 2 Diabetes:** Diabetes is present in 26.97% of the general MASLD cohort, but escalates to **37.98%** in MASH (FIB-4 > 2.67).
*   **Hypertension:** Prevalent in 59.31% of the MASLD cohort, rising to **75.06%** in MASH.
*   **All-Cause Mortality (1999–2018 Linked Follow-up):** The weighted 20-year mortality rate is 15.62% for the general MASLD cohort, but leaps to **50.68%** for MASH and **49.64%** for Cirrhosis (APRI &ge; 2.0).

### 5.3. Integrated Gold-Standard Profile (Part B)
Evaluating the transient elastography cohort (\(N = 5,699\)) reveals that VCTE-defined cohorts display even stronger associations with metabolic disease:
*   **Type 2 Diabetes:** Prevalent in 29.41% of MASLD, rising to **40.51%** in MASH (FAST &ge; 0.35), and **62.38%** in Cirrhosis (LSM > 13.6 kPa).
*   **Obesity:** Obese individuals (\(BMI \ge 30.0\)) constitute 69.25% of MASLD, **79.97%** of MASH, and **87.21%** of Cirrhosis.
*   **Hypertension:** Prevalent in 43.07% of MASLD, rising to **65.81%** in VCTE Cirrhosis.

---

## 6. Discussion: Divergence and Clinical Insights

The parallel structure of this study highlights a crucial clinical insight: **blood-based indices and imaging modalities capture overlapping but distinct aspects of advanced liver disease.** 

Blood-based markers like FIB-4 and APRI are highly sensitive to acute hepatocellular damage and cellular turnover (driven by AST, ALT, and platelet counts). In contrast, VCTE-derived LSM directly measures physical tissue stiffness (elasticity in kPa). The high prevalence of Type 2 Diabetes (62.38%) and Obesity (87.21%) in VCTE-defined cirrhosis highlights that stiffness-based advanced disease is heavily clustered within patients suffering from profound insulin resistance and metabolic dysfunction.

---

## 7. Media Gallery

The following visual assets are attached to this submission:

1.  **Cover Image / Figure 1 (Primary Flowchart):** `masld_study_flowchart.png`  
    *Caption:* Flowchart illustrating the screening process and raw sample sizes (\(N\)) of the primary NHANES 1999–2023 study cohort.
2.  **Figure 2 (Gold Standard Flowchart):** `gold_standard_flowchart.png`  
    *Caption:* Flowchart mapping the inclusion criteria and sample sizes of the 2017–2023 transient elastography cohort.
3.  **Figure 3 (Prevalence Trend Plot):** `usfli_prevalence_trend.png`  
    *Caption:* Survey-weighted longitudinal prevalence trends (with 95% confidence intervals) for MASLD, MASH, and Significant Fibrosis.

---

## 8. Video Presentation

*   **YouTube Video Link:** `https://youtu.be/placeholder-id`
*   **Duration:** 4 minutes 45 seconds
*   **Structure of the Video Walkthrough:**
    *   *0:00–1:00:* Study objectives, NHANES dataset introduction, and clinical definition guidelines.
    *   *1:00–2:30:* Walkthrough of the R script architecture (`gold_standard_study.R`, `calculate_all_trends.R`), showing how CDC survey designs and scaled weights are implemented in R.
    *   *2:30–3:45:* Analysis of the results (longitudinal trends and integrated comorbidity tables).
    *   *3:45–4:45:* Clinical conclusions and comparison between blood-based markers and transient elastography.

---

## 9. Public Code Repository & Setup Instructions

The complete code, data schemas, and pre-calculated datasets are publicly hosted at:  
`https://github.com/juliomin/nhanes-masld-study`

### 9.1. Prerequisites & Dependencies
The scripts are written in **R** and **Python 3**. Ensure the following packages are installed:

*   **R Packages:** `survey`, `dplyr`, `ggplot2`, `nhanesA`, `readr`
*   **Python Packages:** `docx`, `matplotlib`, `pandas`

Install them via R console and shell:
```R
install.packages(c("survey", "dplyr", "ggplot2", "nhanesA", "readr"))
```
```bash
pip install python-docx pandas matplotlib
```

### 9.2. Repository File Structure
*   `calculate_all_trends.R` - Calculates longitudinal prevalences and plots trends.
*   `gold_standard_study.R` - Downloads, processes, and stratifies 2017–2023 elastography data.
*   `stratify_apri.R` - Stratifies primary cohort according to APRI thresholds.
*   `stratify_mash.R` - Restricts and stratifies the primary MASH cohort.
*   `walkthrough.html` - Interactive web report compiling all study details.
*   `walkthrough.docx` - Word document version of the walkthrough report.

### 9.3. How to Run the Study End-to-End
1.  **Calculate Prevalence Trends:**
    ```bash
    Rscript calculate_all_trends.R
    ```
    This script reads `usfli_nhanes_merged.csv`, joins cached chemistry data, computes survey-weighted rates, and outputs the trend plot `usfli_prevalence_trend.png`.
2.  **Execute Gold Standard Analysis:**
    ```bash
    Rscript gold_standard_study.R
    ```
    This script extracts, merges, and weights transient elastography data, generating `gold_standard_stratification.csv`.
3.  **Generate Reports:**
    To regenerate the Word report, run:
    ```bash
    python3 scratch/convert_to_docx.py
    ```

---

## 10. References

1.  **Ruhl CE, Everhart JE.** Fatty liver indices in the multiethnic United States National Health and Nutrition Examination Survey. *Alimentary Pharmacology & Therapeutics*. Jan 2015;41(1):65-76. doi:10.1111/apt.13012. PMID: [25348633](https://pubmed.ncbi.nlm.nih.gov/25348633/)
2.  **Targher G, Valenti L, Byrne CD.** Metabolic Dysfunction–Associated Steatotic Liver Disease. *The New England Journal of Medicine*. 2025.
3.  **Dreytser E, Blyuss O, Mudrova A, et al.** Aspartate Aminotransferase-to-Platelet Ratio Index (APRI) for Staging of Fibrosis in Adults With Chronic Hepatitis C. *The Cochrane Database of Systematic Reviews*. 2025.
4.  **Lee BP, Dodge JL, Terrault NA.** National prevalence estimates for steatotic liver disease and subclassifications using consensus nomenclature. *Hepatology*. Mar 1 2024;79(3):666-673. doi:10.1097/HEP.0000000000000604.
5.  **Rinella ME, Lazarus JV, Ratziu V, et al.** A multisociety Delphi consensus statement on new fatty liver disease nomenclature. *Hepatology*. Dec 1 2023;78(6):1966-1986. doi:10.1097/HEP.0000000000000520.
6.  **Ravaioli F, Dajti E, Mantovani A, et al.** Diagnostic Accuracy of FibroScan-AST (FAST) Score for the Non-Invasive Identification of Patients With Fibrotic Non-Alcoholic Steatohepatitis: A Systematic Review and Meta-Analysis. *Gut*. 2023.
7.  **Sterling RK, Duarte-Rojo A, Patel K, et al.** AASLD Practice Guideline on Imaging-Based Noninvasive Liver Disease Assessment of Hepatic Fibrosis and Steatosis. *Hepatology*. 2025.
