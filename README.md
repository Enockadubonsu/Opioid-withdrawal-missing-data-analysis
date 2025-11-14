# Lofexidine for Opioid Withdrawal: Rigorous Reanalysis with Comprehensive Missing Data Sensitivity Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-4.3.1-blue.svg)](https://www.r-project.org/)

## Overview

This repository contains a complete reanalysis of the Yu et al. (2008) phase III randomized controlled trial evaluating lofexidine versus placebo for opioid withdrawal management, with emphasis on rigorous missing data handling methods.

**Key Features:**
- ✅ Linear mixed models using all available data (valid under MAR)
- ✅ Six comprehensive sensitivity analyses
- ✅ Multiple imputation (direct & passive component derivation)
- ✅ MNAR tipping point analysis
- ✅ Fully reproducible R code
- ✅ manuscript

## Study Context

**Original Study:** Yu E, Miotto K, Akerele E, et al. (2008). A phase 3 placebo-controlled, double-blind, multi-site trial of the alpha-2-adrenergic agonist, lofexidine, for opioid withdrawal. *Drug Alcohol Depend*, 97(1-2):158-168.

**Data Source:** National Institute on Drug Abuse (NIDA) Data Share Program  
**Analysis Sample:** N=68 participants (subset of original N=264)

## Key Findings

**Primary Result:** Lofexidine reduced withdrawal severity by **11.84 MHOWS points** (95% CI: -18.79 to -4.90; p=0.001) compared to placebo during Days 4-8.

**Sensitivity Analysis Highlights:**
- Treatment effects remained significant across all 6 sensitivity analyses
- Multiple imputation approaches (direct vs. passive) yielded nearly identical results (difference = 0.28 points, 4%)
- No tipping point identified under MNAR assumptions up to δ=1.5 SD
- Complete-case analysis showed 38% precision loss but similar point estimates

**Missing Data:** 75% by Day 8; 84% dropout (placebo) vs. 66% (lofexidine)

## Repository Structure
```
├── code/               # R analysis scripts (01-07)
├── data/               # Raw data files (CSV)
├── output/             # Generated tables, figures, results
├── manuscript/         # Final report (Word/PDF)
└── environment/        # Session info for reproducibility
```

## Reproducibility

### Software Requirements

- **R version:** 4.3.1
- **Required packages:**
```r
  tidyverse, nlme, lme4, emmeans, mice, 
  broom.mixed, gt, gtsummary, ggpubr, 
  naniar, VIM, gridExtra, cowplot
```

### Running the Analysis

1. **Clone the repository:**
```bash
   git clone https://github.com/YOUR_USERNAME/lofexidine-withdrawal-missing-data-analysis.git
   cd lofexidine-withdrawal-missing-data-analysis
```

2. **Install required packages:**
```r
   install.packages(c("tidyverse", "nlme", "emmeans", "mice", 
                      "broom.mixed", "gt", "gtsummary"))
```

3. **Run scripts sequentially:**
```r
   source("code/01_data_cleaning.R")
   source("code/02_descriptive_stats.R")
   source("code/03_primary_analysis.R")
   source("code/04_sensitivity_analyses.R")
   source("code/05_multiple_imputation.R")
   source("code/06_passive_imputation.R")
   source("code/07_visualization.R")
```

4. **View outputs:**
   - Tables: `output/tables/`
   - Figures: `output/figures/`
   - Manuscript: `manuscript/lofexidine_report.pdf`

## Analytical Methods

### Primary Analysis
- **Model:** Linear mixed-effects model (LMM) with random intercepts and slopes
- **Estimation:** Restricted maximum likelihood (REML)
- **Assumptions:** Missing at random (MAR)
- **Software:** `nlme` package in R

### Sensitivity Analyses

| Analysis | Purpose | Key Finding |
|----------|---------|-------------|
| **SA1:** Alternative covariances | Test robustness to correlation structure | Estimates differ by <3% |
| **SA2:** MI Direct | Alternative MAR implementation | TE = -6.42 (FMI=46%) |
| **SA3:** Complete-case | Quantify efficiency loss | 38% wider CI |
| **SA4:** Covariate adjustment | Test precision gain | Minimal impact (6%) |
| **SA5:** MNAR tipping point | Test MNAR robustness | No tipping point <1.5 SD |
| **SA6:** MI Passive derivation | Component-level imputation | TE = -6.70 (FMI=36%) |

## Key Statistical Concepts

**Missing Data Mechanisms:**
- **MCAR:** Missingness unrelated to observed/unobserved data
- **MAR:** Missingness related to observed data only
- **MNAR:** Missingness related to unobserved data

**Why Multiple Methods?**
- Primary LMM uses all available observed data (efficient under MAR)
- Multiple imputation fills in missing values then analyzes (alternative MAR)
- Tipping point analysis tests MNAR scenarios (robustness check)

**Fraction of Missing Information (FMI):** Proportion of total variance due to missingness
- Direct MI: 46.2%
- Passive MI: 35.9% (more efficient)

## Tables and Figures

### Main Tables
- **Table 1:** Baseline characteristics by treatment arm
- **Table 2:** Missing data pattern by timepoint
- **Table 3:** Treatment effects from primary LMM
- **Table 4:** Comprehensive sensitivity analysis summary

### Main Figures
- **Figure 1:** Participant flow diagram (CONSORT)
- **Figure 2:** Individual MHOWS trajectories (spaghetti plot)
- **Figure 3:** Mean MHOWS trajectories with 95% CI

## Citation

If you use this code or analysis approach, please cite:
```
[Your Name(s)]. (2024). Lofexidine for Opioid Withdrawal: Rigorous Reanalysis 
with Comprehensive Missing Data Sensitivity Analysis. GitHub repository. 
https://github.com/YOUR_USERNAME/lofexidine-withdrawal-missing-data-analysis
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**Data License:** Original data from NIDA Data Share Program. Please review NIDA data sharing policies.

## Authors

- **[Your Name]** - *Primary Analyst* - [Your Email/Website]
- **[Collaborator Names if applicable]**

## Acknowledgments

- Original study investigators: Yu et al. (2008)
- Data source: National Institute on Drug Abuse (NIDA) Data Share Program
- Course instructor: [Professor Name] (if academic project)

## Contact

For questions or collaborations, please open an issue or contact [your email].

**Keywords:** opioid withdrawal, lofexidine, missing data, multiple imputation, sensitivity analysis, longitudinal analysis, linear mixed models, MNAR, tipping point analysis, randomized controlled trial, biostatistics
```

## 🔒 LICENSE FILE

Create a `LICENSE` file with this content (MIT License):
```
MIT License

Copyright (c) 2024 [Your Name]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## 📊 .gitignore FILE

Create a `.gitignore` to exclude unnecessary files:
```
# R files
.Rhistory
.RData
.Rproj.user
*.Rproj

# Output files (optional - include if you want outputs tracked)
# output/
# *.png
# *.pdf
# *.csv

# OS files
.DS_Store
Thumbs.db

# Temporary files
~$*.docx
*.tmp

# Large files
*.RData
*.rda
