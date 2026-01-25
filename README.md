# Time Mastery Project, Sense of Time Battery

This repository hosts the code and materials for a research project validating psychometric tools for assessing **time-related abilities in children**, and for providing a **Shiny app** that automates scoring and standardization of the *Sense of Time* battery.

- Project website: https://www.timemasteryproject.com/  
- Shiny app (battery scoring): https://timemastery.shinyapps.io/TPTM-battery/

---

## Overview

The *Time Mastery Project* focuses on the assessment of children’s temporal processing abilities, combining experimental tasks and questionnaires into a unified psychometric battery.  
This repository contains:

- the full **norming and table-export pipeline** used during development,
- **auditable conversion tables** for raw-to-standard-score mapping,
- and a **Shiny application** for automated scoring, interpretation, and reporting.

A core design principle is that **end users do not need access to GAMLSS models or fitted objects**. All standardization in the Shiny app is performed via explicit lookup tables that can also be used manually by professionals.

---

## Contents of the repository

### 1. Shiny app for scoring and standardization

The Shiny app provides automated scoring for the Sense of Time battery, including:

**Inputs**
- Age (years, months)
- Raw scores for questionnaire-based measures
- Upload of task output files (CSV) for:
  - **Time Reproduction (TR)** from PsychoPy/OpenSesame
  - **Time Discrimination (TD)** from PsychoPy

**Outputs**
- Profile plot with **scaled scores 1–19** (mean 10, SD 3)
- Results table reporting:
  - raw score  
  - aligned z score  
  - scaled score (1–19)
- Optional **Sense of Time Quotient (SoTQ)**, computed only when all required subtests are available

The Shiny app code is located in:

```text
shiny/app.R
```

---

### 2. Norming pipeline and exported conversion tables

Normative models are developed in R using GAMLSS for internal research purposes only.  
These models are then exported into **explicit CSV lookup tables** used by the Shiny app.

Key exported files:

- `scoring/norms_lookup.csv`  
  Main raw → aligned z → scaled score lookup table, indexed by age.

- `scoring/ss_to_raw_intervals.csv`  
  Human-readable intervals allowing manual conversion from scaled scores to raw values.

- `scoring/SoT_totalConversionTable.csv`  
  Conversion table from aggregated scaled scores to the Sense of Time Quotient.

- `scoring/TR_total_norm_params.csv`  
  Age-specific parameters for standardizing the TR total score computed from 11 duration-specific deviations.

Optional interpretive metadata:
- `scoring/SoTQ_CI_params.csv`  
  Parameters for computing confidence intervals around the SoTQ.
- `scoring/SoTQ_profileSpreadThreshold.csv`  
  Threshold for flagging unusually heterogeneous profiles.

---

### 3. Task scoring logic (TR and TD)

Task scoring implemented in the Shiny app mirrors the original experimental scripts.

**Time Discrimination (TD)**
- Reversal points are extracted from the adaptive staircase.
- The first two reversals are discarded.
- The mean of the last six available reversals is computed (or fewer if fewer remain).
- The resulting ratio is a threshold measure, lower values indicate better performance.

**Time Reproduction (TR)**
- Invalid trials are excluded (<100 ms or >36 s).
- Absolute percent deviation from the target duration is computed.
- Deviations are calculated separately for each of the 11 target durations (2–12 s).
- These 11 values form the basis for TR total standardization.

---

## Repository structure

Main directories:

```text
shiny/        # Shiny app code and deployed scoring tables
scoring/      # Exported lookup tables and measure specifications
R/            # Norming, export, sanity-check, and composite scripts
materials/    # Task materials and questionnaires
data/         # Research datasets (not necessarily redistributable)
```

See `repo_tree.txt` for the complete directory tree.

---

## Developer workflow: rebuilding norms and tables

The following scripts reproduce the full norming and export pipeline:

1. Fit normative models (internal use):
```r
source("R/01_build_normPercTabs_full.R")
```

2. Export lookup tables for the Shiny app:
```r
source("R/02_export_norms_to_tables.R")
```

3. Run sanity checks on exported tables:
```r
source("R/03_sanity_check_norms.R")
```

4. Build the SoTQ conversion table and TR-total parameters:
```r
source("R/04_build_SoTQ_table.R")
```

5. Estimate reliability and CI parameters for the SoTQ:
```r
source("R/05_estimate_SoTQ_CI_params.R")
```

---

## Scoring conventions

- All reported profile scores use **scaled scores 1–19** (mean 10, SD 3).
- Scaled scores are derived from aligned z scores as:

```text
SS = clamp(round(10 + 3 * z_aligned), 1, 19)
```

- Directionality is harmonized so that:
  - Higher scores always reflect **better performance**.
  - Error-based or threshold measures (e.g., TD ratio, TR deviation) are aligned accordingly before plotting or aggregation.

---

## Citation

If you use this codebase, the Shiny app, or derived scores in academic work, please cite the relevant *Time Mastery Project* publications (to be added) and acknowledge the project website: https://www.timemasteryproject.com/

---

## License

Creative Commons Attribution–NonCommercial 4.0 International (CC BY-NC 4.0)

---

## Contact

For scientific questions, validation details, collaboration inquiries, or anything, please contact the project maintainers via the Time Mastery Project website: https://www.timemasteryproject.com/contatti
