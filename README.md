# Replication Package

---

## Overview

This repository contains all code, small data files, and output figures needed to replicate the analysis. The main datasets (CPS, ACS, OEWS) are too large to host here; instructions for obtaining them are below. Once you have the data, you run the seven scripts in `code/` in order.

---

## Repository Structure

```
GeoOccGit/
├── code/
│   ├── 00 prelude.R          — Raw data cleaning (run once)
│   ├── 01 data.R             — Data loading; sourced by every analysis script
│   ├── 02 functions.R        — Core helper functions; sourced by every analysis script
│   ├── 03 ACS Analysis.R     — ACS/Census cosine similarity (Figs 2–5, SI S2, S4–S6)
│   ├── 04 OEWS Analysis.R    — OEWS cosine similarity (Figs 2–3, SI S4)
│   ├── 05 CPS Analysis.R     — CPS cosine similarity (SI S2, S4, S5)
│   └── 06 SI Appendix.R      — Robustness checks (SI S1, S3)
├── data/
│   ├── misc_data/            — Crosswalk files (included in repo; see data/README_data.md)
│   └── intermediate/         — Small precomputed results (included in repo)
├── figures/
│   ├── main/                 — Main paper figures (output)
│   └── si/                   — SI appendix figures (output)
│       ├── Cosine Similarity/
│       └── SVD Basics/
├── tex/
│   ├── si_appendix.tex       — SI appendix LaTeX source
│   └── soc11_crosswalk_table.tex
└── README.md
```

---

## Setup: Two Paths to Configure

At the top of each analysis script you will find a short configuration block:

```r
repo_root <- "~/Dropbox/GeoOcc/GeoOccGit"   # path to this repository
data_dir  <- "~/Dropbox/GeoOcc/Analysis/data" # path to large data files (not in repo)
```

Change `repo_root` to wherever you cloned this repository. Change `data_dir` to wherever you store the large data files described below.

---

## Data You Need to Obtain

### 1. IPUMS microdata (requires free registration at ipums.org)

- **CPS**: An IPUMS CPS extract. The variables needed are: `year`, `month`, `cpsid`, `cpsidp`, `serial`, `pernum`, `wtfinl`, `statefip`, `metarea`, `occ1990`, `ind1990`, `egp`, `da`, `age`, `sex`, `race`, `hispan`, `bpl`, `educ`, `empstat`, `classwkr`, `uhrsworkt`, `earnweek`, `labforce`. Save as `cps_cleaned_2.RData` (or run `00 prelude.R`).

- **ACS/Census**: An IPUMS USA extract covering 1980–2021 decennial and ACS. Variables: `year`, `datanum`, `serial`, `pernum`, `perwt`, `statefip`, `puma`, `czone` (commuting zone from IPUMS crosswalk), `occ1990`, `egp`, `da`, `age`, `sex`, `race`, `hispan`, `bpl`, `educ`, `uhrswork`, `incwage`, `afactor`. Save as `acs_cleaned_2.RData`.

- **ACS Tract**: ACS 5-year tract-level occupational counts, 2007–2021, downloaded via the Census API using `tidycensus`. Save as `acs_tract.RData`. This is produced by Section 5 of `03 ACS Analysis.R`.

### 2. BLS Occupational Employment and Wage Statistics (OEWS)

Download MSA-level OEWS files for 1999–2021 from the [BLS website](https://www.bls.gov/oes/tables.htm). The script expects them in `data_dir/raw_data/OEWS_97-21/` with subdirectory names like `oes99ma/`, `oesm03ma/`, etc. The loading code in `01 data.R` handles the varying column formats across vintages.

### 3. CPI deflator (optional)

Download the BLS CPI-U-RS annual research series from [bls.gov](https://www.bls.gov/cpi/research-series/home.htm) and save as `data/misc_data/CPI-U-RS_1978_to_2020.dta`. If this file is absent, the script falls back to nominal wages with a warning.

### 4. MSA shapefiles (for OEWS geographic standardization only)

Download NHGIS CBSA/PMSA shapefiles if you need to regenerate the OEWS metropolitan-area crosswalk. These are used only by the spatial crosswalk step in `04 OEWS Analysis.R` and are not needed if you use the precomputed crosswalk already embedded in the OEWS loading code.

---

## Workflow

Run the scripts in this order. Each one sources `01 data.R` and `02 functions.R` automatically.

### Step 0 — Clean raw data (run once)

```r
source("code/00 prelude.R")
```

Reads raw IPUMS extracts, applies crosswalks and variable recoding, and saves `cps_cleaned_2.RData` and `acs_cleaned_2.RData` to `data_dir`. This takes about 30–60 minutes and only needs to run once.

### Step 1 — ACS analysis (`03 ACS Analysis.R`)

Produces the main ACS/Census cosine similarity figures and tables. Contains six sections:

| Section | Output | Paper reference |
|---------|--------|-----------------|
| 1 | `ACS_DA_PM_base.png/tiff` | Fig. 2, SI S2 |
| 2 | Causal table (printed) | Table 1 |
| 3 | `ACS_DA_occ_heterogeneity.png`, `ACS_DA_occ_ols.png` | Fig. 4 lower, Fig. 5 |
| 4 | `ACS_EGP_PM_base.png` | SI S2 |
| 5 | `acs_M_base_3_tract_*.png` | SI S6 |

Runtime: Section 3 (occupation LOO) is the most computationally intensive, but is made feasible by computing the SVD embedding once per year across all occupations and then dropping each occupation only at the centroid step — avoiding a full SVD recomputation per occupation. Expected runtime is under 10 minutes on a modern laptop.

### Step 2 — OEWS analysis (`04 OEWS Analysis.R`)

Produces the OEWS cosine similarity figures and the OEWS LOO location scatter plot.

| Section | Output | Paper reference |
|---------|--------|-----------------|
| 1 | `OEWS_M_base_3.png/tiff` | Fig. 2–3 |
| 2 | `LOO_location_OEWS.png` | SI S4 |

### Step 3 — CPS analysis (`05 CPS Analysis.R`)

Produces CPS cosine similarity and occupation heterogeneity figures.

| Section | Output | Paper reference |
|---------|--------|-----------------|
| 1 | `CPS_DA_PM_base.png` | SI S2 |
| 2 | `CPS_EGP_PM_base.png`, `CPS_EGP_DA_combined.png` | SI S2 |
| 3 | `LOO_location_CPS.png` | SI S4 |
| 4 | `CPS_occ_heterogeneity.png`, `CPS_occ_ols.png` | SI S5 |

### Step 4 — SI robustness checks (`06 SI Appendix.R`)

Produces SI S1 (raw vs SVD comparison), SI S3 (placebo test and exclusion robustness).

| Section | Output | Paper reference |
|---------|--------|-----------------|
| SI S1 | `figures/si/SVD Basics/raw_vs_svd_cosine*.png` | SI S1 |
| SI S3 | `placebo_combined.png`, `excl_pmpc_combined.png` | SI S3 |
| SI S5 | `OEWS_occ_ols.png` | SI S5 (OEWS only) |

The SI S2 k-robustness figures (`ACS_DA_PM_base_k5.png`, etc.) require re-running Sections 1–3 of the respective analysis scripts with `nu=5` or `nu=15`.

---

## Precomputed intermediate results

`data/intermediate/placebo_results.RData` contains the Monte Carlo placebo permutation results used for SI S3 (100 permutations × 3 datasets). This is included in the repository because the computation takes several hours. If you want to regenerate it from scratch, set `file.exists(placebo_rdata)` to `FALSE` in `06 SI Appendix.R` and allow it to run the compute step.

---

## Software

All analysis is in R. The scripts were developed and tested on R 4.3–4.4. Key packages:

- **Data**: `dplyr`, `tidyr`, `haven`, `readxl`, `zoo`, `sjlabelled`
- **Analysis**: `lfe`, `fixest`, `AER`, `questionr`, `lsa`
- **Spatial**: `sf`, `nngeo`
- **Visualization**: `ggplot2`, `ggrepel`, `gridExtra`, `patchwork`, `scales`
- **Tract SVD**: `irlba`, `tidycensus`

Install all packages with:

```r
install.packages(c("dplyr","tidyr","haven","readxl","zoo","sjlabelled",
                   "lfe","fixest","AER","questionr","lsa",
                   "sf","nngeo","ggplot2","ggrepel","gridExtra",
                   "patchwork","scales","irlba","tidycensus"))
```

---

## Questions

Contact Wenhao Jiang at wenhao.jiang@duke.edu.
