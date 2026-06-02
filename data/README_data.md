# Data Directory

## What is included here

`misc_data/` contains small crosswalk and classification files committed to the repository:

| File | Description |
|------|-------------|
| `occegp.csv` | Occupation–EGP class crosswalk |
| `Autor_Dorn_occ_class_scheme.dta` | Autor–Dorn occupational classification |
| `occ1990_titles.dta` | 1990 Census occupation titles |
| `occ1990_occ1990dd.dta` | 1990 occupation harmonized crosswalk |
| `occ2010_occ1990dd.dta` | 2010→1990 occupation crosswalk |
| `soc2010_to_occ2010.xlsx` | SOC 2010 to OCC 2010 crosswalk |
| `soc_2000_to_2010_crosswalk.xls` | SOC 2000→2010 crosswalk |
| `socoes98.xls` | OES 1998 occupation classification |
| `soc_structure_2010.xls` | SOC 2010 hierarchy |
| `1990CZ_LMA.txt` | 1990 commuting zone labor market areas |
| `cz1990_shapefile/` | Shapefile for 1990 commuting zones |
| `PUMA_CZ/` | PUMA-to-CZ crosswalks (1980, 1990, 2000, 2010) |
| `macro_meso_micro_occ/` | Grusky occupational class scheme |

`intermediate/` contains small precomputed results:

| File | Description |
|------|-------------|
| `placebo_results.RData` | Monte Carlo placebo permutation results (SI S3) |

## What is NOT included (too large for git)

The following large files must be obtained separately and placed in your local `data_dir` (set in each script's path configuration block at the top):

### Cleaned analysis datasets (~1 GB each; built by `code/00 prelude.R`)
- `cps_cleaned_2.RData` — CPS individual-level panel, 1994–2022
- `acs_cleaned_2.RData` — ACS/Census individual-level panel, 1980–2021
- `acs_tract.RData` — ACS tract-level occupational counts, 2007–2021

### Raw IPUMS extracts (required only to run `00 prelude.R`)
- `raw_data/cps_00006.dta` — IPUMS CPS extract
- `raw_data/usa_00063.xml` + `.dat.gz` — IPUMS ACS extract
- `raw_data/OEWS_97-21/` — BLS OEWS MSA files, 1999–2021

### Intermediate LOO computation results (built by analysis scripts)
- `intermediate/s4_loo_results.rds` — CPS LOO location results
- `intermediate/oews_occ_loo_m1.rds` — OEWS occupation LOO results
- `intermediate/onet_soc39.rds` — O*NET SOC-39 scores

### Not redistributable
- `msa_shapefile/` — NHGIS MSA shapefiles (922 MB; download from nhgis.org)
- `soc_code_strat.dta` — SOC stratification file (472 MB)
- `CPI-U-RS_1978_to_2020.dta` — BLS CPI-U-RS deflator (download from bls.gov)

See `README.md` for full download and setup instructions.
