## 01 data.R
## Loads all cleaned datasets and shared objects used by the analysis scripts.
## Sourced at the top of each analysis script via source(file.path(code_dir, "01 data.R")).
##
## Prerequisites (see README.md for download instructions):
##   - 00 prelude.R must be run once to produce cps_cleaned_2.RData and acs_cleaned_2.RData
##   - CPI deflator: data/misc_data/CPI-U-RS_1978_to_2020.dta
##     (BLS CPI-U-RS series: https://www.bls.gov/cpi/research-series/home.htm)
##
## Objects produced (available to sourcing scripts):
##   occegp, occda         — occupation classification crosswalks
##   covariate             — occupation-level baseline covariates (≤1990)
##   cz_density            — CZ-level density and demographics (≤1990)
##   cps                   — CPS individual-level panel, 1994–2022
##   acs                   — ACS/Census individual-level panel, 1980–2021
##   oews                  — OEWS MSA-occupation-year panel, 1999–2021
##   onet_soc39            — O*NET scores for SOC 39 (personal care) occupations

## ── Path configuration ────────────────────────────────────────────────────────
## repo_root: root of this repository (adjust to your local clone path)
## data_dir:  directory containing large cleaned data files (NOT in repository)
##            See README.md and data/README_data.md for download instructions.

repo_root <- "~/Dropbox/GeoOcc/GeoOccGit"
data_dir  <- "~/Dropbox/GeoOcc/Analysis/data"
misc_dir  <- file.path(repo_root, "data", "misc_data")

## ── Packages ────────────────────────────────────────────────────────────────
library(dplyr)
library(haven)
library(sjlabelled)
library(ggplot2)
library(zoo)
library(readxl)
library(sf)

## ── Crosswalks ──────────────────────────────────────────────────────────────
occegp <- read.csv(file.path(misc_dir, "occegp.csv")) %>% arrange(occ1990)
occegp[which(occegp$egp %in% c("IIIa","IIIb")), "egp"] <- "III"
occegp[which(occegp$egp %in% c("VIIa","VIIb")), "egp"] <- "VII"

occgw <- read_dta(file.path(misc_dir, "macro_meso_micro_occ/Grusky_occ_class_scheme.dta")) %>%
  dplyr::select(-occ1990_labels) %>%
  rename(gw_macro = macro_adj, gw_meso = meso_adj, gw_micro = micro_adj)

occda <- read_dta(file.path(misc_dir, "Autor_Dorn_occ_class_scheme.dta")) %>%
  rename(occ1990 = occ)
occda[which(occda$occ2_exec == 1 | occda$occ2_mgmtrel == 1),              "da"] <- "managers/executives"
occda[which(occda$occ2_prof == 1),                                         "da"] <- "professionals"
occda[which(occda$occ2_tech == 1),                                         "da"] <- "technicians"
occda[which(occda$occ2_finsales == 1 | occda$occ2_retsales == 1),          "da"] <- "sales"
occda[which(occda$occ2_cleric == 1 | occda$occ2_firepol == 1),             "da"] <- "administrative/office"
occda[which(occda$occ2_product == 1),                                      "da"] <- "production"
occda[which(occda$occ2_operator == 1 | occda$occ1_transmechcraft == 1),    "da"] <- "laborers"
occda[which(occda$occ3_clean == 1 | occda$occ3_protect == 1 |
              occda$occ3_guard == 1 | occda$occ3_janitor == 1),            "da"] <- "clean and protect services"
occda[which(occda$occ3_beauty == 1 | occda$occ3_recreation == 1 |
              occda$occ3_child == 1 | occda$occ3_othpers == 1 |
              occda$occ3_shealth == 1 | occda$occ3_food == 1),             "da"] <- "personal services"
occda <- occda %>% dplyr::select(occ1990, da)
## numeric DA codes: 1=PM, 3=tech, 4=sales, 5=admin, 6=prod, 7=labor, 8=clean/protect, 9=personal svc
occda[which(occda$da == "managers/executives"),   "da"] <- "1"
occda[which(occda$da == "professionals"),          "da"] <- "1"
occda[which(occda$da == "technicians"),            "da"] <- "3"
occda[which(occda$da == "sales"),                  "da"] <- "4"
occda[which(occda$da == "administrative/office"),  "da"] <- "5"
occda[which(occda$da == "production"),             "da"] <- "6"
occda[which(occda$da == "laborers"),               "da"] <- "7"
occda[which(occda$da == "clean and protect services"), "da"] <- "8"
occda[which(occda$da == "personal services"),      "da"] <- "9"
occda <- rbind(occda, data.frame(occ1990 = 349, da = "5"))
occda$da <- as.integer(occda$da)

## ── CPS ─────────────────────────────────────────────────────────────────────
load(file.path(data_dir, "cps_cleaned_2.RData"))
cps <- cps %>% filter(!is.na(occ1990) & !is.na(metarea) & !is.na(wtfinl))
cps$occ1990  <- as.character(cps$occ1990)
cps$metarea  <- as.character(cps$metarea)
cps[which(cps$egp == "VIIa"), "egp"] <- "VII"
cps[which(cps$egp == "IIIa"), "egp"] <- "III"
cps[which(cps$egp == "IIIb"), "egp"] <- "III"
cps[which(cps$da == "managers/executives"),    "da"] <- "1"
cps[which(cps$da == "professionals"),           "da"] <- "1"
cps[which(cps$da == "technicians"),             "da"] <- "3"
cps[which(cps$da == "sales"),                   "da"] <- "4"
cps[which(cps$da == "administrative/office"),   "da"] <- "5"
cps[which(cps$da == "production"),              "da"] <- "6"
cps[which(cps$da == "laborers"),                "da"] <- "7"
cps[which(cps$da == "clean and protect services"), "da"] <- "8"
cps[which(cps$da == "personal services"),       "da"] <- "9"
cps$da <- as.integer(cps$da)
cps[which(cps$earnweek > 9999),  "earnweek"]  <- NA
cps[which(cps$bpl == 99999),     "bpl"]        <- NA
cps[which(cps$uhrsworkt == 999), "uhrsworkt"]  <- NA
cps <- cps %>%
  dplyr::select(-c(age, labforce, classwkr, gw_macro, gw_meso, gw_micro,
                   race, hispan, empstat))

## ── ACS ─────────────────────────────────────────────────────────────────────
load(file.path(data_dir, "acs_cleaned_2.RData"))
acs$wtfinl  <- acs$perwt * acs$afactor
acs <- acs %>% filter(!is.na(occ1990) & !is.na(czone) & !is.na(wtfinl))
acs$occ1990 <- as.character(acs$occ1990)
acs$czone   <- as.character(acs$czone)
acs <- acs %>% dplyr::select(-puma) %>% rename(metarea = czone)
acs[which(acs$egp == "VIIa"), "egp"] <- "VII"
acs[which(acs$egp == "IIIa"), "egp"] <- "III"
acs[which(acs$egp == "IIIb"), "egp"] <- "III"
acs[which(acs$da == "managers/executives"),    "da"] <- "1"
acs[which(acs$da == "professionals"),           "da"] <- "1"
acs[which(acs$da == "technicians"),             "da"] <- "3"
acs[which(acs$da == "sales"),                   "da"] <- "4"
acs[which(acs$da == "administrative/office"),   "da"] <- "5"
acs[which(acs$da == "production"),              "da"] <- "6"
acs[which(acs$da == "laborers"),                "da"] <- "7"
acs[which(acs$da == "clean and protect services"), "da"] <- "8"
acs[which(acs$da == "personal services"),       "da"] <- "9"
acs$da <- as.integer(acs$da)
acs[which(acs$occ1990 == 349), "da"] <- 5
acs <- acs %>% dplyr::select(-c(gw_macro, gw_meso, gw_micro, afactor))

## CPI deflation (base = 2020 dollars, index value 399.2)
## Source: BLS CPI-U-RS; file should be in misc_data/CPI-U-RS_1978_to_2020.dta
cpi_path <- file.path(misc_dir, "CPI-U-RS_1978_to_2020.dta")
if (file.exists(cpi_path)) {
  cpi <- read_dta(cpi_path) %>% rename(year = income_year)
  acs <- merge(acs, cpi, by = "year", all.x = TRUE)
  acs$incwage_cpi <- acs$incwage / acs$cpi_u_rs_index * 399.2
} else {
  warning("CPI file not found at ", cpi_path, ". incwage_cpi will equal nominal incwage.")
  acs$cpi_u_rs_index <- 1
  acs$incwage_cpi    <- acs$incwage
}

## Occupation-level covariates (baseline ≤1990, from ACS/Census)
covariate <- acs %>%
  filter(incwage > 0 & uhrswork < 99 & uhrswork > 0, year <= 1990) %>%
  mutate(
    immigrant  = if_else(bpl > 120, 1L, 0L, missing = 0L),
    college    = if_else(educ >= 10 & educ <= 11, 1L, 0L),
    overwork   = if_else(uhrswork > 40, 1L, 0L),
    parttime   = if_else(uhrswork < 35, 1L, 0L),
    incwage    = incwage_cpi
  ) %>%
  group_by(occ1990) %>%
  summarize(
    incwage          = weighted.mean(incwage,                   wtfinl, na.rm = TRUE),
    uhrswork_male    = weighted.mean(uhrswork[sex == 1],        wtfinl[sex == 1], na.rm = TRUE),
    uhrswork_female  = weighted.mean(uhrswork[sex == 2],        wtfinl[sex == 2], na.rm = TRUE),
    uhrswork         = weighted.mean(uhrswork,                  wtfinl, na.rm = TRUE),
    share_female     = weighted.mean(sex,                       wtfinl, na.rm = TRUE) - 1,
    share_immigration= weighted.mean(immigrant,                 wtfinl, na.rm = TRUE),
    overwork_share   = weighted.mean(overwork,                  wtfinl, na.rm = TRUE),
    overwork_male    = weighted.mean(overwork[sex == 1],        wtfinl[sex == 1], na.rm = TRUE),
    overwork_female  = weighted.mean(overwork[sex == 2],        wtfinl[sex == 2], na.rm = TRUE),
    parttime_share   = weighted.mean(parttime,                  wtfinl, na.rm = TRUE),
    parttime_male    = weighted.mean(parttime[sex == 1],        wtfinl[sex == 1], na.rm = TRUE),
    parttime_female  = weighted.mean(parttime[sex == 2],        wtfinl[sex == 2], na.rm = TRUE),
    share_college    = weighted.mean(college,                   wtfinl, na.rm = TRUE),
    share_white      = weighted.mean(race == 0,                 wtfinl, na.rm = TRUE),
    share_black      = weighted.mean(race == 1,                 wtfinl, na.rm = TRUE),
    share_asian      = weighted.mean(race == 2,                 wtfinl, na.rm = TRUE),
    share_hispanic   = weighted.mean(race == 3,                 wtfinl, na.rm = TRUE),
    share_otherrace  = weighted.mean(race == 4,                 wtfinl, na.rm = TRUE),
    weight           = sum(wtfinl, na.rm = TRUE),
    .groups          = "drop"
  )

## CZ-level density and demographic covariates (baseline ≤1990)
lma <- read.table(file.path(misc_dir, "1990CZ_LMA.txt"),
                  as.is = TRUE, sep = "", head = TRUE,
                  strip.white = TRUE, fill = TRUE) %>%
  filter(FIPS != "Market" & !is.na(Labor_Force)) %>%
  rename(CZ = LMA.CZ) %>%
  mutate(LMA = floor(CZ / 100))

cz1990_shp <- read_sf(file.path(misc_dir, "cz1990_shapefile/cz1990.shp"))
cz1990_shp$size <- as.numeric(st_area(cz1990_shp)) / 1e6  # km²
cz1990_shp <- st_drop_geometry(cz1990_shp[, c("size", "cz")])

cz_density <- lma %>%
  mutate(Labor_Force = as.numeric(Labor_Force)) %>%
  group_by(CZ) %>%
  summarize(total_pop = sum(Labor_Force), .groups = "drop") %>%
  rename(metarea = CZ)

cz_immigrant <- acs %>%
  filter(incwage > 0 & uhrswork < 99 & uhrswork > 0, year <= 1990) %>%
  mutate(
    incwage_cpi = incwage_cpi / uhrswork,
    immigrant   = if_else(bpl > 120, 1L, 0L, missing = 0L),
    college     = if_else(educ >= 10 & educ <= 11, 1L, 0L),
    overwork    = if_else(uhrswork > 40, 1L, 0L)
  ) %>%
  group_by(metarea) %>%
  summarize(
    immigrant_share = sum(wtfinl[immigrant == 1], na.rm = TRUE) / sum(wtfinl, na.rm = TRUE),
    incwage_cpi     = weighted.mean(incwage_cpi, wtfinl, na.rm = TRUE),
    college         = weighted.mean(college, wtfinl, na.rm = TRUE),
    overwork        = weighted.mean(overwork, wtfinl, na.rm = TRUE),
    .groups = "drop"
  )

cz_density <- merge(cz_density, cz_immigrant, by = "metarea", all.x = TRUE)
cz_density <- merge(cz_density, cz1990_shp, by.x = "metarea", by.y = "cz", all.x = TRUE)
cz_density$density <- cz_density$total_pop / cz_density$size

## ── OEWS ────────────────────────────────────────────────────────────────────
oews <- data.frame(area = character(), occ_code = character(), group = character(),
                   h_median = character(), tot_emp = character(), year = character())
setwd(file.path(data_dir, "raw_data/OEWS_97-21"))

## 1999-2000
for (file in c("oes99ma/msa_1999_dl_1.xls", "oes99ma/msa_1999_dl_2.xls",
               "oes00ma/msa_2000_dl_1.xls", "oes00ma/msa_2000_dl_2.xls")) {
  msa <- read_excel(file) %>% filter(!is.na(!!sym(colnames(.)[1]))) %>%
    `colnames<-`(.[1, ]) %>% .[-1, ] %>%
    dplyr::select(area, occ_code, group, h_median, tot_emp)
  msa$year <- substr(file, 4, 5)
  oews <- rbind(oews, msa[is.na(msa$group), ])
}

## 2001-2002
for (file in c("oes01ma/msa_2001_dl_1.xls", "oes01ma/msa_2001_dl_2.xls",
               "oes01ma/msa_2001_dl_3.xls", "oes02ma/msa_2002_dl_1.xls",
               "oes02ma/msa_2002_dl_2.xls")) {
  msa <- read_excel(file) %>%
    dplyr::select(area, occ_code, group, h_median, tot_emp) %>% filter(is.na(group))
  msa$year <- substr(file, 4, 5); oews <- rbind(oews, msa)
}

## 2003-2011
files_0311 <- c(
  "oesm03ma/msa_may2003_dl_1.xls", "oesm03ma/msa_may2003_dl_2.xls",
  "oesm04ma/MSA_may2004_dl_1.xls", "oesm04ma/MSA_may2004_dl_2.xls", "oesm04ma/MSA_may2004_dl_3.xls",
  "oesm05ma/MSA_may2005_dl_1.xls", "oesm05ma/MSA_may2005_dl_2.xls", "oesm05ma/MSA_may2005_dl_3.xls",
  "oesm05ma/aMSA_may2005_dl.xls",
  "oesm06ma/MSA_may2006_dl_1.xls", "oesm06ma/MSA_may2006_dl_2.xls", "oesm06ma/MSA_may2006_dl_3.xls",
  "oesm06ma/aMSA_may2006_dl.xls",
  "oesm07ma/MSA_May2007_dl_1.xls", "oesm07ma/MSA_May2007_dl_2.xls", "oesm07ma/MSA_May2007_dl_3.xls",
  "oesm07ma/aMSA_May2007_dl.xls",
  "oesm08ma/MSA__M2008_dl_1.xls",  "oesm08ma/MSA_M2008_dl_2.xls",  "oesm08ma/MSA_M2008_dl_3.xls",
  "oesm08ma/aMSA__M2008_dl.xls",
  "oesm09ma/MSA_dl_1.xls", "oesm09ma/MSA_dl_2.xls", "oesm09ma/MSA_dl_3.xls",
  "oesm09ma/aMSA_M2009_dl.xls",
  "oesm10ma/MSA_M2010_dl_1.xls", "oesm10ma/MSA_M2010_dl_2.xls", "oesm10ma/MSA_M2010_dl_3.xls",
  "oesm10ma/aMSA_M2010_dl.xls",
  "oesm11ma/MSA_M2011_dl_1_AK_IN.xls", "oesm11ma/MSA_M2011_dl_2_KS_NY.xls",
  "oesm11ma/MSA_M2011_dl_3_OH_WY.xls", "oesm11ma/aMSA_M2011_dl.xls"
)
for (file in files_0311) {
  msa <- read_excel(file) %>% filter(!("Division" %in% AREA_NAME)) %>%
    dplyr::select(AREA, OCC_CODE, GROUP, H_MEDIAN, TOT_EMP) %>%
    rename(area = AREA, occ_code = OCC_CODE, group = GROUP, h_median = H_MEDIAN, tot_emp = TOT_EMP) %>%
    filter(is.na(group))
  msa$year <- substr(file, 5, 6)
  if (unique(msa$year) %in% c("03", "04")) oews <- rbind(oews, msa)
  else oews <- rbind(oews, msa %>% filter(substr(area, 5, 5) == "0"))
}

## 2012-2018
files_1218 <- c(
  "oesm12ma/MSA_M2012_dl_1_AK_IN.xls", "oesm12ma/MSA_M2012_dl_2_KS_NY.xls",
  "oesm12ma/MSA_M2012_dl_3_OH_WY.xls", "oesm12ma/aMSA_M2012_dl.xls",
  "oesm13ma/MSA_M2013_dl_1_AK_IN.xls", "oesm13ma/MSA_M2013_dl_2_KS_NY.xls",
  "oesm13ma/MSA_M2013_dl_3_OH_WY.xls", "oesm13ma/aMSA_M2013_dl.xls",
  "oesm14ma/MSA_M2014_dl.xlsx", "oesm14ma/aMSA_M2014_dl.xlsx",
  "oesm15ma/MSA_M2015_dl.xlsx", "oesm15ma/aMSA_M2015_dl.xlsx",
  "oesm16ma/MSA_M2016_dl.xlsx", "oesm16ma/aMSA_M2016_dl.xlsx",
  "oesm17ma/MSA_M2017_dl.xlsx", "oesm17ma/aMSA_M2017_dl.xlsx",
  "oesm18ma/MSA_M2018_dl.xlsx"
)
for (file in files_1218) {
  msa <- read_excel(file) %>% filter(!("Division" %in% AREA_NAME)) %>%
    dplyr::select(AREA, OCC_CODE, OCC_GROUP, H_MEDIAN, TOT_EMP) %>%
    rename(area = AREA, occ_code = OCC_CODE, group = OCC_GROUP, h_median = H_MEDIAN, tot_emp = TOT_EMP) %>%
    filter(group == "detailed")
  msa$year <- substr(file, 5, 6)
  oews <- rbind(oews, msa %>% filter(substr(area, 5, 5) == "0"))
}

## 2019
msa <- read_excel("oesm19ma/MSA_M2019_dl.xlsx") %>%
  dplyr::select(area, occ_code, o_group, h_median, tot_emp) %>%
  rename(group = o_group) %>% filter(group == "detailed") %>%
  filter(substr(area, 5, 5) == "0")
msa$year <- "19"; oews <- rbind(oews, msa)

## 2020-2021
for (file in c("oesm20ma/MSA_M2020_dl.xlsx", "oesm21ma/MSA_M2021_dl.xlsx")) {
  msa <- read_excel(file) %>%
    dplyr::select(AREA, OCC_CODE, O_GROUP, H_MEDIAN, TOT_EMP) %>%
    rename(area = AREA, occ_code = OCC_CODE, h_median = H_MEDIAN, group = O_GROUP, tot_emp = TOT_EMP) %>%
    filter(group == "detailed", substr(area, 5, 5) == "0")
  msa$year <- substr(file, 5, 6); oews <- rbind(oews, msa)
}

oews <- oews %>% dplyr::select(-group) %>%
  filter(occ_code != "00-0000", !is.na(occ_code)) %>%
  rename(metarea = area)

## year standardization
yr_map <- c("99"="1999","00"="2000","01"="2001","02"="2002","03"="2003","04"="2004",
            "05"="2005","06"="2006","07"="2007","08"="2008","09"="2009","10"="2010",
            "11"="2011","12"="2012","13"="2013","14"="2014","15"="2015","16"="2016",
            "17"="2017","18"="2018","19"="2019","20"="2020","21"="2021")
oews$year <- as.numeric(dplyr::recode(oews$year, !!!yr_map))

## OES99 → SOC2000 crosswalk
oes99 <- read_excel(file.path(misc_dir, "socoes98.xls")) %>%
  dplyr::select(soccode, oes99code) %>%
  filter(!oes99code %in% c("na.", "na"), !is.na(oes99code)) %>%
  distinct(oes99code, .keep_all = TRUE)
oews <- rbind(
  merge(oews[oews$year <= 2003, ], oes99, by.x = "occ_code", by.y = "oes99code", all.x = TRUE) %>%
    dplyr::select(-occ_code) %>% rename(occ_code = soccode),
  oews[oews$year > 2003, ]
)

## SOC2000 → SOC2010 crosswalk
soc0010 <- read_excel(file.path(misc_dir, "soc_2000_to_2010_crosswalk.xls")) %>%
  dplyr::select(`2000 SOC code`, `2010 SOC code`) %>%
  setNames(c("soc2000", "soc2010")) %>%
  distinct(soc2000, .keep_all = TRUE)
oews <- rbind(
  merge(oews[oews$year <= 2009, ], soc0010, by.x = "occ_code", by.y = "soc2000", all.x = TRUE) %>%
    dplyr::select(-occ_code) %>% rename(occ_code = soc2010),
  oews[oews$year >= 2010, ]
)

oews <- oews %>%
  filter(tot_emp != "**") %>%
  dplyr::select(metarea, tot_emp, year, occ_code, h_median) %>%
  filter(!is.na(occ_code))

oews$occ_code  <- as.character(oews$occ_code)
oews$metarea   <- as.character(oews$metarea)
oews$tot_emp   <- as.numeric(oews$tot_emp)
oews[which(oews$h_median %in% c("#", "*")), "h_median"] <- NA
oews$h_median  <- as.numeric(oews$h_median)

rm(msa, oes99, soc0010, files_0311, files_1218, file, yr_map,
   cz_immigrant)

cat("01 data.R: all objects loaded.\n")
