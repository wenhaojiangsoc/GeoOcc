## 00 prelude.R
## Raw data cleaning — RUN ONCE to produce cps_cleaned_2.RData and acs_cleaned_2.RData
##
## This script reads raw IPUMS extracts and crosswalk files, cleans and recodes variables,
## and saves cleaned datasets to GeoOcc/Analysis/data/. It only needs to be run once;
## all analysis scripts load from the saved .RData files.
##
## Prerequisites:
##   - IPUMS CPS extract:  GeoOcc/Analysis/raw_data/cps_00006.dta
##   - IPUMS ACS extract:  GeoOcc/Analysis/raw_data/usa_00063.xml + .dat
##   - PUMA-CZ crosswalks: GeoOcc/Analysis/misc_data/PUMA_CZ/
##   - EGP/Grusky/Autor-Dorn crosswalks in GeoOcc/Analysis/misc_data/
##
## Outputs:
##   GeoOcc/Analysis/data/cps_cleaned_2.RData
##   GeoOcc/Analysis/data/acs_cleaned_2.RData

## ── Path configuration ──────────────────────────────────────────────────────
data_dir <- "~/Dropbox/GeoOcc/Analysis"
misc_dir  <- file.path(data_dir, "misc_data")
setwd(data_dir)

## ── Packages ────────────────────────────────────────────────────────────────
library(dplyr)
library(haven)
library(sjlabelled)
library(ipumsr)

## ── Crosswalks ──────────────────────────────────────────────────────────────
occegp <- read.csv("misc_data/occegp.csv") %>% arrange(occ1990)

occgw <- read_dta("misc_data/macro_meso_micro_occ/Grusky_occ_class_scheme.dta") %>%
  dplyr::select(-occ1990_labels) %>%
  rename(gw_macro = macro_adj, gw_meso = meso_adj, gw_micro = micro_adj)

occda <- read_dta("misc_data/Autor_Dorn_occ_class_scheme.dta") %>%
  rename(occ1990 = occ)
occda[which(occda$occ2_exec == 1 | occda$occ2_mgmtrel == 1), "da"] <- "managers/executives"
occda[which(occda$occ2_prof == 1),                             "da"] <- "professionals"
occda[which(occda$occ2_tech == 1),                             "da"] <- "technicians"
occda[which(occda$occ2_finsales == 1 | occda$occ2_retsales == 1), "da"] <- "sales"
occda[which(occda$occ2_cleric == 1 | occda$occ2_firepol == 1),    "da"] <- "administrative/office"
occda[which(occda$occ2_product == 1),                              "da"] <- "production"
occda[which(occda$occ2_operator == 1 | occda$occ1_transmechcraft == 1), "da"] <- "laborers"
occda[which(occda$occ3_clean == 1 | occda$occ3_protect == 1 |
              occda$occ3_guard == 1 | occda$occ3_janitor == 1), "da"] <- "clean and protect services"
occda[which(occda$occ3_beauty == 1 | occda$occ3_recreation == 1 | occda$occ3_child == 1 |
              occda$occ3_othpers == 1 | occda$occ3_shealth == 1 | occda$occ3_food == 1), "da"] <- "personal services"
occda <- occda %>% dplyr::select(occ1990, da)

## ── CPS cleaning ────────────────────────────────────────────────────────────
cps <- read_dta("raw_data/cps_00006.dta")
cps <- cps %>%
  dplyr::select(year, month, region, statefip, metarea, county, wtfinl, age, sex, race, hispan,
                empstat, bpl, educ, labforce, occ1990, ind1990, earnweek, uhrsworkt, classwkr)
cps <- cps %>% filter(age >= 25 & age <= 64 & labforce == 2)
cps <- cps %>% filter(occ1990 != 905 & occ1990 != 999 & occ1990 < 900)
cps <- merge(cps, occegp, by = "occ1990", all.x = TRUE)

## recode EGP based on self-employment
cps[which(cps$egp %in% c("IIIa","IIIb","V","VI","VIIa") & cps$classwkr == 14), "egp"] <- "IVa"
cps[which(cps$egp %in% c("IIIa","IIIb","V","VI","VIIa") & cps$classwkr == 13), "egp"] <- "IVb"
cps[which(cps$egp %in% c("VIIb")                         & cps$classwkr == 13), "egp"] <- "IVc"

cps <- merge(cps, occgw,  by = "occ1990", all.x = TRUE)
cps <- merge(cps, occda,  by = "occ1990", all.x = TRUE)
cps <- cps %>%
  mutate(metarea = floor(metarea / 10)) %>%
  filter(metarea < 999 & !is.na(metarea))

save(cps, file = "data/cps_cleaned_2.RData")
cat("Saved: cps_cleaned_2.RData\n")

## ── ACS cleaning ────────────────────────────────────────────────────────────
ddi <- read_ipums_ddi("raw_data/usa_00063.xml")
acs <- read_ipums_micro(ddi)
acs <- acs %>%
  rename_with(tolower) %>%
  dplyr::select(year, statefip, cntygp98, incwage, uhrswork, bpl, educ,
                puma, perwt, sex, age, labforce, classwkr, occ1990, occ2010, race, hispan)

acs <- acs %>% rename(race_2 = race)
acs[which(acs$race_2 == 1 & acs$hispan == 0),                                  "race"] <- 0  # white
acs[which(acs$race_2 == 2 & acs$hispan == 0),                                  "race"] <- 1  # black
acs[which(acs$race_2 %in% c(4,5,6) & acs$hispan == 0),                         "race"] <- 2  # Asian
acs[which(acs$hispan %in% c(1,2,3,4)),                                          "race"] <- 3  # Hispanic
acs[which(is.na(acs$race)),                                                      "race"] <- 4
acs <- acs %>% dplyr::select(-c(race_2, hispan))

acs <- acs %>% filter(age >= 25 & age <= 64 & labforce == 2 & year <= 2021)
acs <- acs %>% filter(occ1990 != 905 & occ1990 != 999 & occ1990 < 900)
acs <- merge(acs, occegp, by = "occ1990", all.x = TRUE)
acs[which(acs$egp %in% c("IIIa","IIIb","V","VI","VIIa") & acs$classwkr == 14), "egp"] <- "IVa"
acs[which(acs$egp %in% c("IIIa","IIIb","V","VI","VIIa") & acs$classwkr == 13), "egp"] <- "IVb"
acs[which(acs$egp %in% c("VIIb")                         & acs$classwkr == 13), "egp"] <- "IVc"
acs <- merge(acs, occgw,  by = "occ1990", all.x = TRUE)
acs <- merge(acs, occda,  by = "occ1990", all.x = TRUE)

## PUMA→CZ crosswalks
cz1980 <- read_dta("misc_data/PUMA_CZ/cw_ctygrp1980_czone/cw_ctygrp1980_czone_corr.dta")
acs[which(acs$year == 1980), "ctygrp1980"] <-
  acs[which(acs$year == 1980), "cntygp98"] + 1000 * acs[which(acs$year == 1980), "statefip"]
acs1980 <- merge(acs[which(acs$year == 1980), ], cz1980, by = "ctygrp1980") %>%
  dplyr::select(-c(ctygrp1980, cntygp98, statefip)) %>% remove_all_labels()

cz1990 <- read_dta("misc_data/PUMA_CZ/cw_puma1990_czone/cw_puma1990_czone.dta")
acs1990 <- acs[which(acs$year == 1990), ]
acs1990$puma1990 <- paste0(as.character(acs1990$statefip), sprintf("%04d", acs1990$puma))
acs1990 <- merge(acs1990, cz1990, by = "puma1990") %>%
  dplyr::select(-c(puma1990, cntygp98, statefip, ctygrp1980)) %>% remove_all_labels()

cz2000 <- read_dta("misc_data/PUMA_CZ/cw_puma2000_czone/cw_puma2000_czone.dta")
acs2000_11 <- acs[which(acs$year %in% c(2000, 2005:2011)), ]
acs2000_11$puma2000 <- paste0(as.character(acs2000_11$statefip), sprintf("%04d", acs2000_11$puma))
acs2000_11 <- merge(acs2000_11, cz2000, by = "puma2000") %>%
  dplyr::select(-c(puma2000, cntygp98, statefip, ctygrp1980)) %>% remove_all_labels()

cz2010 <- read_dta("misc_data/PUMA_CZ/cw_puma2010_czone/cw_puma2010_czone.dta")
acs2012_21 <- acs[which(acs$year >= 2012 & acs$year <= 2021), ]
acs2012_21$puma2010 <- paste0(as.character(acs2012_21$statefip), sprintf("%05d", acs2012_21$puma))
acs2012_21 <- merge(acs2012_21, cz2010, by = "puma2010") %>%
  dplyr::select(-c(puma2010, cntygp98, statefip, ctygrp1980)) %>% remove_all_labels()

acs <- rbind(acs1980, acs1990, acs2000_11, acs2012_21)
rm(acs1980, acs1990, acs2000_11, acs2012_21, cz1980, cz1990, cz2000, cz2010)
save(acs, file = "data/acs_cleaned_2.RData")
cat("Saved: acs_cleaned_2.RData\n")
