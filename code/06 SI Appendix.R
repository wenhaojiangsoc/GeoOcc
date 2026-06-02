## 06 SI Appendix.R
## Cross-dataset SI robustness checks:
##   SI S1: Raw vs SVD cosine similarity comparison (all datasets)
##   SI S2: k=5 (ACS/CPS) and k=15 (OEWS) robustness checks
##   SI S3: Placebo test (random geographic reallocation)
##   SI S3: Exclude PM workers in personal care industry
##   SI S4: LOO location regression table (reference note)
##   SI S5: OEWS occupation OLS
##
## Sources: 13_robust_raw_svd.R, 12_robust_k5.R, 14_robust_placebo_compute.R,
##          14_robust_placebo_plot.R, 15_robust_excl_pmpc.R, 10_S5_CPS_OEWS_occ_ols.R (OEWS parts)

repo_root <- "~/Dropbox/GeoOcc/GeoOccGit"
code_dir  <- file.path(repo_root, "code")
source(file.path(code_dir, "01 data.R"))
source(file.path(code_dir, "02 functions.R"))
select <- dplyr::select; filter <- dplyr::filter

library(reshape2)
library(ggrepel)
library(gridExtra)
library(grid)
library(tidyr)
library(patchwork)
library(scales)
library(ragg)

fig_dir <- file.path(repo_root, "figures", "si", "Cosine Similarity")


## ── SI S1: Raw vs SVD cosine similarity (all datasets) ───────────────────────
## → Output: Analysis/figures/raw_vs_svd_cosine_CPS.png
##           Analysis/figures/raw_vs_svd_cosine_ACS.png
##           Analysis/figures/raw_vs_svd_cosine_OEWS.png
##           Analysis/figures/raw_vs_svd_cosine.png
##           (SI Appendix S1)

pairwise_cosine_s1 <- function(mat) {
  n <- nrow(mat)
  pairs <- expand.grid(i = 1:n, j = 1:n)
  pairs <- pairs[pairs$i < pairs$j, ]
  pairs$cos <- mapply(function(a, b) cosine(mat[a, ], mat[b, ])[1, 1],
                      pairs$i, pairs$j)
  return(pairs)
}

compare_raw_svd_s1 <- function(olm, k = 50, dataset_label = "") {
  olm[is.na(olm)] <- 0
  raw_pairs <- pairwise_cosine_s1(olm)
  sv <- svd(olm, nu = k, nv = k)
  occvec_svd <- olm %*% sv$v[, 1:k]
  rownames(occvec_svd) <- rownames(olm)
  svd_pairs <- pairwise_cosine_s1(occvec_svd)
  data.frame(raw = raw_pairs$cos, svd = svd_pairs$cos, dataset = dataset_label)
}

k_s1 <- 50
i_s1 <- 2010

## CPS
load(file.path(data_dir, "cps_cleaned_2.RData"))
cps_s1 <- cps %>% filter(!is.na(occ1990) & !is.na(metarea) & !is.na(wtfinl))
cps_s1$occ1990 <- as.character(cps_s1$occ1990)
cps_s1$metarea <- as.character(cps_s1$metarea)
olm_cps_s1 <- with(five_rd(cps_s1, i_s1),
                   questionr::wtd.table(occ1990, metarea, wtfinl)) %>%
  prop.table(margin = 2)
olm_cps_s1 <- scale(olm_cps_s1, center = TRUE, scale = TRUE)
cps_compare_s1 <- compare_raw_svd_s1(olm_cps_s1, k = k_s1, dataset_label = "CPS (2010)")
cat("CPS done:", nrow(cps_compare_s1), "pairs\n")
rm(cps_s1)

## ACS
load(file.path(data_dir, "acs_cleaned_2.RData"))
acs_s1 <- acs %>%
  mutate(wtfinl = perwt * afactor) %>%
  filter(!is.na(occ1990) & !is.na(czone) & !is.na(wtfinl))
acs_s1$occ1990 <- as.character(acs_s1$occ1990)
acs_s1$czone   <- as.character(acs_s1$czone)
acs_s1 <- acs_s1 %>% rename(metarea = czone)
olm_acs_s1 <- with(five_rd(acs_s1, i_s1),
                   questionr::wtd.table(occ1990, metarea, wtfinl)) %>%
  prop.table(margin = 2)
olm_acs_s1 <- scale(olm_acs_s1, center = TRUE, scale = TRUE)
acs_compare_s1 <- compare_raw_svd_s1(olm_acs_s1, k = k_s1, dataset_label = "ACS (2010)")
cat("ACS done:", nrow(acs_compare_s1), "pairs\n")
rm(acs_s1, acs)

## OEWS (already loaded from 01_data.R)
oews_clean_s1 <- oews %>% filter(!is.na(occ_code) & !is.na(metarea))
oews_clean_s1$occ_code <- as.character(oews_clean_s1$occ_code)
oews_clean_s1$metarea  <- as.character(oews_clean_s1$metarea)
oews_clean_s1$tot_emp  <- as.numeric(oews_clean_s1$tot_emp)
olm_oews_s1 <- with(five_rd(oews_clean_s1, i_s1),
                    questionr::wtd.table(occ_code, metarea, weights = tot_emp)) %>%
  prop.table(margin = 2)
olm_oews_s1 <- scale(olm_oews_s1, center = TRUE, scale = TRUE)
oews_compare_s1 <- compare_raw_svd_s1(olm_oews_s1, k = k_s1, dataset_label = "OEWS (2010)")
cat("OEWS done:", nrow(oews_compare_s1), "pairs\n")

all_compare_s1 <- rbind(cps_compare_s1, acs_compare_s1, oews_compare_s1)

## individual plots
plot_raw_svd_s1 <- function(df, filename) {
  p <- ggplot(df, aes(x = raw, y = svd)) +
    geom_point(alpha = 0.1, size = 0.3, color = "grey30") +
    geom_abline(slope = 1, intercept = 0, lty = 2, color = "red") +
    xlab("Raw cosine similarity") +
    ylab(paste0("SVD-projected cosine similarity (k = ", k_s1, ")")) +
    ggtitle(unique(df$dataset)) +
    theme_classic() +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, hjust = 0.5))
  ggsave(filename, plot = p, width = 12, height = 12, units = "cm", dpi = 800)
  cat("Saved:", filename, "\n")
}

raw_fig_dir <- file.path(repo_root, "figures", "si")
plot_raw_svd_s1(cps_compare_s1,  file.path(raw_fig_dir, "raw_vs_svd_cosine_CPS.png"))
plot_raw_svd_s1(acs_compare_s1,  file.path(raw_fig_dir, "raw_vs_svd_cosine_ACS.png"))
plot_raw_svd_s1(oews_compare_s1, file.path(raw_fig_dir, "raw_vs_svd_cosine_OEWS.png"))

## combined faceted plot — k=50
svd_basics_dir <- file.path(repo_root, "figures", "si", "SVD Basics")
ggplot(all_compare_s1, aes(x = raw, y = svd)) +
  geom_point(alpha = 0.1, size = 0.3, color = "grey30") +
  geom_abline(slope = 1, intercept = 0, lty = 2, color = "red") +
  facet_wrap(~ dataset, ncol = 3, scales = "free") +
  xlab("Raw cosine similarity") +
  ylab(paste0("SVD-projected cosine similarity (k = ", k_s1, ")")) +
  theme_classic() +
  theme(strip.text = element_text(size = 14), axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggsave(file.path(svd_basics_dir, "raw_vs_svd_cosine.png"),
       width = 30, height = 10, units = "cm", dpi = 800)
cat("SI S1: raw_vs_svd_cosine.png (k=50) saved.\n")

## combined faceted plot — k=5
k_s1 <- 5
compare_raw_svd_k5 <- function(olm, k = 5, dataset_label = "") {
  olm[is.na(olm)] <- 0
  raw_pairs <- pairwise_cosine_s1(olm)
  sv <- svd(olm, nu = k, nv = k)
  occvec_svd <- olm %*% sv$v[, 1:k]
  rownames(occvec_svd) <- rownames(olm)
  svd_pairs <- pairwise_cosine_s1(occvec_svd)
  data.frame(raw = raw_pairs$cos, svd = svd_pairs$cos, dataset = dataset_label)
}
cps_compare_k5  <- compare_raw_svd_k5(olm_cps_s1,  k = 5, dataset_label = "CPS (2010)")
acs_compare_k5  <- compare_raw_svd_k5(olm_acs_s1,  k = 5, dataset_label = "ACS (2010)")
oews_compare_k5 <- compare_raw_svd_k5(olm_oews_s1, k = 5, dataset_label = "OEWS (2010)")
all_compare_k5  <- rbind(cps_compare_k5, acs_compare_k5, oews_compare_k5)

ggplot(all_compare_k5, aes(x = raw, y = svd)) +
  geom_point(alpha = 0.1, size = 0.3, color = "grey30") +
  geom_abline(slope = 1, intercept = 0, lty = 2, color = "red") +
  facet_wrap(~ dataset, ncol = 3, scales = "free") +
  xlab("Raw cosine similarity") +
  ylab(paste0("SVD-projected cosine similarity (k = ", k_s1, ")")) +
  theme_classic() +
  theme(strip.text = element_text(size = 14), axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggsave(file.path(svd_basics_dir, "raw_vs_svd_cosine_k5.png"),
       width = 30, height = 10, units = "cm", dpi = 800)
cat("SI S1: raw_vs_svd_cosine_k5.png (k=5) saved.\n")
k_s1 <- 50  ## restore for any downstream use


## ── SI S2: k=5 (ACS/CPS) and k=15 (OEWS) robustness ─────────────────────────
## → Output: Analysis/figures/Cosine Similarity/ACS_DA_PM_base_k5.png
##           Analysis/figures/Cosine Similarity/CPS_DA_PM_base_k5.png
##           Analysis/figures/Cosine Similarity/OEWS_M_base_k15.png
##           (SI Appendix S2)
##
## Full code is in Analysis/code/12_robust_k5.R (800+ lines).
## Run that file directly for complete k=5/k=15 robustness figures.
cat("SI S2 (k=5/k=15 robustness): Run the extended k-robustness pipeline (not included; contact authors).\n")
cat("Key outputs: ACS_DA_PM_base_k5.png, CPS_DA_PM_base_k5.png, OEWS_M_base_k15.png\n")


## ── SI S3: Placebo test (random geographic reallocation) ─────────────────────
## → Output: Analysis/figures/Cosine Similarity/placebo_ACS_DA_k50.png
##           Analysis/figures/Cosine Similarity/placebo_combined.png
##           (SI Appendix S3)

## Phase 1-3: compute (slow if not pre-saved)
## Load pre-saved placebo results if available; otherwise compute.
placebo_rdata <- file.path(repo_root, "data", "intermediate", "placebo_results.RData")

if (file.exists(placebo_rdata)) {
  cat("Loading pre-saved placebo results...\n")
  load(placebo_rdata)
  ## placebo_results.RData provides: years_acs, cos_obs_acs, cos_plac_acs,
  ##                                 years_cps, cos_obs_cps, cos_plac_cps,
  ##                                 years_oews, cos_obs_oews, cos_plac_oews
} else {
  cat("placebo_results.RData not found. Running placebo compute...\n")
  cat("(This may take several minutes — see README.md for compute instructions.)\n")
  source(file.path(code_dir, "14_robust_placebo_compute.R"))
  cat("Compute done. Re-run this script to produce the combined plot.\n")
  stop("Placebo compute completed. Reload and re-run to produce combined figure.")
}

## Plot combined figure from saved results
placebo_panel <- function(years, cos_obs, cos_plac, x_breaks, x_labels, title,
                          show_legend = FALSE) {
  plac_lo   <- apply(cos_plac, 1, quantile, 0.025, na.rm = TRUE)
  plac_hi   <- apply(cos_plac, 1, quantile, 0.975, na.rm = TRUE)
  plac_mean <- rowMeans(cos_plac, na.rm = TRUE)
  all_vals  <- c(cos_obs, plac_lo, plac_hi)
  lo   <- floor(min(all_vals, na.rm = TRUE) / 0.2) * 0.2 - 0.2
  hi   <- ceiling(max(all_vals, na.rm = TRUE) / 0.2) * 0.2 + 0.2
  step <- if ((hi - lo) <= 2.5) 0.4 else 0.7

  plot_df <- data.frame(year = years, observed = cos_obs,
                        plac_mean = plac_mean, plac_lo = plac_lo, plac_hi = plac_hi)

  ggplot(plot_df, aes(x = year)) +
    geom_ribbon(aes(ymin = plac_lo, ymax = plac_hi, fill = "Permutation (95% CI)"), alpha = 0.4) +
    geom_line(aes(y = plac_mean, color = "Permutation (95% CI)"), lty = 2, linewidth = 0.8) +
    geom_point(aes(y = observed, color = "Observed"), shape = 21, fill = "#82243b", size = 3.5) +
    geom_line(aes(y = observed, color = "Observed"), linewidth = 1) +
    scale_color_manual(name = NULL,
                       values = c("Observed" = "#82243b", "Permutation (95% CI)" = "grey40"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(1,2), shape = c(21,NA), fill = c("#82243b","grey70"),
                         linewidth = c(1,0.8), alpha = c(1,0.6)))) +
    scale_fill_manual(name = NULL, values = c("Permutation (95% CI)" = "grey70"), guide = "none") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(lo, hi), breaks = seq(lo, hi, step)) +
    ylab("cosine similarity (z-standardized)") + ggtitle(title) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5),
      axis.title.x = element_text(size = 16, hjust = 1, vjust = 8, margin = margin(t = -12)),
      axis.title.y = element_text(size = 16, angle = 90),
      plot.title   = element_text(size = 16, hjust = 0.5),
      legend.position = if (show_legend) "bottom" else "none",
      legend.direction = "vertical", legend.text = element_text(size = 16),
      legend.key.width = unit(1.5, "cm")
    )
}

get_legend_grob_p <- function(p) {
  gt  <- ggplot_gtable(ggplot_build(p))
  idx <- which(sapply(gt$grobs, function(x) x$name) == "guide-box")
  gt$grobs[[idx]]
}

p_acs_plac <- placebo_panel(
  years_acs, cos_obs_acs, cos_plac_acs,
  x_breaks = c(2002, 2005, 2008, 2011, 2014, 2017, 2020),
  x_labels = c("1980", "2005", "2008", "2011", "2014", "2017", "2020"),
  title    = "ACS/Census, 1980–2021"
)
p_cps_plac <- placebo_panel(
  years_cps, cos_obs_cps, cos_plac_cps,
  x_breaks = seq(1994, 2022, 4),
  x_labels = as.character(seq(1994, 2022, 4)),
  title    = "CPS, 1994–2022"
)
p_oews_plac <- placebo_panel(
  years_oews, cos_obs_oews, cos_plac_oews,
  x_breaks = seq(1999, 2020, 3),
  x_labels = as.character(seq(1999, 2020, 3)),
  title    = "OEWS, 1999–2020"
)

legend_grob_p <- get_legend_grob_p(
  placebo_panel(years_acs, cos_obs_acs, cos_plac_acs,
                c(2002, 2005, 2008, 2011, 2014, 2017, 2020),
                c("1980","2005","2008","2011","2014","2017","2020"),
                title = "", show_legend = TRUE)
)

panels_p  <- arrangeGrob(p_acs_plac, p_cps_plac, p_oews_plac, ncol = 3)
combined_p <- arrangeGrob(panels_p, legend_grob_p, nrow = 2, heights = unit(c(10, 1.5), "cm"))

ggsave(file.path(fig_dir, "placebo_combined.png"),
       plot = combined_p, width = 33, height = 11.5, units = "cm", dpi = 800)
cat("SI S3: placebo_combined.png saved.\n")


## ── SI S3: Exclude PM workers in personal care industry ───────────────────────
## → Output: Analysis/figures/Cosine Similarity/excl_pmpc_combined.png
##           (SI Appendix S3)

nu_ex <- 50
vec_cols_ex <- paste0("V", 1:nu_ex)

make_centroid_ex <- function(df, vc, group_col, weight_col) {
  df %>%
    filter(!is.na(.data[[group_col]]), !is.na(.data[[weight_col]])) %>%
    group_by(.data[[group_col]]) %>%
    summarize(across(all_of(vc), ~weighted.mean(., w = .data[[weight_col]], na.rm = TRUE)),
              .groups = "drop") %>%
    arrange(.data[[group_col]]) %>%
    select(all_of(vc)) %>% as.matrix()
}

zstd_pair_ex <- function(centroid, row_hi, row_lo) {
  n <- nrow(centroid); coslist <- c()
  for (k in 1:n) for (m in 1:n)
    if (k != m) coslist <- c(coslist, cosine(centroid[k, ], centroid[m, ])[1, 1])
  mu <- mean(coslist, na.rm = TRUE); sg <- sd(coslist, na.rm = TRUE)
  (cosine(centroid[row_hi, ], centroid[row_lo, ])[1, 1] - mu) / sg
}

## -- 1. ACS proxy exclusion (occ1990 15 and 19) --------------------------------
cat("Loading ACS for excl-PMPC analysis...\n")
load(file.path(data_dir, "acs_cleaned_2.RData"))
acs_ex <- acs %>%
  mutate(wtfinl = perwt * afactor, occ1990 = as.integer(acs$occ1990),
         czone = as.character(czone)) %>%
  filter(!is.na(occ1990) & !is.na(czone) & !is.na(wtfinl) & !is.na(da)) %>%
  rename(metarea = czone)
for (yr in list(c(1980,2002), c(1990,2003), c(2000,2004)))
  acs_ex[acs_ex$year == yr[1], "year"] <- yr[2]

pc_occ_proxy <- c(15L, 19L)
acs_excl_ex <- acs_ex %>% filter(!(da == 1 & occ1990 %in% pc_occ_proxy))
acs_agg_excl_ex <- acs_excl_ex %>%
  group_by(occ1990, year, da, metarea) %>%
  summarize(wtfinl = sum(wtfinl, na.rm = TRUE), .groups = "drop")
Nocc_acs_ex <- acs_ex %>%
  group_by(year, occ1990) %>%
  summarize(Nocc = sum(wtfinl, na.rm = TRUE), .groups = "drop")
occda_acs_ex <- acs_agg_excl_ex %>% distinct(occ1990, da)
rm(acs_ex, acs_excl_ex, acs)

years_acs_ex <- 2002:2021
cos_sim <- data.frame(year = years_acs_ex)
for (i in years_acs_ex) {
  occvec <- create_occvec(data = acs_agg_excl_ex, nu = 50, crosswalk = occda, cps = FALSE)
  occvec <- occvec_centroid(data = occvec)
  cosine_similarity(occvec, category = 8, standardize = TRUE)
  cat("  ACS excl year", i, "\n")
}
cos_restr_acs_ex <- cos_sim[, 9]  ## col 9 = DA=1 vs DA=9 (f=1, s=8 with category=8)

## -- 2. CPS direct ind1990 exclusion -------------------------------------------
cat("Loading CPS for excl-PMPC analysis...\n")
load(file.path(data_dir, "cps_cleaned_2.RData"))
cps_ex <- cps
cps_ex$occ1990 <- as.integer(as.numeric(cps_ex$occ1990))
cps_ex$ind1990 <- as.integer(as.numeric(cps_ex$ind1990))
cps_ex$metarea <- as.character(cps_ex$metarea)
cps_ex <- cps_ex %>% filter(!is.na(occ1990) & !is.na(metarea) & !is.na(wtfinl) & !is.na(da))

pc_ind_cps <- c(761, 771, 772, 780, 781, 791, 832, 862, 863, 870)
cps_excl_ex <- cps_ex %>% filter(!(da == 1 & !is.na(ind1990) & ind1990 %in% pc_ind_cps))
cps_agg_excl_ex <- cps_excl_ex %>%
  group_by(occ1990, year, da, metarea) %>%
  summarize(wtfinl = sum(wtfinl, na.rm = TRUE), .groups = "drop")
Nocc_cps_ex <- cps_ex %>%
  group_by(year, occ1990) %>%
  summarize(Nocc = sum(wtfinl, na.rm = TRUE), .groups = "drop")
occda_cps_ex <- cps_agg_excl_ex %>% distinct(occ1990, da)
rm(cps_ex, cps_excl_ex)

years_cps_ex <- 1994:2022
cos_sim <- data.frame(year = years_cps_ex)
for (i in years_cps_ex) {
  occvec <- create_occvec(data = cps_agg_excl_ex, nu = 50, crosswalk = occda, cps = TRUE)
  occvec <- occvec_centroid(data = occvec)
  cosine_similarity(occvec, category = 8, standardize = TRUE)
  cat("  CPS excl year", i, "\n")
}
cos_restr_cps_raw_ex <- cos_sim[, 9]
cos_restr_cps_ex <- c(cos_restr_cps_raw_ex[1],
                      rollmean(cos_restr_cps_raw_ex, 3),
                      cos_restr_cps_raw_ex[length(cos_restr_cps_raw_ex)])

## -- 3. OEWS exclusion of PM-in-PC SOC codes ----------------------------------
cat("Loading OEWS for excl-PMPC analysis...\n")
pc_mgmt_soc <- c("11-9111", "11-9151", "11-9061", "11-9171")
oews_excl_ex <- oews %>% filter(!occ_code %in% pc_mgmt_soc)
occind_excl_ex <- data.frame(occ_code = unique(oews_excl_ex$occ_code), stringsAsFactors = FALSE)
occind_excl_ex$ind <- substr(occind_excl_ex$occ_code, 1, 2)
Nocc_oews_ex <- oews %>%
  group_by(year, occ_code) %>%
  summarize(Nocc = sum(tot_emp, na.rm = TRUE), .groups = "drop")

cos_restr_oews_ex <- numeric(length(years_oews))
for (yi in seq_along(years_oews)) {
  i <- years_oews[yi]
  dat_ex <- oews_excl_ex %>%
    filter(year >= i - 2 & year <= i + 2 & !is.na(occ_code) & !is.na(metarea))
  olm_ex <- questionr::wtd.table(dat_ex$occ_code, dat_ex$metarea, weights = dat_ex$tot_emp) %>%
    prop.table(margin = 2)
  olm_ex <- scale(olm_ex, center = TRUE, scale = TRUE); olm_ex[is.na(olm_ex)] <- 0
  sv_ex  <- svd(olm_ex, nu = nu_ex, nv = nu_ex)
  ovm_ex <- sv_ex$u %*% diag(sv_ex$d[1:nu_ex])
  ovm_ex <- t(apply(ovm_ex, 1, l2.norm))
  df_ex  <- as.data.frame(ovm_ex); colnames(df_ex) <- vec_cols_ex
  df_ex$occ_code <- rownames(olm_ex)
  df_ex <- df_ex %>%
    left_join(occind_excl_ex, by = "occ_code") %>%
    left_join(Nocc_oews_ex %>% filter(year %in% (i-2):(i+2)) %>%
                group_by(occ_code) %>% summarize(Nocc = sum(Nocc, na.rm=TRUE)), by = "occ_code")
  ct_ex <- make_centroid_ex(df_ex, vec_cols_ex, "ind", "Nocc")
  ind_levels_ex <- sort(unique(df_ex$ind[!is.na(df_ex$ind)]))
  mgmt_row_ex   <- which(ind_levels_ex == "11")
  pc_row_ex     <- which(ind_levels_ex == "39")
  if (length(mgmt_row_ex) == 1 & length(pc_row_ex) == 1)
    cos_restr_oews_ex[yi] <- zstd_pair_ex(ct_ex, mgmt_row_ex, pc_row_ex)
  cat("  OEWS excl year", i, ":", round(cos_restr_oews_ex[yi], 3), "\n")
}

## -- Combined 3-panel comparison figure ----------------------------------------
compare_panel_ex <- function(years, cos_orig, cos_restr, x_breaks, x_labels, title,
                              show_legend = FALSE) {
  all_vals <- c(cos_orig, cos_restr)
  lo   <- floor(min(all_vals, na.rm = TRUE) / 0.2) * 0.2 - 0.2
  hi   <- ceiling(max(all_vals, na.rm = TRUE) / 0.2) * 0.2 + 0.2
  step <- if ((hi - lo) <= 2.5) 0.4 else 0.7

  plot_df <- data.frame(year = rep(years, 2),
                        value = c(cos_orig, cos_restr),
                        series = rep(c("Original", "Excl. PM in personal care (Autor-Dorn)"),
                                     each = length(years)))
  plot_df$series <- factor(plot_df$series,
                           levels = c("Original", "Excl. PM in personal care (Autor-Dorn)"))

  ggplot(plot_df, aes(x = year, y = value, color = series, linetype = series, shape = series)) +
    geom_point(aes(fill = series), size = 3.5) + geom_line(linewidth = 0.9) +
    scale_color_manual(name = NULL,
                       values = c("Original"="#82243b",
                                  "Excl. PM in personal care (Autor-Dorn)"="grey40")) +
    scale_fill_manual(name = NULL,
                      values = c("Original"="#82243b",
                                 "Excl. PM in personal care (Autor-Dorn)"="grey40")) +
    scale_shape_manual(name = NULL, values = c("Original"=21,
                                               "Excl. PM in personal care (Autor-Dorn)"=21)) +
    scale_linetype_manual(name = NULL, values = c("Original"=1,
                                                  "Excl. PM in personal care (Autor-Dorn)"=2)) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(lo, hi), breaks = seq(lo, hi, step)) +
    ylab("cosine similarity (z-standardized)") + ggtitle(title) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5),
      axis.title.x = element_text(size = 16, hjust = 1, vjust = 8, margin = margin(t = -12)),
      axis.title.y = element_text(size = 16, angle = 90),
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.position = if (show_legend) "bottom" else "none",
      legend.direction = "vertical", legend.text = element_text(size = 16),
      legend.key.width = unit(1.5, "cm")
    )
}

get_legend_grob_ex <- function(p) {
  gt  <- ggplot_gtable(ggplot_build(p))
  idx <- which(sapply(gt$grobs, function(x) x$name) == "guide-box")
  gt$grobs[[idx]]
}

p_acs_ex  <- compare_panel_ex(years_acs_ex, cos_obs_acs,  cos_restr_acs_ex,
                               c(2002,2005,2008,2011,2014,2017,2020),
                               c("1980","2005","2008","2011","2014","2017","2020"),
                               "ACS/Census, 1980–2021")
p_cps_ex  <- compare_panel_ex(years_cps_ex, cos_obs_cps,  cos_restr_cps_ex,
                               seq(1994,2022,4), as.character(seq(1994,2022,4)),
                               "CPS, 1994–2022")
p_oews_ex <- compare_panel_ex(years_oews, cos_obs_oews, cos_restr_oews_ex,
                               seq(1999,2020,3), as.character(seq(1999,2020,3)),
                               "OEWS, 1999–2020")

legend_grob_ex <- get_legend_grob_ex(
  compare_panel_ex(years_acs_ex, cos_obs_acs, cos_restr_acs_ex,
                   c(2002), c("1980"), "", show_legend = TRUE)
)

panels_ex   <- arrangeGrob(p_acs_ex, p_cps_ex, p_oews_ex, ncol = 3)
combined_ex <- arrangeGrob(panels_ex, legend_grob_ex, nrow = 2, heights = unit(c(10, 1.5), "cm"))
ggsave(file.path(fig_dir, "excl_pmpc_combined.png"),
       plot = combined_ex, width = 33, height = 11.5, units = "cm", dpi = 800)
cat("SI S3: excl_pmpc_combined.png saved.\n")


## ── SI S4: LOO location regression table ─────────────────────────────────────
## SI S4 LOO regression table — requires ACS LOO results from 03 ACS Analysis.R
## The table in si_appendix.tex is hardcoded; rerun 03 ACS Analysis.R to regenerate.
## LOO scatter figures: LOO_location_OEWS.png (04 OEWS Analysis.R),
##                     LOO_location_CPS.png (05 CPS Analysis.R),
##                     ACS_DA_location_heterogeneity.png (03 ACS Analysis.R).


## ── SI S5: OEWS occupation OLS ───────────────────────────────────────────────
## → Output: Analysis/figures/Cosine Similarity/OEWS_occ_ols.png
##           (SI Appendix S5)

oews_loo_rds <- file.path(data_dir, "intermediate/oews_occ_loo_m1.rds")
if (!file.exists(oews_loo_rds)) {
  cat("SI S5 OEWS OLS: oews_occ_loo_m1.rds not found.\n")
  cat("data/intermediate/oews_occ_loo_m1.rds not found — see README.md for how to generate it.\n")
} else {
  library(haven)
  library(readxl)

  oews_m1 <- readRDS(oews_loo_rds) %>%
    mutate(occ_code = as.character(occ_code),
           da = if_else(substr(occ_code, 1, 2) == "11", 1L, 9L))

  onet_soc39 <- readRDS(file.path(data_dir, "intermediate/onet_soc39.rds")) %>%
    mutate(soc_code = as.character(soc_code))

  ## SOC crosswalk
  xw1 <- read_excel(file.path(misc_dir, "soc2010_to_occ2010.xlsx")) %>%
    dplyr::select(OCC2010, SOC2010) %>% mutate(occ_num = as.numeric(OCC2010))
  xw2 <- read_dta(file.path(misc_dir, "occ2010_occ1990dd.dta")) %>%
    mutate(occ_num = as.numeric(occ)) %>%
    dplyr::select(occ_num, occ1990 = occ1990dd) %>%
    filter(!is.na(occ1990), occ1990 < 900)
  xwalk <- xw1 %>% left_join(xw2, by = "occ_num") %>%
    dplyr::select(SOC2010, occ1990) %>% filter(!is.na(occ1990)) %>%
    mutate(occ1990 = as.character(occ1990)) %>% distinct(SOC2010, .keep_all = TRUE)

  ## OEWS wages 2002-2003
  oews_wages_early <- oews %>%
    filter(year %in% c(2002, 2003), !is.na(h_median), h_median > 0, !is.na(tot_emp), tot_emp > 0) %>%
    group_by(occ_code) %>%
    summarize(incwage = weighted.mean(h_median, tot_emp, na.rm=TRUE) * 2080,
              weight  = sum(tot_emp, na.rm=TRUE), .groups = "drop")

  ## attach covariates
  cov_df_s5 <- covariate %>% mutate(occ1990 = as.character(occ1990))

  oews_ols_df <- oews_m1 %>%
    left_join(oews_wages_early, by = "occ_code") %>%
    left_join(xwalk, by = c("occ_code" = "SOC2010")) %>%
    left_join(cov_df_s5 %>% dplyr::select(-incwage, -any_of("weight")), by = "occ1990") %>%
    left_join(onet_soc39 %>% rename(occ_code = soc_code), by = "occ_code")

  oews_ols_df <- oews_ols_df %>%
    group_by(da) %>%
    mutate(
      influence_z       = scale(influence)[,1],
      log_incwage_z     = scale(log(incwage))[,1],
      overwork_male_z   = scale(overwork_male)[,1],
      overwork_female_z = scale(overwork_female)[,1],
      fulltime_z        = scale(1 - parttime_share)[,1],
      uhrswork_male_z   = scale(uhrswork_male)[,1],
      uhrswork_female_z = scale(uhrswork_female)[,1],
      share_female_z    = scale(share_female)[,1],
      share_immigration_z = scale(share_immigration)[,1],
      physical_prox_z   = scale(physical_prox)[,1],
      face2face_z       = scale(face2face)[,1],
      interpersonal_z   = scale(interpersonal)[,1],
      ext_customers_z   = scale(ext_customers)[,1],
      responsible_z     = scale(responsible_health)[,1],
      share_white_z     = scale(share_white)[,1],
      share_black_z     = scale(share_black)[,1],
      share_asian_z     = scale(share_asian)[,1],
      share_hispanic_z  = scale(share_hispanic)[,1]
    ) %>% ungroup() %>%
    mutate(weight = if_else(da == 9L & (is.na(weight) | weight <= 0), 1, weight))

  pm_vars_s5 <- c("log_incwage_z","overwork_male_z","overwork_female_z",
                  "uhrswork_male_z","uhrswork_female_z","fulltime_z",
                  "share_white_z","share_black_z","share_asian_z","share_hispanic_z")
  pm_labels_s5 <- c("Annual income (log)","Share of overwork (male)","Share of overwork (female)",
                    "    Weekly working hours (male)","    Weekly working hours (female)",
                    "Share of full-time working","Share of white","Share of Black",
                    "Share of Asian","Share of Hispanic")
  svc_vars_s5 <- c("share_white_z","share_black_z","share_asian_z","share_hispanic_z",
                   "share_immigration_z","physical_prox_z","face2face_z","interpersonal_z",
                   "ext_customers_z","responsible_z")
  svc_labels_s5 <- c("Share of white","Share of Black","Share of Asian","Share of Hispanic",
                     "Share of immigrants","Physical proximity","Face-to-face discussions",
                     "Interpersonally oriented","Serves external clients","Care responsibility")
  pm_features_s5 <- c("demand","demand","demand","demand","demand","demand",
                      "race/ethnicity","race/ethnicity","race/ethnicity","race/ethnicity")
  svc_features_s5 <- c("race/ethnicity","race/ethnicity","race/ethnicity","race/ethnicity","",
                       "work content/style","work content/style","work content/style",
                       "work content/style","work content/style")

  results_oews_s5 <- data.frame(
    group = c(rep("PM", length(pm_vars_s5)), rep("Service", length(svc_vars_s5))),
    variable = c(pm_labels_s5, svc_labels_s5),
    feature = c(pm_features_s5, svc_features_s5),
    coef = NA_real_, se = NA_real_
  )

  d_pm_s5  <- oews_ols_df %>% filter(da == 1L)
  d_svc_s5 <- oews_ols_df %>% filter(da == 9L)

  for (ii in seq_along(pm_vars_s5)) {
    fit <- tryCatch(lm(as.formula(paste("influence_z ~", pm_vars_s5[ii])),
                       data = d_pm_s5, weights = weight), error = function(e) NULL)
    if (!is.null(fit) && nrow(summary(fit)$coefficients) >= 2) {
      results_oews_s5$coef[ii] <- summary(fit)$coefficients[2, "Estimate"]
      results_oews_s5$se[ii]   <- summary(fit)$coefficients[2, "Std. Error"]
    }
  }
  n_pm_s5 <- length(pm_vars_s5)
  for (ii in seq_along(svc_vars_s5)) {
    fit <- tryCatch(lm(as.formula(paste("influence_z ~", svc_vars_s5[ii])),
                       data = d_svc_s5, weights = weight), error = function(e) NULL)
    if (!is.null(fit) && nrow(summary(fit)$coefficients) >= 2) {
      results_oews_s5$coef[n_pm_s5 + ii] <- summary(fit)$coefficients[2, "Estimate"]
      results_oews_s5$se[n_pm_s5 + ii]   <- summary(fit)$coefficients[2, "Std. Error"]
    }
  }

  ## plot (same lollipop style)
  reorder_wf_s5 <- function(df) {
    df %>%
      group_by(feature) %>% arrange(desc(coef)) %>% mutate(.ord = row_number()) %>% ungroup() %>%
      group_by(feature) %>% mutate(fm = mean(coef, na.rm=TRUE)) %>% ungroup() %>%
      arrange(desc(fm), feature, .ord) %>%
      mutate(variable = factor(variable, levels = unique(variable)))
  }
  sep_ys5 <- function(d) {
    d %>% mutate(.y = as.numeric(variable)) %>%
      group_by(feature) %>% summarize(max_y = max(.y), .groups = "drop") %>%
      arrange(max_y) %>% pull(max_y) %>%
      { if (length(.) > 1) head(., -1) + 0.5 else numeric(0) }
  }

  df_ow_s5 <- results_oews_s5 %>% mutate(lo = coef - 1.96 * se, hi = coef + 1.96 * se)
  pm_ow  <- df_ow_s5 %>% filter(group == "PM") %>% mutate(variable = as.character(variable)) %>% reorder_wf_s5()
  svc_ow <- df_ow_s5 %>% filter(group == "Service") %>% mutate(variable = as.character(variable)) %>% reorder_wf_s5()
  x_m_pm_ow  <- max(abs(c(pm_ow$lo,  pm_ow$hi)),  na.rm = TRUE)
  x_m_svc_ow <- max(abs(c(svc_ow$lo, svc_ow$hi)), na.rm = TRUE)

  p_pm_ow <- ggplot(pm_ow, aes(x = coef, y = variable)) +
    { sep_ys5(pm_ow) %>% { if(length(.)>0) geom_hline(yintercept=., color="grey85", linewidth=0.35) } } +
    geom_vline(xintercept=0, linetype="dashed", color="grey60") +
    geom_errorbarh(aes(xmin=lo, xmax=hi), height=0, color="grey40") +
    geom_point(size=3.7, shape=21, fill="#82243b", color="grey40") +
    scale_y_discrete(limits=levels(pm_ow$variable)) +
    coord_cartesian(xlim=c(-x_m_pm_ow-0.2, x_m_pm_ow+0.2), expand=TRUE) +
    theme_classic() +
    labs(x=NULL, y=NULL, title="Cosine Change (OEWS, 1999-2017) and Occupational Features")

  p_svc_ow <- ggplot(svc_ow, aes(x = coef, y = variable)) +
    { sep_ys5(svc_ow) %>% { if(length(.)>0) geom_hline(yintercept=., color="grey85", linewidth=0.35) } } +
    geom_vline(xintercept=0, linetype="dashed", color="grey60") +
    geom_errorbarh(aes(xmin=lo, xmax=hi), height=0, color="grey40") +
    geom_point(size=3.7, shape=21, fill="#ffe4ec", color="grey40") +
    scale_y_discrete(limits=levels(svc_ow$variable)) +
    scale_x_continuous(labels=label_number(accuracy=0.1)) +
    coord_cartesian(xlim=c(-x_m_svc_ow-1, x_m_svc_ow+1), expand=TRUE) +
    theme_classic() +
    labs(x="coefficient (bivariate)", y=NULL, title=NULL) +
    theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=13))

  p_oews_s5 <- p_pm_ow / p_svc_ow +
    plot_layout(heights = c(nrow(pm_ow), nrow(svc_ow))) +
    theme(plot.margin = margin(t=0, r=0, b=0, l=0.5, unit="cm"))

  ggsave(file.path(fig_dir, "OEWS_occ_ols.png"),
         plot = p_oews_s5, width = 14, height = 17.5, units = "cm", dpi = 800)
  cat("Saved OEWS_occ_ols.png\n")
}

## ── Cleanup ───────────────────────────────────────────────────────────────────
rm(list = c(
  "all_compare_s1", "cps_compare_s1", "acs_compare_s1", "oews_compare_s1",
  "cps_compare_k5", "acs_compare_k5", "oews_compare_k5", "all_compare_k5",
  "olm_cps_s1", "olm_acs_s1", "olm_oews_s1",
  "oews_clean_s1", "cps_s1",
  "k_s1", "i_s1", "svd_basics_dir", "raw_fig_dir",
  "acs_agg_excl_ex", "cps_agg_excl_ex", "oews_excl_ex",
  "Nocc_acs_ex", "Nocc_cps_ex", "Nocc_oews_ex",
  "occda_acs_ex", "occda_cps_ex", "occind_excl_ex",
  "nu_ex", "vec_cols_ex",
  "cos_restr_acs_ex", "cos_restr_cps_raw_ex", "cos_restr_cps_ex", "cos_restr_oews_ex",
  "p_acs_plac", "p_cps_plac", "p_oews_plac", "panels_p", "combined_p",
  "p_acs_ex", "p_cps_ex", "p_oews_ex", "panels_ex", "combined_ex",
  "legend_grob_p", "legend_grob_ex",
  "years_acs_ex", "years_cps_ex"
)[sapply(c(
  "all_compare_s1", "cps_compare_s1", "acs_compare_s1", "oews_compare_s1",
  "cps_compare_k5", "acs_compare_k5", "oews_compare_k5", "all_compare_k5",
  "olm_cps_s1", "olm_acs_s1", "olm_oews_s1",
  "oews_clean_s1", "cps_s1",
  "k_s1", "i_s1", "svd_basics_dir", "raw_fig_dir",
  "acs_agg_excl_ex", "cps_agg_excl_ex", "oews_excl_ex",
  "Nocc_acs_ex", "Nocc_cps_ex", "Nocc_oews_ex",
  "occda_acs_ex", "occda_cps_ex", "occind_excl_ex",
  "nu_ex", "vec_cols_ex",
  "cos_restr_acs_ex", "cos_restr_cps_raw_ex", "cos_restr_cps_ex", "cos_restr_oews_ex",
  "p_acs_plac", "p_cps_plac", "p_oews_plac", "panels_p", "combined_p",
  "p_acs_ex", "p_cps_ex", "p_oews_ex", "panels_ex", "combined_ex",
  "legend_grob_p", "legend_grob_ex",
  "years_acs_ex", "years_cps_ex"
), exists)])

## ── Interactive ──────────────────────────────────────────────────────────────
cat("\nSI Appendix outputs summary:\n")
cat("  SI S1: Analysis/figures/SVD Basics/raw_vs_svd_cosine.png (k=50)\n")
cat("         Analysis/figures/SVD Basics/raw_vs_svd_cosine_k5.png (k=5)\n")
cat("  SI S2: *_k5.png, *_k15.png — run 12_robust_k5.R directly\n")
cat("  SI S3: placebo_combined.png, excl_pmpc_combined.png\n")
cat("  SI S4: LOO_location_OEWS.png, LOO_location_CPS.png (from dataset scripts)\n")
cat("  SI S5: OEWS_occ_ols.png, CPS_occ_ols.png (from 05 CPS Analysis.R)\n")
