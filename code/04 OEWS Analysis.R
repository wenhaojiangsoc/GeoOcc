## 04 OEWS Analysis.R
## OEWS main cosine similarity figure (Fig 3 right panel), location LOO (SI S4),
## Table 1 OEWS column (Panel A lead-lag + Panel B Bartik IV),
## and LaTeX causal_table.tex (Table 1 assembly).
##
## Sources: 09_OEWS_main.R, 12h_oews_panel_comparison.R

repo_root <- "~/Dropbox/GeoOcc/GeoOccGit"
code_dir  <- file.path(repo_root, "code")
source(file.path(code_dir, "01 data.R"))
source(file.path(code_dir, "02 functions.R"))
select <- dplyr::select; filter <- dplyr::filter

library(reshape2)
library(ggrepel)
library(sf)
library(readxl)
library(gridExtra)
library(grid)
library(fixest)
library(tidyr)
library(ragg)

fig_dir <- file.path(repo_root, "figures", "si", "Cosine Similarity")

## Harmonize OEWS area codes to PMSA 1999 definitions across all years.
## This spatial crosswalk maps CBSA 2004-2014 and CBSA 2015-2021 codes to
## PMSA 1999 codes via nearest-neighbor join, enabling a balanced 1999-2021 panel.
library(nngeo)

## keep non-NA occ-codes and create industry crosswalk
## Merge SOC 13 (business/financial) into SOC 11 (management) for consistency with lead-lag and IV analyses
oews <- oews[which(!is.na(oews$occ_code)), ]
occind <- data.frame(occ_code = unique(oews$occ_code))
occind$ind <- substr(occind$occ_code, 1, 2)
occind$ind[occind$ind == "13"] <- "11"


## ── Section 1: OEWS main cosine similarity figure (Fig 3 right panel) ─────────
## → Output: Analysis/figures/Cosine Similarity/OEWS_M_base_3.png
##           (Main paper Fig 3, right panel)

cos_sim <- data.frame(year = 1999:2021)
for (i in 1999:2021) {
  occvec <- create_occvec_owes(data = oews, nu = 50, crosswalk = occind)
  occvec <- occvec_centroid(data = occvec)
  cosine_similarity(occvec, category = 21, standardize = TRUE)
  cat(sprintf("year %d done\n", i))
}

## SOC 2-digit groups sorted (21 total after merging 13→11):
## position 1 = "11+13" (management+business/financial), position 14 = "39" (personal care)
## Column for baseline vs group: 1 + position
groups_sorted <- sort(unique(occind$ind))
ind_to_col    <- setNames(1 + seq_along(groups_sorted), groups_sorted)

## six comparison groups; Personal Care first (primary highlight)
grps_fig   <- c("39","33","37","31","41","43")
labels_fig <- c("39"="Personal Care and Service","33"="Protective Service",
                "37"="Building Maintenance","31"="Healthcare Support",
                "41"="Sales and Related","43"="Office and Admin Support")
level_order <- c("Personal Care and Service","Protective Service","Building Maintenance",
                 "Healthcare Support","Sales and Related","Office and Admin Support")

## Extract columns for the 6 comparison groups (vs management baseline, SOC 11)
cols_fig <- ind_to_col[grps_fig]
cos_fig  <- cos_sim[, c(1, unname(cols_fig)), drop = FALSE]
names(cos_fig) <- c("year", grps_fig)

cos_fig_long <- cos_fig %>%
  pivot_longer(-year, names_to = "ind", values_to = "value") %>%
  mutate(variable = factor(labels_fig[ind], levels = level_order))

ggplot(cos_fig_long,
       aes(x = year, y = value, group = variable, color = variable,
           fill = variable, shape = variable, lty = variable)) +
  geom_point(aes(size = variable)) +
  geom_line() +
  scale_x_continuous(limits = c(1999, 2020),breaks = seq(1999, 2020, 3)) +
  scale_y_continuous(limits = c(0.15, 1.45), breaks = seq(0.2, 1.4, 0.2)) +
  theme_classic() +
  scale_colour_manual(values = c("#82243b","grey30","grey30","#82243b","#ed5278","grey30")) +
  scale_fill_manual(  values = c("#82243b","#ed5278","#ffb6c1",NA,NA,NA)) +
  scale_shape_manual( values = rep(21, 6)) +
  scale_linetype_manual(values = c(2,3,3,3,3,3)) +
  scale_size_manual(  values = rep(3.5, 6)) +
  ylab("cosine similarity") +
  ggtitle("OEWS, 1999–2021") +
  theme(
    axis.text.x     = element_text(size = 14),
    axis.text.y     = element_text(size = 14, angle = 90, hjust = 0.5),
    axis.title.x    = element_text(size = 16, hjust = 1, vjust = 8, margin = margin(t = -12)),
    axis.title.y    = element_text(size = 16, angle = 90),
    plot.title      = element_text(size = 18, hjust = 0.5),
    legend.title    = element_blank(),
    legend.text     = element_text(size = 16),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(ncol = 1), size = guide_legend(direction = "vertical"))

ggsave(file.path(fig_dir, "OEWS_M_base_3.png"),
       width = 11.5, height = 16, units = "cm", dpi = 800)
ggsave(file.path(fig_dir, "OEWS_M_base_3.tiff"),
       device = ragg::agg_tiff,
       width = 11.5, height = 16, units = "cm", dpi = 800, compression = "lzw")
cat("Saved OEWS_M_base_3.png\n")

cos_sim_full <- cos_sim   ## save full-sample result


## ── Section 2: OEWS location LOO (SI S4 Fig) ──────────────────────────────────
## → Output: Analysis/figures/Cosine Similarity/LOO_location_OEWS.png
##           (SI Appendix S4)
if (FALSE) {
pmsa <- intersect(
  unique(oews$metarea[oews$year %in% 1999]),
  unique(oews$metarea[oews$year %in% 2019])
)
cos_sim_pmsa_dropped <- data.frame()

msa_dir <- file.path(data_dir, "misc_data/msa_shapefile")
pmsa99  <- read_sf(file.path(msa_dir, "pm_sa_99_shp/pma_us99.shp")) %>%
  st_set_crs(4326) %>% st_transform(crs = 29101)
pmsa99$area <- st_area(pmsa99)

for (l in seq_along(pmsa)) {
  ## reset cos_sim for this LOO iteration so cosine_similarity() fills it correctly
  cos_sim <- data.frame(year = 1999:2021)

  for (i in 1999:2021) {
    occvec <- create_occvec_owes(data = oews %>% filter(metarea != pmsa[l]),
                                 nu = 50, crosswalk = occind)
    occvec <- occvec_centroid(data = occvec)
    cosine_similarity(occvec, category = 21, standardize = TRUE)
  }

  cos_sim$pmsa <- pmsa[l]
  cos_sim_pmsa_dropped <- rbind(cos_sim_pmsa_dropped, cos_sim)
  cat(sprintf("LOO %d/%d done\n", l, length(pmsa)))
}

## column 15 = baseline (pos 1, SOC 11+13) vs personal care (pos 14, SOC 39): 1+14=15
change_pmsa <- cos_sim_pmsa_dropped[, c(1, 15, ncol(cos_sim_pmsa_dropped))] %>%
  setNames(c("year", "cos_pm_ps", "pmsa")) %>%
  group_by(pmsa) %>%
  summarize(cos_change = cos_pm_ps[year == 2021] - cos_pm_ps[year == 1999],
            .groups = "drop")

## subtract full-sample change as baseline
full_change <- cos_sim_full[cos_sim_full$year == 2021, 15] -
               cos_sim_full[cos_sim_full$year == 1999, 15]
change_pmsa$cos_change <- change_pmsa$cos_change - full_change

density_pmsa <- merge(
  oews %>% group_by(metarea) %>% summarize(tot_emp = sum(tot_emp, na.rm = TRUE)),
  pmsa99[, c("area","P_MSA","NAME")] %>% st_drop_geometry(),
  by.x = "metarea", by.y = "P_MSA", all.x = TRUE
) %>% mutate(area = as.numeric(area), density = tot_emp / area * 10^6)

change_pmsa <- merge(density_pmsa, change_pmsa, by.x = "metarea", by.y = "pmsa", all.y = TRUE)
change_pmsa <- change_pmsa %>%
  mutate(
    dominant = case_when(
      ((cos_change <= -0.011 & log(density) > 6.5) | (NAME == "Las Vegas, NV-AZ MSA")) ~ "dominant",
      TRUE ~ "others"
    ),
    dominant = factor(dominant, levels = c("dominant","others"))
  )

make_panel_oews <- function(data, xvar, xlabel, label_var = "dominant", title = NULL) {
  data <- data %>%
    mutate(NAME = stringr::str_remove(NAME, " MSA"),
           NAME = stringr::str_remove(NAME, " PMSA"),
           x    = log(.data[[xvar]]))
  lm_fit  <- lm(cos_change ~ x, data = data, weight = tot_emp)
  lm_coef <- coef(lm_fit)[2]
  lm_pval <- summary(lm_fit)$coefficients[2, 4]
  plab    <- if (lm_pval < 0.001) "p < 0.001" else sprintf("p = %.3f", lm_pval)
  reg_lab <- sprintf("b = %.4f\n%s", lm_coef, plab)

  ggplot(data, aes(x = x, y = cos_change)) +
    geom_point(aes(size = tot_emp, shape = .data[[label_var]],
                   color = .data[[label_var]], alpha = .data[[label_var]])) +
    geom_smooth(method = "loess", se = FALSE, color = "#82243b",
                aes(weight = log(tot_emp))) +
    annotate("text", x = max(data$x, na.rm = TRUE), y = 0.04,
             label = reg_lab, hjust = 1, vjust = 1, size = 3.5, color = "grey20") +
    scale_shape_manual(values = c(dominant = 5, others = 1)) +
    scale_color_manual(values = c(dominant = "#82243b", others = "grey30")) +
    scale_alpha_manual(values = c(dominant = 1, others = 0.5)) +
    scale_size_continuous(range = c(1, 3)) +
    scale_y_continuous(limits = c(-0.06, 0.04), breaks = seq(-0.06, 0.04, 0.02)) +
    ylab("cosine change") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    xlab(xlabel) + ggtitle(title) +
    geom_text_repel(
      data = ~ subset(.x, .x[[label_var]] == "dominant"),
      aes(label = NAME),
      nudge_y = 0.01, size = 4, color = "#82243b",
      box.padding = 0.35, point.padding = 0.5,
      segment.color = "#82243b", max.overlaps = Inf, force = 8, show.legend = FALSE
    ) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12, angle = 90, hjust = 0.5),
      axis.title.x = element_text(size = 14, hjust = 1, vjust = 8, margin = margin(t = -10)),
      axis.title.y = element_text(size = 14, angle = 90),
      plot.title   = element_text(size = 13, hjust = 0.5),
      panel.grid   = element_blank(), legend.position = "none"
    )
}

q_inf  <- quantile(change_pmsa$cos_change, 0.05, na.rm = TRUE)
q_size <- quantile(change_pmsa$tot_emp,    0.75, na.rm = TRUE)
change_pmsa <- change_pmsa %>%
  mutate(dominant_size = factor(
    ifelse((cos_change <= q_inf & tot_emp >= q_size) |
             grepl("Las Vegas", NAME, ignore.case = TRUE), "dominant", "others"),
    levels = c("dominant","others")))

p_density <- make_panel_oews(change_pmsa, "density", "labor force density (log)",
                              label_var = "dominant",      title = "OEWS, 1999–2021")
p_size    <- make_panel_oews(change_pmsa, "tot_emp",  "labor force size (log)",
                              label_var = "dominant_size", title = "OEWS, 1999–2021")

ggsave(file.path(fig_dir, "LOO_location_OEWS.png"),
       plot = arrangeGrob(p_density, p_size, nrow = 1, ncol = 2),
       width = 22, height = 11, units = "cm", dpi = 800)
cat("Saved LOO_location_OEWS.png\n")
}  ## end if(FALSE) — Section 2


## ── Section 3: Table 1 OEWS column (Panel A lead-lag + Panel B Bartik IV) ──────
## → Built entirely from scratch using standardized OEWS from 01_data.R
## → PM=SOC 11+13, 3yr pooled SVD, nu=20, 6-window stacked lead-lag, Bartik IV

pm_ind  <- c("11","13")
ps_ind  <- "39"
nu_o    <- 20
vec_o   <- paste0("V", 1:nu_o)
oews_t0 <- 1999:2001
oews_t1 <- 2019:2021

## Standardize OEWS area codes to PMSA 1999 definitions
source(file.path(code_dir, "04a_standardize_OEWS_areas.R"))  ## spatial crosswalk for OEWS area codes
oews_b <- oews %>%
  filter(!is.na(occ_code), !is.na(metarea), !is.na(tot_emp), tot_emp > 0) %>%
  mutate(ind = substr(occ_code, 1, 2))

## ── Build PMSA-year panel (3yr pooled SVD) ───────────────────────────────────
all_years_o  <- sort(unique(oews_b$year))
panel_list_o <- list()

for (yr in all_years_o) {
  pool_yrs <- (yr-1):(yr+1); pool_yrs <- pool_yrs[pool_yrs %in% all_years_o]
  sub <- oews_b %>% filter(year %in% pool_yrs) %>%
    group_by(metarea, occ_code) %>%
    summarize(emp = sum(tot_emp), .groups = "drop") %>%
    mutate(ind = substr(occ_code, 1, 2))

  olm <- questionr::wtd.table(sub$occ_code, sub$metarea, sub$emp) %>%
    prop.table(margin = 2)
  olm <- scale(olm, center = TRUE, scale = TRUE); olm[is.na(olm)] <- 0

  sv       <- svd(olm, nu = nu_o, nv = nu_o)
  occ_vecs <- t(apply(sv$u %*% diag(sv$d[1:nu_o]), 1, l2.norm))
  rownames(occ_vecs) <- rownames(olm)
  occ_df   <- as.data.frame(occ_vecs); colnames(occ_df) <- vec_o
  occ_df$occ_code <- rownames(occ_vecs); occ_df$ind <- substr(occ_df$occ_code, 1, 2)

  area_occ  <- sub %>% left_join(occ_df %>% select(-ind), by = "occ_code")
  centroids <- area_occ %>%
    filter(ind %in% c(pm_ind, ps_ind)) %>%
    mutate(grp = if_else(ind %in% pm_ind, "PM", "PS")) %>%
    group_by(metarea, grp) %>%
    summarize(across(all_of(vec_o), ~weighted.mean(., w = emp, na.rm = TRUE)),
              emp_total = sum(emp), n_occ = n(), .groups = "drop") %>%
    filter(n_occ >= 2)

  pm_ct  <- centroids %>% filter(grp == "PM", emp_total >= 5) %>% arrange(metarea)
  ps_ct  <- centroids %>% filter(grp == "PS", emp_total >= 3) %>% arrange(metarea)
  common <- intersect(pm_ct$metarea, ps_ct$metarea)
  if (length(common) < 10) next
  pm_ct  <- pm_ct %>% filter(metarea %in% common) %>% arrange(metarea)
  ps_ct  <- ps_ct %>% filter(metarea %in% common) %>% arrange(metarea)

  pm_mat   <- as.matrix(pm_ct[, vec_o]); ps_mat <- as.matrix(ps_ct[, vec_o])
  cos_vals <- rowSums(pm_mat * ps_mat) / (sqrt(rowSums(pm_mat^2)) * sqrt(rowSums(ps_mat^2)))

  all_emp <- area_occ %>% filter(metarea %in% common) %>% group_by(metarea) %>%
    summarize(total_emp = sum(emp), pm_emp = sum(emp[ind %in% pm_ind]),
              ps_emp = sum(emp[ind == ps_ind]), .groups = "drop") %>% arrange(metarea)

  panel_list_o[[length(panel_list_o) + 1]] <- data.frame(
    metarea  = common, year = yr, cos_sim = cos_vals,
    pm_share = all_emp$pm_emp / all_emp$total_emp,
    ps_share = all_emp$ps_emp / all_emp$total_emp)
}

oews_panel <- bind_rows(panel_list_o) %>%
  mutate(metarea = as.character(metarea), year = as.integer(year)) %>%
  arrange(metarea, year)

## ── Panel A: 6-window stacked lead-lag ───────────────────────────────────────
bal_ids  <- oews_panel %>% dplyr::count(metarea) %>%
  filter(n == n_distinct(oews_panel$year)) %>% pull(metarea)
oews_bal <- oews_panel %>% filter(metarea %in% bal_ids) %>% arrange(metarea, year)

all_yrs_o  <- sort(unique(oews_bal$year))
nwin_o     <- 6
win_size_o <- floor(length(all_yrs_o) / nwin_o)
windows_o  <- vector("list", nwin_o)
for (i in 1:(nwin_o-1)) windows_o[[i]] <- all_yrs_o[((i-1)*win_size_o+1):(i*win_size_o)]
windows_o[[nwin_o]] <- all_yrs_o[((nwin_o-1)*win_size_o+1):length(all_yrs_o)]

stacked_oews <- do.call(rbind, lapply(1:nwin_o, function(i) {
  oews_bal %>% filter(year %in% windows_o[[i]]) %>%
    group_by(metarea) %>%
    summarize(cos_sim  = mean(cos_sim), pm_share = mean(pm_share),
              ps_share = mean(ps_share), .groups = "drop") %>%
    mutate(window = i)
})) %>% arrange(metarea, window) %>% group_by(metarea) %>%
  mutate(d_cos = cos_sim - lag(cos_sim), d_pm = pm_share - lag(pm_share),
         d_ps  = ps_share - lag(ps_share), l_d_pm = lag(d_pm), l_d_ps = lag(d_ps)) %>%
  ungroup()

oews_ll_pm <- feols(d_cos ~ l_d_pm | metarea + window, data = stacked_oews, vcov = ~metarea)
oews_ll_ps <- feols(d_cos ~ l_d_ps | metarea + window, data = stacked_oews, vcov = ~metarea)
oews_n_ll  <- n_distinct(stacked_oews$metarea[!is.na(stacked_oews$l_d_pm)])
oews_ll    <- oews_ll_pm

## ── Panel B: Bartik shift-share IV ───────────────────────────────────────────
## d_pm and bartik from oews_b; d_cos from oews_panel built above

oews_s <- oews_b %>%
  filter(year %in% c(oews_t0, oews_t1)) %>%
  mutate(period = if_else(year %in% oews_t0, "t0", "t1")) %>%
  group_by(metarea, period) %>%
  summarize(pm = sum(tot_emp[ind %in% pm_ind]), tot = sum(tot_emp), .groups = "drop") %>%
  mutate(pm_sh = pm / tot) %>%
  select(metarea, period, pm_sh) %>%
  pivot_wider(names_from = period, values_from = pm_sh) %>%
  filter(!is.na(t0), !is.na(t1)) %>%
  mutate(d_pm = t1 - t0)

pm0_occ <- oews_b %>% filter(ind %in% pm_ind, year %in% oews_t0) %>%
  group_by(metarea, occ_code) %>% summarize(e0 = sum(tot_emp), .groups = "drop") %>%
  group_by(metarea) %>% mutate(sh0 = e0 / sum(e0)) %>% ungroup()
pm1_occ <- oews_b %>% filter(ind %in% pm_ind, year %in% oews_t1) %>%
  group_by(metarea, occ_code) %>% summarize(e1 = sum(tot_emp), .groups = "drop")
pm_nat  <- pm0_occ %>% left_join(pm1_occ, by = c("metarea","occ_code")) %>%
  mutate(e1 = if_else(is.na(e1), 0, e1)) %>%
  group_by(occ_code) %>% summarize(n0 = sum(e0), n1 = sum(e1), .groups = "drop")
bartik_oews <- pm0_occ %>%
  left_join(pm1_occ, by = c("metarea","occ_code")) %>%
  mutate(e1 = if_else(is.na(e1), 0, e1)) %>%
  left_join(pm_nat, by = "occ_code") %>%
  mutate(loo_g = ((n1 - e1) - (n0 - e0)) / (n0 - e0)) %>%
  filter(is.finite(loo_g)) %>%
  group_by(metarea) %>%
  summarize(bartik_pm = sum(sh0 * loo_g, na.rm = TRUE), .groups = "drop")

oews_cos_ld <- oews_panel %>%
  filter(year %in% c(oews_t0, oews_t1)) %>%
  mutate(period = if_else(year %in% oews_t0, "t0", "t1")) %>%
  group_by(metarea, period) %>% summarize(cos = mean(cos_sim), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = cos) %>%
  filter(!is.na(t0), !is.na(t1)) %>% mutate(d_cos = t1 - t0)

oews_iv_df <- oews_s %>%
  inner_join(oews_cos_ld %>% select(metarea, d_cos), by = "metarea") %>%
  inner_join(bartik_oews, by = "metarea")

oews_fs  <- feols(d_pm ~ bartik_pm,             data = oews_iv_df, vcov = "HC1")
oews_ols <- feols(d_cos ~ d_pm,                 data = oews_iv_df, vcov = "HC1")
oews_iv  <- feols(d_cos ~ 1 | d_pm ~ bartik_pm, data = oews_iv_df, vcov = "HC1")

cat("\n=== TABLE 1 OEWS RESULTS ===\n")
cat("Panel A PM:\n"); print(oews_ll_pm)
cat("\nPanel A PS:\n"); print(oews_ll_ps)
cat("\nPanel B OLS:\n"); print(oews_ols)
cat("\nPanel B IV:\n"); print(oews_iv)
cat(sprintf("\nF: %.1f  |  Panel A N: %d  |  Panel B N: %d\n",
            fitstat(oews_fs,"f")$f$stat, oews_n_ll, nrow(oews_iv_df)))

## ── Clean up intermediate objects ────────────────────────────────────────────
rm(occind, cos_sim, cos_sim_full, cos_fig, cos_fig_long,
   groups_sorted, ind_to_col, grps_fig, labels_fig, level_order, cols_fig,
   oews_b, panel_list_o, oews_panel, oews_bal, bal_ids, stacked_oews,
   oews_s, pm0_occ, pm1_occ, pm_nat, bartik_oews, oews_cos_ld, oews_iv_df,
   pm_ind, ps_ind, nu_o, vec_o, oews_t0, oews_t1,
   all_years_o, all_yrs_o, win_size_o, windows_o, nwin_o,
   oews_ll_pm, oews_ll_ps, oews_ols, oews_iv, oews_fs, oews_ll, oews_n_ll)
gc()
