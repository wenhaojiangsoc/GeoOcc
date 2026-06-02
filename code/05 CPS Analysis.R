## 05 CPS Analysis.R
## All CPS analysis: CPS DA cosine similarity (SI S2), CPS EGP cosine similarity (SI S2),
## CPS location LOO plot (SI S4), and CPS occupation OLS (SI S5).
##
## Sources: 05_CPS_cosine.R, 06_CPS_egp.R, 08_CPS_OEWS_LOO_plot.R,
##          10_S5_CPS_OEWS_occ_ols.R (CPS parts)

repo_root <- "~/Dropbox/GeoOcc/GeoOccGit"
code_dir  <- file.path(repo_root, "code")
source(file.path(code_dir, "01 data.R"))
source(file.path(code_dir, "02 functions.R"))
select <- dplyr::select; filter <- dplyr::filter

library(reshape2)
library(ggrepel)
library(sf)
library(gridExtra)
library(grid)
library(patchwork)
library(scales)

fig_dir <- file.path(repo_root, "figures", "si", "Cosine Similarity")


## ── Section 1: CPS DA cosine similarity (SI S2) ───────────────────────────────
## → Output: Analysis/figures/Cosine Similarity/CPS_DA_PM_base.png
##           (SI Appendix S2, CPS Autor-Dorn figure)

## simplify CPS
cps_agg <- cps %>%
  group_by(occ1990, year, egp, da, metarea) %>%
  summarize(wtfinl=sum(wtfinl,na.rm=T))

## recode barbers and hairdressers
occegp[which(occegp$occ1990==457),"egp"] <- "VII"
occegp[which(occegp$occ1990==458),"egp"] <- "VII"
cps_agg[which(cps_agg$occ1990==457),"egp"] <- "VII"
cps_agg[which(cps_agg$occ1990==458),"egp"] <- "VII"

## empty df to store cosine similarity in each year
cos_sim <-
  data.frame(year=1994:2022)

for (i in 1994:2022){

  ## create occupation vector
  occvec <-
    create_occvec(data=cps_agg %>% filter(!is.na(da)),
                  nu=50,
                  crosswalk=occda,
                  cps=T)

  ## calculate centroid and weighted by individual counts
  occvec <- occvec_centroid(data = occvec)

  ## calculate cosine similarity
  cosine_similarity(occvec,category=8,standardize=TRUE)

  ## monitor progress
  print(paste("year", as.character(i), "is done!"))
}

## plot
cos_sim_PM_base <- cos_sim[, c(1,3:9)]
for (col in colnames(cos_sim_PM_base)[2:8])
  cos_sim_PM_base[, col] <- three_ma(cos_sim_PM_base[, col])
cos_sim_PM_base <- cos_sim_PM_base %>%
  pivot_longer(cols = -year, names_to = "variable", values_to = "value")
cos_sim_PM_base$variable <- as.character(cos_sim_PM_base$variable)
cos_sim_PM_base[cos_sim_PM_base$variable=="V3","variable"] <- "Technicians"
cos_sim_PM_base[cos_sim_PM_base$variable=="V4","variable"] <- "Sales Occupations"
cos_sim_PM_base[cos_sim_PM_base$variable=="V5","variable"] <- "Administrative and Office"
cos_sim_PM_base[cos_sim_PM_base$variable=="V6","variable"] <- "Production"
cos_sim_PM_base[cos_sim_PM_base$variable=="V7","variable"] <- "Laborers"
cos_sim_PM_base[cos_sim_PM_base$variable=="V8","variable"] <- "Clean and Protective"
cos_sim_PM_base[cos_sim_PM_base$variable=="V9","variable"] <- "Personal Services"

cos_sim_PM_base <- cos_sim_PM_base %>%
  filter(variable %in% c("Personal Services","Clean and Protective",
                          "Administrative and Office","Sales Occupations",
                          "Laborers","Production")) %>%
  mutate(variable = factor(variable, levels = c("Personal Services","Clean and Protective",
                                                "Administrative and Office","Sales Occupations",
                                                "Laborers","Production")),
         value = if_else(variable == "Production", value + 1.5, value))

p_s1 <- ggplot(cos_sim_PM_base,
       aes(x=year, y=value, group=variable, color=variable, fill=variable,
           shape=variable, lty=variable)) +
  geom_point(size=3.5) +
  geom_line() +
  scale_x_continuous(breaks = seq(1994, 2022, 4)) +
  theme_classic() +
  scale_colour_manual(values = c("Personal Services"="#82243b","Clean and Protective"="grey30",
                                 "Administrative and Office"="grey30","Sales Occupations"="#82243b",
                                 "Laborers"="#ed5278","Production"="grey60")) +
  scale_fill_manual(values = c("Personal Services"="#82243b","Clean and Protective"="#ed5278",
                                "Administrative and Office"="#ffb6c1","Sales Occupations"=NA,
                                "Laborers"=NA,"Production"=NA), na.value = NA) +
  scale_shape_manual(values = rep(21, 6)) +
  scale_linetype_manual(values = c("Personal Services"=2,"Clean and Protective"=3,
                                   "Administrative and Office"=3,"Sales Occupations"=3,
                                   "Laborers"=3,"Production"=4)) +
  ylab("cosine similarity") +
  ggtitle("CPS, 1994-2022") +
  theme(axis.text.x     = element_text(size=14),
        axis.text.y     = element_text(size=14, angle=90, hjust=0.5),
        axis.title.x    = element_text(size=16, hjust=1, vjust=8, margin=margin(t=-12)),
        axis.title.y    = element_text(size=16, angle=90),
        plot.title      = element_text(size=18, hjust=0.5),
        legend.title    = element_blank(),
        legend.text     = element_text(size=16),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.key.width = unit(1.2, "cm")) +
  guides(color=guide_legend(ncol=1))

## save with fixed legend height (consistent with Section 2 save_fixed_c style)
get_leg <- function(p) {
  gt <- ggplot_gtable(ggplot_build(p))
  gt$grobs[[which(sapply(gt$grobs, function(x) x$name) == "guide-box")]]
}
p_s1_combined <- arrangeGrob(p_s1 + theme(legend.position="none"),
                              get_leg(p_s1),
                              nrow=2, heights=unit(c(10, 6), "cm"))
ggsave(file.path(fig_dir, "CPS_DA_PM_base.png"),
       plot=p_s1_combined, width=11.5, height=16, units="cm", dpi=800)
cat("Saved CPS_DA_PM_base.png\n")


## ── Section 2: CPS EGP cosine similarity (SI S2) ──────────────────────────────
## → Output: Analysis/figures/Cosine Similarity/CPS_EGP_PM_base.png
##           Analysis/figures/Cosine Similarity/CPS_EGP_DA_combined.png
##
## Uses same create_occvec/cosine_similarity pipeline as ACS Section 4.

## save DA cos_sim from Section 1 before overwriting
cos_sim_da <- cos_sim

## collapse EGP sub-classes in cps_agg (reuses data from Section 1)
cps_agg[cps_agg$egp %in% c("IIIa","IIIb"), "egp"] <- "III"
cps_agg[cps_agg$egp %in% c("VIIa"),         "egp"] <- "VII"
egp_keep <- c("I","II","III","V","VI","VII")

cos_sim <- data.frame(year = 1994:2022)
for (i in 1994:2022) {
  occvec <- create_occvec(data = cps_agg[which(cps_agg$egp %in% egp_keep), ],
                          nu = 50, crosswalk = occegp, cps = TRUE)
  occvec <- occvec_centroid(data = occvec)
  cosine_similarity(occvec, category = 6, standardize = TRUE)
  cat(sprintf("EGP year %d done\n", i))
}

## smooth and label: V3=I_II, V4=I_III, V5=I_V, V6=I_VI, V7=I_VII
for (col in colnames(cos_sim)[2:ncol(cos_sim)])
  cos_sim[, col] <- three_ma(cos_sim[, col])

egp_level_order <- c("Class VII (Semi+Unskilled)", "Class II (Lower Service)",
                     "Class III (Routine Non-Manual)", "Class VI (Skilled Manual)",
                     "Class V (Manual Supervisors)")
egp_df <- cos_sim[, c(1, 3:7)] %>%
  pivot_longer(-year, names_to = "pair", values_to = "value") %>%
  mutate(label = recode(pair,
    V3 = "Class II (Lower Service)",    V4 = "Class III (Routine Non-Manual)",
    V5 = "Class V (Manual Supervisors)", V6 = "Class VI (Skilled Manual)",
    V7 = "Class VII (Semi+Unskilled)"),
    label = factor(label, levels = egp_level_order),
    value = case_when(
      label == "Class V (Manual Supervisors)" ~ value + 1.5,
      label == "Class VI (Skilled Manual)"    ~ value + 1.2,
      TRUE ~ value
    ))

egp_colors <- c("Class VII (Semi+Unskilled)"    = "#82243b",
                "Class II (Lower Service)"       = "grey30",
                "Class III (Routine Non-Manual)" = "grey30",
                "Class VI (Skilled Manual)"      = "#ed5278",
                "Class V (Manual Supervisors)"   = "grey60")
egp_fills  <- c("Class VII (Semi+Unskilled)"    = "#82243b",
                "Class II (Lower Service)"       = "#ed5278",
                "Class III (Routine Non-Manual)" = "#ffb6c1",
                "Class VI (Skilled Manual)"      = NA,
                "Class V (Manual Supervisors)"   = NA)
egp_ltypes <- c("Class VII (Semi+Unskilled)"    = 2,
                "Class II (Lower Service)"       = 3,
                "Class III (Routine Non-Manual)" = 3,
                "Class VI (Skilled Manual)"      = 3,
                "Class V (Manual Supervisors)"   = 4)
lvls_egp <- levels(egp_df$label)

## EGP standalone figure (same layout as Section 1 DA figure)
p_egp <- ggplot(egp_df,
       aes(x = year, y = value, group = label,
           color = label, fill = label, shape = label, lty = label)) +
  geom_point(size = 3.5) + geom_line() +
  scale_color_manual(values = egp_colors) +
  scale_fill_manual(values = egp_fills, na.value = NA) +
  scale_shape_manual(values = setNames(rep(21L, length(lvls_egp)), lvls_egp)) +
  scale_linetype_manual(values = egp_ltypes) +
  scale_x_continuous(breaks = seq(1994, 2022, 4)) +
  ylab("cosine similarity (z-standardized)") +
  ggtitle("CPS, 1994-2022 (EGP)") +
  theme_classic() +
  theme(axis.text.x      = element_text(size = 14),
        axis.text.y      = element_text(size = 14, angle = 90, hjust = 0.5),
        axis.title.x     = element_text(size = 16, hjust = 1, vjust = 8, margin = margin(t = -12)),
        axis.title.y     = element_text(size = 16, angle = 90),
        plot.title       = element_text(size = 18, hjust = 0.5),
        legend.title     = element_blank(),
        legend.text      = element_text(size = 16),
        legend.position  = "bottom",
        legend.direction = "vertical",
        legend.key.width = unit(1.2, "cm")) +
  guides(color = guide_legend(ncol = 1))

p_egp_combined <- arrangeGrob(p_egp + theme(legend.position = "none"),
                               get_leg(p_egp),
                               nrow = 2, heights = unit(c(10, 6), "cm"))
ggsave(file.path(fig_dir, "CPS_EGP_PM_base.png"),
       plot = p_egp_combined, width = 11.5, height = 16, units = "cm", dpi = 800)
cat("Saved CPS_EGP_PM_base.png\n")

## combined EGP + DA figure (2 panels side by side, shared legend row)
p_egp_bare <- p_egp + theme(legend.position = "none")
p_da_bare  <- p_s1  + theme(legend.position = "none")
legend_egp <- get_leg(p_egp)
legend_da  <- get_leg(p_s1)
combined_c <- arrangeGrob(
  arrangeGrob(p_egp_bare, p_da_bare, ncol = 2),
  arrangeGrob(legend_egp, legend_da, ncol = 2),
  nrow = 2, heights = unit(c(10, 6), "cm"))
ggsave(file.path(fig_dir, "CPS_EGP_DA_combined.png"),
       plot = combined_c, width = 22, height = 16, units = "cm", dpi = 800)
cat("Saved CPS_EGP_DA_combined.png\n")


## ── Section 3: CPS location LOO plot (SI S4) ──────────────────────────────────
## → Output: (see Analysis/code/08_CPS_OEWS_LOO_plot.R — uses pre-saved RDS files)
##           data/intermediate/s4_loo_results.rds must exist (see README.md)
##
## The LOO compute step is slow; this section loads the pre-saved results.

loo_rds <- file.path(data_dir, "intermediate/s4_loo_results.rds")
if (!file.exists(loo_rds)) {
  cat("SI S4 CPS LOO: RDS not found at", loo_rds, "\n")
  cat("data/intermediate/s4_loo_results.rds not found — see README.md for how to generate it.\n")
} else {
  res      <- readRDS(loo_rds)
  cps_loo  <- res$cps_loo

  ## load PMSA shapefile for names
  msa_dir  <- file.path(misc_dir, "msa_shapefile")
  pmsa_shp <- read_sf(file.path(msa_dir, "pm_sa_99_shp/pma_us99.shp")) %>%
    st_set_crs(4326) %>% st_transform(5070)
  pmsa_meta <- pmsa_shp %>%
    mutate(pmsa_code = as.character(P_MSA),
           area_sqmi = as.numeric(st_area(geometry)) / 2589988,
           short_name = gsub(",.*", "", NAME, useBytes = TRUE)) %>%
    as.data.frame() %>%
    dplyr::select(pmsa_code, short_name, area_sqmi)

  ## CPS LOO plot (density scatter, analogous to OEWS)
  cps_loo_plot <- cps_loo %>%
    left_join(pmsa_meta, by = c("metarea" = "pmsa_code")) %>%
    mutate(
      dominant = factor(ifelse(cos_change <= quantile(cos_change, 0.05, na.rm=TRUE), "dominant", "others"),
                        levels = c("dominant","others"))
    )

  ggplot(cps_loo_plot, aes(x = log(density), y = cos_change)) +
    geom_point(aes(size = tot_emp, shape = dominant, color = dominant, alpha = dominant)) +
    geom_smooth(method = "loess", se = FALSE, color = "#82243b", aes(weight = log(tot_emp))) +
    scale_shape_manual(values = c(dominant = 5, others = 1)) +
    scale_color_manual(values = c(dominant = "#82243b", others = "grey30")) +
    scale_alpha_manual(values = c(dominant = 1, others = 0.5)) +
    scale_size_continuous(range = c(1, 3)) +
    geom_text_repel(
      data = ~ subset(.x, dominant == "dominant"),
      aes(label = short_name), size = 4, color = "#82243b",
      box.padding = 0.35, point.padding = 0.5, segment.color = "#82243b",
      max.overlaps = Inf, force = 8, show.legend = FALSE
    ) +
    ylab("cosine change") + xlab("labor force density (log)") +
    theme_classic() +
    theme(axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12, angle=90, hjust=0.5),
          axis.title.x = element_text(size=14, hjust=1, vjust=8, margin=margin(t=-10)),
          axis.title.y = element_text(size=14, angle=90),
          legend.position = "none")

  ggsave(file.path(fig_dir, "LOO_location_CPS.png"),
         width = 12, height = 11, units = "cm", dpi = 800)
  cat("Saved LOO_location_CPS.png\n")
}


## ── Section 4: CPS occupation OLS (SI S5) ─────────────────────────────────────
## → Output: Analysis/figures/Cosine Similarity/CPS_occ_ols.png
##           (SI Appendix S5, CPS occupation OLS)
##
## Step-by-step mirror of ACS Section 3:
##   1. Compute occupation vectors once (margin=1 OLM, same SVD approach)
##   2. LOO: drop each occupation from vectors, recompute centroids + cos_sim
##   3. Compute change in affinity attributable to each occupation
##   4. Merge manual O*NET + covariates, run bivariate OLS, plot

library(haven)

## Step 1: local create_occvec (margin=1, keeps occ1990 label) — mirrors ACS Section 3
create_occvec_loo_cps <- function(data, nu = 50, crosswalk = occda) {
  olm <- with(five_rd(data, i), questionr::wtd.table(occ1990, metarea, wtfinl)) %>%
    prop.table(margin = 1)
  olm <- scale(olm, center = TRUE, scale = TRUE)
  Nocc <- data %>% dplyr::group_by(year, occ1990) %>%
    dplyr::summarize(Nocc = sum(wtfinl, na.rm = TRUE), .groups = "drop")
  sv     <- svd(olm, nu = nu)
  occvec <- sv$u %*% diag(sv$d)[1:nu, 1:nu]
  occvec <- t(apply(occvec, 1, l2.norm))
  occvec <- as.data.frame(occvec)
  occvec$occ1990 <- rownames(olm)
  occvec <- merge(occvec, crosswalk, by = "occ1990", all.x = TRUE)
  occvec <- merge(occvec,
                  Nocc[which(Nocc$year %in% five_mw(i)), ] %>%
                    dplyr::group_by(occ1990) %>%
                    dplyr::summarize(Nocc = sum(Nocc, na.rm = TRUE), .groups = "drop"),
                  by = "occ1990", all.x = TRUE)
  return(occvec)
}

cps_loo_data <- cps_agg[which(!is.na(cps_agg$da)), ]
occvecs <- data.frame()
for (i in 1994:2022) {
  ov       <- create_occvec_loo_cps(data = cps_loo_data, nu = 50, crosswalk = occda)
  ov$year  <- i
  occvecs  <- rbind(occvecs, ov)
  cat(sprintf("occvec year %d done\n", i))
}

## Step 2: LOO — drop each occupation, recompute centroids + cos_sim
occ <- unique(cps_loo_data %>% distinct(occ1990) %>% arrange(occ1990) %>% pull(occ1990))
cos_sim_occ_dropped <- data.frame()

for (l in seq_along(occ)) {
  cos_sim <- data.frame(year = 1994:2022)

  for (i in 1994:2022) {
    occvec <- occvecs[which(occvecs$occ1990 != occ[l]), ] %>%
      dplyr::select(-occ1990) %>%
      filter(year == i) %>%
      dplyr::select(-year)
    occvec <- occvec_centroid(data = occvec)
    cosine_similarity(occvec, category = 8, standardize = TRUE)
    cos_sim$occ <- occ[l]
  }

  if (l %% 20 == 0) cat("Finished", l, "of", length(occ), "\n")
  cos_sim_occ_dropped <- rbind(cos_sim_occ_dropped, cos_sim)
}

## Step 3: compute cos_change — mirrors ACS lines 466–520
## with personal service (DA = 1 as baseline)
t <- merge(cos_sim_occ_dropped,
           occda %>% mutate(occ1990 = as.character(occ1990)),
           by.x = "occ", by.y = "occ1990", all.x = TRUE) %>%
  mutate(da = da - 1, da = case_when(da == 0 ~ 1, .default = da)) %>%
  filter(da != "8")

change_occs <- data.frame()
for (g in seq(1, 7)) {
  change_occ <- melt(t[t$da == g, c(1, 2, 1+1+8*g)], id.vars = c("year","occ")) %>%
    group_by(occ) %>%
    summarize(cos_change = value[year == 2022] - value[year == 1994], .groups = "drop")
  change_occ$cos_change <- change_occ$cos_change -
    (cos_sim[, c(1, 1+8*g)] %>%
       summarize(cos_change = .data[[paste0("V", 1+8*g)]][year == 2022] -
                               .data[[paste0("V", 1+8*g)]][year == 1994]) %>%
       pull())
  change_occs <- rbind(change_occs, change_occ)
}
change_occ_pm <- change_occs

## with professional and managerial (DA = 9 as baseline)
t <- merge(cos_sim_occ_dropped,
           occda %>% mutate(occ1990 = as.character(occ1990)),
           by.x = "occ", by.y = "occ1990", all.x = TRUE) %>%
  mutate(da = da - 1, da = case_when(da == 0 ~ 1, .default = da)) %>%
  filter(da != "1")

change_occs <- data.frame()
for (g in seq(2, 8)) {
  change_occ <- melt(t[t$da == g, c(1, 2, 1+1+1+8*(g-1))], id.vars = c("year","occ")) %>%
    group_by(occ) %>%
    summarize(cos_change = value[year == 2022] - value[year == 1994], .groups = "drop")
  change_occ$cos_change <- change_occ$cos_change -
    (cos_sim[, c(1, 2+8*(g-1))] %>%
       summarize(cos_change = .data[[paste0("V", 2+8*(g-1))]][year == 2022] -
                               .data[[paste0("V", 2+8*(g-1))]][year == 1994]) %>%
       pull())
  change_occs <- rbind(change_occs, change_occ)
}
change_occ_service <- change_occs

change_occ <- rbind(
  merge(change_occ_pm, occda %>% mutate(da = as.character(da)),
        by.x = "occ", by.y = "occ1990", all.x = TRUE) %>% filter(da == "1"),
  merge(change_occ_service, occda %>% mutate(da = as.character(da)),
        by.x = "occ", by.y = "occ1990", all.x = TRUE) %>% filter(da == "9")
)

## Step 4: merge covariates + O*NET, run OLS, plot

## manual O*NET coding (same as ACS Section 3)
  onet_2023 <- data.frame(
    occ = c(175, 404, 406, 433, 434, 435, 436, 438, 439, 443, 444, 445, 446,
            447, 456, 457, 458, 459, 461, 462, 464, 466, 468, 469, 487, 773),
    physical_prox = c(72.0, 80.8, 81.7, 95.2, 78.5, 77.8, 80.8, 73.0, 69.2,
                      79.2, 62.7, 98.5, 84.0, 90.8, 84.2, 94.8, 92.3, 83.5,
                      72.0, 81.5, 78.2, 81.7, 81.7, 71.0, 54.8, 36.8),
    face2face     = c(90.8, 75.8, 88.8, 87.0, 80.0, 76.2, 75.8, 90.3, 80.8,
                      86.8, 62.0, 91.5, 80.2, 92.8, 98.0, 82.2, 95.8, 87.8,
                      85.8, 92.3, 78.8, 88.8, 88.8, 94.8, 89.0, 76.2),
    interpersonal = c(79.0, 69.8, 79.0, 65.2, 77.2, 76.5, 69.8, 72.0, 73.2,
                      72.5, 43.8, 86.8, 78.8, 86.8, 78.8, 78.0, 71.8, 68.0,
                      79.2, 79.2, 58.2, 79.0, 79.0, 78.8, 47.8, 35.2),
    ext_customers = c(93.3, 52.8, 59.8, 94.2, 81.2, 72.0, 52.8, 86.5, 80.8,
                      73.0, 45.8, 81.2, 58.2, 58.5, 96.2, 59.2, 87.0, 86.8,
                      87.2, 88.5, 83.7, 59.8, 59.8, 60.0, 76.8, 51.5),
    responsible_health = c(75.2, 27.2, 84.2, 42.5, 53.0, 50.2, 27.2, 37.0,
                           52.0, 57.0, 42.0, 89.0, 85.8, 98.2, 60.5, 58.2,
                           75.5, 58.5, 44.0, 59.0, 40.5, 84.2, 84.2, 70.0,
                           52.2, 22.2),
    stringsAsFactors = FALSE
  )
  onet_2023$occ <- as.character(onet_2023$occ)

  change_occ <- change_occ %>%
    merge(onet_2023, by = "occ", all.x = TRUE) %>%
    merge(covariate %>% mutate(occ1990 = as.character(occ1990)),
          by.x = "occ", by.y = "occ1990", all.x = TRUE) %>%
    merge(read_dta(file.path(misc_dir, "occ1990_titles.dta")) %>%
            mutate(occ1990 = as.character(occ1990)),
          by.x = "occ", by.y = "occ1990", all.x = TRUE)

  ## top contributing occupations lollipop (mirrors ACS_DA_occ_heterogeneity.png)
  change_occ %>%
    filter(da %in% c("1","9")) %>%
    filter(!grepl("Health", title) & !grepl("Nursing", title)) %>%
    group_by(da) %>%
    arrange(cos_change) %>%
    slice_head(n = 9) %>%
    ungroup() %>%
    mutate(title = case_when(
      grepl("Kitchen", title)          ~ "Kitchen worker",
      grepl("admin", title)            ~ "Managers and administrators",
      grepl("Health aides", title)     ~ "Health aides",
      grepl("Nursing aides", title)    ~ "Nursing aides",
      grepl("Personal,", title)        ~ "Labor relation specialists",
      grepl("Personal service", title) ~ "Personal service n.e.c",
      grepl("Attendants, amusements", title) ~ "Attendants, amusements/rec.",
      TRUE                             ~ title
    )) %>%
    group_by(da) %>%
    mutate(title_plot = factor(title, levels = unique(title))) %>%
    ungroup() %>%
    ggplot(aes(x = title_plot, y = cos_change)) +
    geom_segment(aes(x = title_plot, xend = title_plot, y = 0, yend = cos_change),
                 color = "grey60", linewidth = 0.6) +
    geom_point(aes(fill = factor(da)), size = 3.7, shape = 21) +
    coord_flip() +
    scale_y_reverse() +
    scale_x_discrete(labels = function(x) x) +
    scale_fill_manual(values = c("1" = "#82243b", "9" = "#ffe4ec"),
                      labels = c("1" = "Professional-Management", "9" = "Personal Service"),
                      name = NULL) +
    facet_grid(da ~ ., scales = "free_y", space = "free") +
    theme_classic() +
    ggtitle("Top Occupations, CPS, 1994-2022") +
    labs(x = "", y = "cosine change") +
    theme(axis.text.x     = element_text(size = 12),
          axis.text.y     = element_text(size = 12, hjust = 1),
          axis.title      = element_text(size = 13),
          legend.text     = element_text(size = 12),
          legend.position = "bottom",
          plot.title      = element_text(hjust = 0.5),
          strip.text      = element_blank(),
          strip.background = element_blank()) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))

  ggsave(file.path(fig_dir, "CPS_occ_heterogeneity.png"),
         width = 16, height = 17.5, units = "cm", dpi = 800)
  cat("Saved CPS_occ_heterogeneity.png\n")

  change_occ <- change_occ %>%
    group_by(da) %>%
    mutate(
      cos_change_z      = scale(cos_change)[,1],
      log_incwage_z     = scale(log(incwage))[,1],
      overwork_male_z   = scale(overwork_male)[,1],
      overwork_female_z = scale(overwork_female)[,1],
      fulltime_z        = scale(1 - parttime_share)[,1],
      uhrswork_male_z   = scale(uhrswork_male)[,1],
      uhrswork_female_z = scale(uhrswork_female)[,1],
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
    ) %>% ungroup()

  results <- data.frame(
    group   = c(rep("PM", 10), rep("Service", 10)),
    variable = c(
      "Annual income (log)", "Share of overwork (male)", "Share of overwork (female)",
      "    Weekly working hours (male)", "    Weekly working hours (female)",
      "Share of full-time working", "Share of white", "Share of Black",
      "Share of Asian", "Share of Hispanic",
      "Share of white", "Share of Black", "Share of Asian", "Share of Hispanic",
      "Share of immigrants", "Physical proximity", "Face-to-face discussions",
      "Interpersonally oriented", "Serves external clients", "Care responsibility"),
    feature  = c(
      "demand","demand","demand","demand","demand","demand",
      "race/ethnicity","race/ethnicity","race/ethnicity","race/ethnicity",
      "race/ethnicity","race/ethnicity","race/ethnicity","race/ethnicity","",
      "work content/style","work content/style","work content/style",
      "work content/style","work content/style"),
    coef = NA_real_, se = NA_real_
  )

  for (i in 1:10) {
    v <- c("log_incwage_z","overwork_male_z","overwork_female_z",
           "uhrswork_male_z","uhrswork_female_z","fulltime_z",
           "share_white_z","share_black_z","share_asian_z","share_hispanic_z")[i]
    fit <- lm(as.formula(paste("cos_change_z ~", v)), data = change_occ %>% filter(da == 1L))
    s <- summary(fit)
    results$coef[i] <- s$coefficients[2, "Estimate"]
    results$se[i]   <- s$coefficients[2, "Std. Error"]
  }
  for (i in 1:10) {
    v <- c("share_white_z","share_black_z","share_asian_z","share_hispanic_z",
           "share_immigration_z","physical_prox_z","face2face_z",
           "interpersonal_z","ext_customers_z","responsible_z")[i]
    fit <- lm(as.formula(paste("cos_change_z ~", v)),
              data = change_occ %>% filter(da == 9L), weight = weight)
    s <- summary(fit)
    results$coef[10 + i] <- s$coefficients[2, "Estimate"]
    results$se[10 + i]   <- s$coefficients[2, "Std. Error"]
  }

  results %>%
    mutate(lo = coef - 1.96 * se, hi = coef + 1.96 * se) %>%
    {
      pm_df  <- filter(., group == "PM")  %>% mutate(variable = as.character(variable))
      svc_df <- filter(., group == "Service") %>% mutate(variable = as.character(variable))

      x_m_pm  <- max(abs(c(pm_df$lo,  pm_df$hi)),  na.rm = TRUE)
      x_m_svc <- max(abs(c(svc_df$lo, svc_df$hi)), na.rm = TRUE)

      reorder_within_feature <- function(df) {
        df %>% group_by(feature) %>% arrange(desc(coef)) %>%
          mutate(.ord = row_number()) %>% ungroup() %>%
          group_by(feature) %>% mutate(feat_mean = mean(coef, na.rm = TRUE)) %>% ungroup() %>%
          arrange(desc(feat_mean), feature, .ord) %>%
          mutate(variable = factor(variable, levels = unique(variable)))
      }
      pm_df  <- reorder_within_feature(pm_df)
      svc_df <- reorder_within_feature(svc_df)

      sep_y <- function(d) {
        d %>% mutate(.y = as.numeric(variable)) %>%
          group_by(feature) %>% summarize(max_y = max(.y), .groups = "drop") %>%
          arrange(max_y) %>% pull(max_y) %>%
          { if (length(.) > 1) head(., -1) + 0.5 else numeric(0) }
      }
      sep_pm  <- sep_y(pm_df)
      sep_svc <- sep_y(svc_df)

      p_pm <- ggplot(pm_df, aes(x = coef, y = variable)) +
        { if (length(sep_pm) > 0) geom_hline(yintercept = sep_pm, color = "grey85", linewidth = 0.35, inherit.aes = FALSE) } +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
        geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, color = "grey40") +
        geom_point(size = 3.7, shape = 21, fill = "#82243b", color = "grey40") +
        scale_y_discrete(limits = levels(pm_df$variable)) +
        coord_cartesian(xlim = c(-x_m_pm - 0.2, x_m_pm + 0.2), expand = TRUE) +
        theme_classic() +
        labs(x = NULL, y = NULL, title = "Cosine Change (CPS, 1994-2022)\nand Occupational Features") +
        theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12, hjust = 1),
              axis.title = element_text(size = 13), plot.title = element_text(hjust = 0.5, size = 11),
              legend.position = "bottom")

      p_svc <- ggplot(svc_df, aes(x = coef, y = variable)) +
        { if (length(sep_svc) > 0) geom_hline(yintercept = sep_svc, color = "grey85", linewidth = 0.35, inherit.aes = FALSE) } +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
        geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, color = "grey40") +
        geom_point(size = 3.7, shape = 21, fill = "#ffe4ec", color = "grey40") +
        scale_y_discrete(limits = levels(svc_df$variable)) +
        scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
        coord_cartesian(xlim = c(-x_m_svc - 1, x_m_svc + 1), expand = TRUE) +
        theme_classic() +
        geom_point(data = data.frame(group = c("Professional-Management","Personal Service"),
                                     coef = NA, variable = NA),
                   aes(fill = group), shape = 21, size = 3.7, color = "black", show.legend = TRUE) +
        scale_fill_manual(values = c("Personal Service" = "#82243b", "Professional-Management" = "#ffe4ec"),
                          labels = c("Professional-Management","Personal Service"), name = NULL) +
        labs(x = "coefficient (bivariate)", y = NULL, title = NULL) +
        theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12, hjust = 1),
              axis.title = element_text(size = 13), legend.text = element_text(size = 12),
              legend.position = "bottom") +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE))

      p_pm / p_svc + plot_layout(heights = c(nrow(pm_df), nrow(svc_df))) +
        theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0.5, unit = "cm"))
    }

ggsave(file.path(fig_dir, "CPS_occ_ols.png"),
       width = 16, height = 17.5, units = "cm", dpi = 800)
cat("Saved CPS_occ_ols.png\n")

## ── Cleanup ──────────────────────────────────────────────────────────────────
## Remove intermediate plot objects and temporary loop accumulators only.
## Keeping: cps_agg, cos_sim_da, change_occ, results, occvecs, onet_2023
##          (may be used in SI Appendix)
rm(cos_sim, cos_sim_PM_base, cos_sim_occ_dropped,
   egp_df, egp_keep, egp_colors, egp_fills, egp_ltypes, egp_level_order, lvls_egp,
   p_s1, p_s1_combined, p_egp, p_egp_combined, p_egp_bare, p_da_bare,
   legend_egp, legend_da, combined_c,
   change_occ_pm, change_occ_service, change_occs, t,
   create_occvec_loo_cps, get_leg)
