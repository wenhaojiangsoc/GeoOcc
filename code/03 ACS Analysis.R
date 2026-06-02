## 03 ACS Analysis.R
## All ACS/Census analysis: cosine similarity trends, location heterogeneity,
## occupation heterogeneity and OLS, ACS EGP typology (SI S2), ACS tract (SI S6),
## and Table 1 ACS column (Panel A lead-lag + Panel B Bartik IV).
##
## Sources: 03_ACS_main.R, 04_ACS_egp.R, 16_ACS_tract.R, 16_causal_ACS_leadlag.R

repo_root <- "~/Dropbox/GeoOcc/GeoOccGit"
code_dir  <- file.path(repo_root, "code")
source(file.path(code_dir, "01 data.R"))
source(file.path(code_dir, "02 functions.R"))

## Re-export masked functions (MASS and other packages loaded by 02_functions.R
## mask dplyr::select and dplyr::filter — restore them explicitly)
select <- dplyr::select
filter <- dplyr::filter

library(reshape2)
library(lme4)
library(matrixStats)
library(ggrepel)
library(patchwork)
library(rlang)
library(ggbrace)
library(ragg)
library(fixest)
library(tidyr)

fig_dir      <- file.path(repo_root, "figures", "si", "Cosine Similarity")
main_fig_dir <- file.path(repo_root, "figures", "main")

## ── Section 1: Main ACS cosine similarity figure (Fig 3 left panel) ──────────
## → Output: Analysis/figures/Cosine Similarity/ACS_DA_PM_base.png
##           (Main paper Fig 3, left panel)

## simplify ACS
acs_agg <- acs %>%
  group_by(occ1990, year, egp, da, metarea) %>%
  summarize(wtfinl=sum(wtfinl,na.rm=T))

acs_agg[which(acs_agg$year==1980),"year"] <- 2002
acs_agg[which(acs_agg$year==1990),"year"] <- 2003
acs_agg[which(acs_agg$year==2000),"year"] <- 2004

## recode barbers and hairdyers (although it does affect results at all)
occegp[which(occegp$occ1990==457),"egp"] <- "VII"
occegp[which(occegp$occ1990==458),"egp"] <- "VII"
acs_agg[which(acs_agg$occ1990==457),"egp"] <- "VII"
acs_agg[which(acs_agg$occ1990==458),"egp"] <- "VII"

## empty df to store cosine similarity in each year (DA typology)
cos_sim <-
  data.frame(year=c(2002:2021))

for (i in c(2002:2021)){

  ## create occupation vector
  occvec <-
    create_occvec(data=acs_agg[which(!is.na(acs_agg$da)),],
                  nu=50,
                  crosswalk=occda,
                  cps=FALSE)

  ## calculate centroid and weighted by individual counts
  occvec <- occvec_centroid(data = occvec)

  ## calculate cosine similarity difference
  cosine_similarity(occvec,category=8,standardize=T)

  ## monitor progress
  print(paste("year", as.character(i), "is done!"))
}

## plot
cos_sim_PM_base <-
  cos_sim[,c(1,4:9)]
cos_sim_PM_base <- cos_sim_PM_base %>%
  pivot_longer(
    cols = -year,
    names_to = "variable",
    values_to = "value"
  )
cos_sim_PM_base$variable <- as.character(cos_sim_PM_base$variable)
cos_sim_PM_base[cos_sim_PM_base$variable=="V4","variable"] <- "Sales Occupations"
cos_sim_PM_base[cos_sim_PM_base$variable=="V5","variable"] <- "Administrative and Office"
cos_sim_PM_base[cos_sim_PM_base$variable=="V6","variable"] <- "Production"
cos_sim_PM_base[cos_sim_PM_base$variable=="V7","variable"] <- "Laborers"
cos_sim_PM_base[cos_sim_PM_base$variable=="V8","variable"] <- "Clean and Protective"
cos_sim_PM_base[cos_sim_PM_base$variable=="V9","variable"] <- "Personal Services"

## ggplot
cos_sim_PM_base$variable <- factor(cos_sim_PM_base$variable, levels=c(
                                                                      "Personal Services",
                                                                      "Clean and Protective",
                                                                      "Administrative and Office",
                                                                      "Sales Occupations",
                                                                      "Laborers",
                                                                      "Production"
                                                                      ))

## shift series with min below -0.2 upward for visual clarity
shift_df <- cos_sim_PM_base %>%
  group_by(variable) %>%
  summarize(mn = min(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(shift = pmax(0, round(-mn + 0.1, 1))) %>%
  filter(shift > 0)
shift_info <- setNames(shift_df$shift, as.character(shift_df$variable))
if (length(shift_info) > 0) {
  cos_sim_PM_base <- cos_sim_PM_base %>%
    mutate(value = value + ifelse(as.character(variable) %in% names(shift_info),
                                  shift_info[as.character(variable)], 0))
}

lo   <- floor(min(cos_sim_PM_base$value, na.rm=TRUE) / 0.2) * 0.2 - 0.05
hi   <- ceiling(max(cos_sim_PM_base$value, na.rm=TRUE) / 0.2) * 0.2 + 0.2
step <- 0.2

ggplot(cos_sim_PM_base,aes(x=year,
                           y=value,
                           group=variable,
                           color=variable,
                           fill=variable,
                           shape=variable,
                           lty=variable)) +
  geom_point(aes(size=variable)) +
  geom_line() +
  scale_x_continuous(
    breaks = seq(2002, 2021, 3),
    labels = c("1980", "2005", "2008", "2011", "2014", "2017", "2020")
  ) +
  scale_y_continuous(limits = c(lo, hi), breaks = seq(lo, hi, step)) +
  theme_classic() +
  scale_colour_manual(values=c("#82243b","grey30","grey30",
                               "#82243b","#ed5278","grey30")) +
  scale_fill_manual(values=c("#82243b","#ed5278","#ffb6c1",
                             NA,NA,NA)) +
  scale_shape_manual(values=c(21,21,21,21,21,21)) +
  scale_linetype_manual(values=c(2,3,3,3,3,3)) +
  scale_size_manual(values=c(3.5,3.5,3.5,3.5,3.5,3.5)) +
  ylab("cosine similarity") +
  ggtitle("ACS and Census, 1980-2021") +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14,angle=90,hjust=0.5),
        axis.title.x = element_text(size=16,hjust=1,vjust=8, margin = margin(t = -12)),
        axis.title.y = element_text(size=16,angle=90),
        plot.title = element_text(size=18,hjust=0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=16),
        legend.position = "bottom") +
  guides(color=guide_legend(ncol=1),
         size = guide_legend(direction = "vertical"))

ggsave(
  file.path(fig_dir, "ACS_DA_PM_base.png"),
  width = 11.5, height = 16, units = "cm", dpi = 800
)
## tiff for submission
ggsave(
  file.path(fig_dir, "ACS_DA_PM_base.tiff"),
  device = ragg::agg_tiff,
  width = 11.5, height = 16, units = "cm", dpi = 800, compression = "lzw"
)


## ── Section 2: Location heterogeneity (Fig 4 upper) ──────────────────────────
## → Output: Analysis/figures/Cosine Similarity/ACS_DA_location_heterogeneity.png
##           (Main paper Fig 4)

## drop each CZ and run the analysis
cz <- unique(cz_density$metarea)
cos_sim_cz_dropped <- data.frame()

Nocc <-
  acs_agg %>% dplyr::group_by(year,occ1990) %>%
  dplyr::summarize(Nocc=sum(wtfinl,na.rm=T))

## define the function
## create occupation vectors
create_occvec <- function(data=cps,
                          nu=200,
                          crosswalk=occegp,
                          cps=TRUE){

  ## create Occupation-Location Matrix in each year

  ## TF
  olm <- with(five_rd(data,i), questionr::wtd.table(occ1990, metarea, wtfinl)) %>%
    prop.table(margin=1)
  olm <- scale(olm,center=TRUE,scale=TRUE) ## Occ-Loc Matrix is centered and scaled by default

  ## calculate number of people (weighted) in each occupation in each year for later centroid weight
  Nocc <-
    data %>% dplyr::group_by(year,occ1990) %>%
      dplyr::summarize(Nocc=sum(wtfinl,na.rm=T))

  ## SVD decomposition and a vector representation of occupations
  svd <- svd(olm,nu=nu)
  occvec <- svd$u %*% diag(svd$d)[1:nu,1:nu]

  ## create class identifier and add occupation weight
  occvec <- t(apply(occvec,1,l2.norm))
  occvec <- as.data.frame(occvec)

  occvec$occ1990 <- rownames(olm)

  occvec <- merge(occvec,crosswalk,by="occ1990",all.x=T)
  occvec <- merge(occvec,Nocc[which(Nocc$year %in% five_mw(i)),]
                  %>% group_by(occ1990) %>%
                    summarize(Nocc=sum(Nocc,na.rm=T)),
                  by="occ1990",all.x=T)
  occvec <- occvec[,-which(names(occvec) %in% c("occ1990"))]

  return(occvec)

}

## loop over CZs
for (l in 1:length(cz)){

  ## empty df to store cosine similarity in each year
  cos_sim <-
    data.frame(year=2002:2021)

  for (i in 2002:2021){

    ## create occupation vector
    occvec <-
      create_occvec(data=acs_agg[which(!is.na(acs_agg$da)),] %>% filter(metarea!=cz[l]),
                         nu=50,
                         crosswalk=occda,
                         cps = FALSE)

    ## calculate centroid and weighted by individual counts
    occvec <- occvec_centroid(data = occvec)

    ## calculate cosine similarity
    cosine_similarity(occvec,category=8,standardize=T)

  }

  ## save results
  cos_sim$cz <- cz[l]
  cos_sim_cz_dropped <- rbind(cos_sim_cz_dropped,cos_sim)

  ## monitor progress
  if(l %% 20 == 0){
    cat("Finished", l, "of", length(cz), "\n")
  }
}

## since the above process takes a long time, read the saved output
cos_sim_cz_dropped <- read.csv(file.path(data_dir, "output data/cos_sim_cz_dropped_acs.csv"))

## select PM-personal service pair
change_cz <- cos_sim_cz_dropped[, c(1, 9, 66)] %>%
  pivot_longer(
    cols = -c(year, cz),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(cz) %>%
  summarize(
    cos_change = value[year == 2021] - value[year == 2002],
    .groups = "drop"
  )
change_cz$cos_change <- change_cz$cos_change - (cos_sim[,c(1,9)] %>%
                                                      summarize(cos_change = V9[year==2021]-V9[year==2002]) %>% ## 2002 is the 1980 data point after re-scale
                                                      pull())

## merge density
change_cz <- merge(
  cz_density,
  change_cz,
  by.x="metarea",
  by.y="cz",
  all.y=T
)
change_cz$density <- as.numeric(change_cz$density)
change_cz <- merge(change_cz,
                   lma %>% group_by(CZ) %>%
                     arrange(desc(Labor_Force)) %>%
                     slice(1) %>%
                     dplyr::select(CZ, County_Name),
                   by.x="metarea",
                   by.y="CZ",
                   all.x=T)

## plot
change_cz <- change_cz %>%
  mutate(
    dominant = case_when(
      cos_change <= -0.022 ~ "dominant",
      TRUE                 ~ "others"
    ),
    dominant = factor(dominant, levels = c("dominant", "others"))
  )

make_panel <- function(data, xvar, logx = FALSE, title = NULL, xlim = NULL, xlabels = NULL, breaks = NULL) {
  data2 <- data %>% mutate(x = if (logx) log(.data[[xvar]]) else .data[[xvar]])

  p <- data2 %>%
    mutate(County_Name = gsub(" County", "", County_Name)) %>%
    ggplot(aes(x = x, y = cos_change)) +
    geom_point(aes(size = total_pop, shape = dominant, color = dominant, alpha = dominant)) +
    geom_smooth(method = "loess", se = FALSE, color = "#82243b") +
    geom_text_repel(
      data = ~ subset(.x, dominant == "dominant" & County_Name != "Kings, NY"),
      aes(label = County_Name),
      size = 4, color = "#82243b", box.padding = 0.35, point.padding = 0.5,
      segment.color = "#82243b", max.overlaps = Inf, min.segment.length = 0.2,
      force = 8, show.legend = FALSE
    ) +
    geom_text_repel(
      data = ~ subset(.x, County_Name == "Kings, NY"),
      aes(label = County_Name),
      nudge_y = 0.01, nudge_x = 0.0,
      size = 4, color = "#82243b", box.padding = 0.35, point.padding = 0.5,
      segment.color = "#82243b", max.overlaps = Inf, force = 8, show.legend = FALSE
    ) +
    scale_shape_manual(values = c(dominant = 5, others = 1)) +
    scale_color_manual(values = c(dominant = "#82243b", others = "grey30")) +
    scale_alpha_manual(values = c(dominant = 1, others = 0.5)) +
    scale_size_continuous(range = c(1, 3)) +
    scale_y_continuous(limits = c(-0.13, 0.04), breaks = seq(-0.13,0.04,0.04)) +
    ylab("cosine change") +
    xlab(title) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
      axis.title.x = element_text(size= 14, hjust=1,vjust=8, margin = margin(t = -10)),
      axis.title.y = element_text(size = 14, angle = 90),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )

  if (is.null(breaks)) {
    data_range <- range(data2$x, na.rm = TRUE)
    plot_range <- if (!is.null(xlim)) xlim else data_range
    brks <- seq(plot_range[1], plot_range[2], length.out = 4)
  } else {
    brks <- breaks
  }

  if (!is.null(xlabels)) {
    p <- p + scale_x_continuous(limits = xlim, breaks = brks, labels = xlabels)
  } else {
    p <- p + scale_x_continuous(limits = xlim, breaks = brks)
  }

  return(p)
}

panels <- list(
  list(var = "density",         logx = TRUE,  title = "labor force density (log)", xlabels = function(x) sprintf("%.2f", x)),
  list(var = "immigrant_share", logx = TRUE,  title = "immigrant share (log)", xlim = c(-5.5, -1), xlabels = function(x) sprintf("%.2f", x)),
  list(var = "incwage_cpi",     logx = TRUE,  title = "income (log)", xlabels = function(x) sprintf("%.2f", x)),
  list(var = "college",         logx = FALSE, title = "college share", xlabels = function(x) sprintf("%.2f", x))
)

plots <- lapply(panels, function(pinfo) {
  make_panel(
    change_cz,
    xvar   = pinfo$var,
    logx   = pinfo$logx,
    title  = pinfo$title,
    xlim   = pinfo$xlim %||% NULL,
    xlabels= pinfo$xlabels %||% NULL
  )
})

layout_plot <- (plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]])
print(layout_plot)

ggsave(
  file.path(fig_dir, "ACS_DA_location_heterogeneity.jpg"),
  width = 16.5*1.2, height = 14*1.2, units = "cm", dpi = 600
)
ggsave(
  file.path(fig_dir, "ACS_DA_location_heterogeneity.tiff"),
  device = ragg::agg_tiff,
  width = 16.5*1.2, height = 14*1.2, units = "cm", dpi = 2500, compression = "lzw"
)


## ── Section 3: Occupation heterogeneity + OLS (Fig 4 lower + Fig 5) ──────────
## → Output: Analysis/figures/Cosine Similarity/ACS_DA_occ_heterogeneity.png (Fig 4 lower)
##           Analysis/figures/Cosine Similarity/ACS_DA_occ_ols.png (Fig 5)

occvecs <- data.frame()

## slightly revise create_occvec to keep occ1990 label in the vector output
create_occvec <- function(data=cps,
                          nu=200,
                          crosswalk=occegp,
                          cps=TRUE){

  olm <- with(five_rd(data,i), questionr::wtd.table(occ1990, metarea, wtfinl)) %>%
    prop.table(margin=1)
  olm <- scale(olm,center=TRUE,scale=TRUE)

  Nocc <-
    data %>% dplyr::group_by(year,occ1990) %>%
      dplyr::summarize(Nocc=sum(wtfinl,na.rm=T))

  svd <- svd(olm,nu=nu)
  occvec <- svd$u %*% diag(svd$d)[1:nu,1:nu]

  occvec <- t(apply(occvec,1,l2.norm))
  occvec <- as.data.frame(occvec)

  occvec$occ1990 <- rownames(olm)

  occvec <- merge(occvec,crosswalk,by="occ1990",all.x=T)
  occvec <- merge(occvec,Nocc[which(Nocc$year %in% five_mw(i)),]
                  %>% group_by(occ1990) %>%
                    summarize(Nocc=sum(Nocc,na.rm=T)),
                  by="occ1990",all.x=T)

  return(occvec)
}

for (i in 2002:2021){
  occvec <-
    create_occvec(data=acs_agg[which(!is.na(acs_agg$da)),],
                  nu=50,
                  crosswalk=occda,
                  cps = FALSE)
  occvec$year <- i
  occvecs <- rbind(occvecs,occvec)
}

occ <- unique(
  acs_agg %>%
    distinct(occ1990) %>%
    arrange(occ1990) %>%
    pull(occ1990)
)

cos_sim_occ_dropped <- data.frame()

for (l in 1:length(occ)){

  cos_sim <-
    data.frame(year=2002:2021)

  for (i in 2002:2021){

    occvec <- occvecs[which(occvecs$occ1990!=occ[l]),] %>%
      dplyr::select(-"occ1990") %>%
      filter(year==i) %>%
      dplyr::select(-"year")
    occvec <- occvec_centroid(data = occvec)

    cosine_similarity(occvec,category=8,standardize=T)
    cos_sim$occ <- occ[l]
  }

  if(l %% 20 == 0){
    cat("Finished", l, "of", length(occ), "\n")}

  cos_sim_occ_dropped <- rbind(cos_sim_occ_dropped, cos_sim)
}

## with personal service
t <- merge(cos_sim_occ_dropped, occda %>% mutate(occ1990=as.character(occ1990)),
      by.x="occ",
      by.y="occ1990",
      all.x=T) %>%
  mutate(da=da-1,
         da=case_when(da==0~1,
                      .default = da)) %>%
  filter(da!="8")

change_occs <- data.frame()
for (g in seq(1,7)) {
  change_occ <- melt(t[t$da==g,c(1,2,1+1+8*g)], id.vars=c("year","occ")) %>%
    group_by(occ) %>%
    summarize(cos_change = value[year==2021]-value[year==2002])
  change_occ$cos_change <- change_occ$cos_change - (cos_sim[,c(1,1+8*g)] %>%
                                                      summarize(cos_change = .data[[paste0("V", 1 + 8*g)]][year==2021]-
                                                                  .data[[paste0("V", 1 + 8*g)]][year==2002]) %>%
                                                      pull())
  change_occs <- rbind(change_occs,change_occ)
}
change_occ_pm <- change_occs

## with professional and managerial
t <- merge(cos_sim_occ_dropped, occda %>% mutate(occ1990=as.character(occ1990)),
           by.x="occ",
           by.y="occ1990",
           all.x=T) %>%
  mutate(da=da-1,
         da=case_when(da==0~1,
                      .default = da)) %>%
  filter(da!="1")

change_occs <- data.frame()
for (g in seq(2,8)) {
  change_occ <- melt(t[t$da==g,c(1,2,1+1+1+8*(g-1))], id.vars=c("year","occ")) %>%
    group_by(occ) %>%
    summarize(cos_change = value[year==2021]-value[year==2002])
  change_occ$cos_change <- change_occ$cos_change - (cos_sim[,c(1,2+8*(g-1))] %>%
                                                      summarize(cos_change = .data[[paste0("V", 2 + 8*(g-1))]][year==2021]-
                                                                  .data[[paste0("V", 2 + 8*(g-1))]][year==2002]) %>%
                                                      pull())
  change_occs <- rbind(change_occs,change_occ)
}
change_occ_service <- change_occs
change_occ <- rbind(
  merge(change_occ_pm,occda %>% mutate(da=as.character(da)),
        by.x="occ",
        by.y="occ1990",
        all.x=T) %>% filter(da=="1"),
  merge(change_occ_service,occda %>% mutate(da=as.character(da)),
        by.x="occ",
        by.y="occ1990",
        all.x=T) %>% filter(da=="9")
)

## merge covariates
change_occ <-
  merge(change_occ,
        covariate,
        by.x="occ",
        by.y="occ1990",
        all.x=T)
change_occ[which(change_occ$occ1990==349),"da"] <- 5

change_occ <-
  merge(change_occ,
        read_dta(file.path(data_dir, "misc_data/occ1990_titles.dta")) %>%
          mutate(occ1990=as.character(occ1990)),
        by.x="occ",
        by.y="occ1990",
        all.x=T)

## plot top PM and service occupations
change_occ %>%
  filter(da %in% c(1, 9)) %>%
  mutate(da = as.character(da)) %>%
  filter(!grepl("Health", title) & !grepl("Nursing", title)) %>%
  group_by(da) %>%
  arrange(cos_change) %>%
  slice_head(n = 9) %>%
  ungroup() %>%
  mutate(title = case_when(
    grepl("Kitchen", title)       ~ "Kitchen worker",
    grepl("admin", title)       ~ "Managers and administrators",
    grepl("Health aides", title)  ~ "Health aides",
    grepl("Nursing aides", title) ~ "Nursing aides",
    grepl("Personal,", title) ~ "Labor relation specialists",
    grepl("Personal service", title) ~ "Personal service n.e.c",
    TRUE                          ~ title
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
                    labels = c("1" = "Professional-Management",
                               "9" = "Personal Service"),
                    name = NULL) +
  facet_grid(da ~ ., scales = "free_y", space = "free") +
  theme_classic() +
  ggtitle("Top Contributing Occupations, Census and ACS, 1980-2021") +
  labs(x = "", y = "cosine change") +
  theme(
    axis.text.x     = element_text(size = 12),
    axis.text.y     = element_text(size = 12, hjust = 1),
    axis.title      = element_text(size = 13),
    legend.title    = element_text(size = 13),
    legend.text     = element_text(size = 12),
    legend.position = "bottom",
    plot.title = element_text(hjust=1.2),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, keyheight = unit(0.4, "lines")))

ggsave(
  file.path(fig_dir, "ACS_DA_occ_heterogeneity.png"),
  width = 14, height = 13, units = "cm", dpi = 800
)

## O*NET occupational characteristics for personal service occupations (da == 9).
## All values normalized to 0-100 scale. Six variables:
##   physical_prox      : Physical Proximity (Work Context), O*NET 23.1
##   face2face          : Face-to-Face Discussions (Work Context), O*NET 23.1
##   interpersonal      : Social Orientation (Work Styles), O*NET 23.1
##   ext_customers      : Deal With External Customers (Work Context), O*NET 23.1
##   responsible_health : Assisting and Caring for Others (Work Activities), O*NET 13.0
##                        (used instead of "Responsible for Others' Health and Safety"
##                         which gives implausibly low values for health/care workers)
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
change_occ <- merge(change_occ, onet_2023, by.x = "occ", by.y = "occ", all.x = TRUE)

## save bivariate regression results
change_occ <- change_occ %>%
  group_by(da) %>%
  mutate(
    cos_change_z          = scale(cos_change)[,1],
    overwork_share_z          = scale(overwork_share)[,1],
    overwork_male_z          = scale(overwork_male)[,1],
    overwork_female_z          = scale(overwork_female)[,1],
    fulltime_z          = scale(1-parttime_share)[,1],
    log_incwage_z       = scale(log(incwage))[,1],
    share_college_z     = scale(share_college)[,1],
    uhrswork_male_z     = scale(uhrswork_male)[,1],
    uhrswork_female_z   = scale(uhrswork_female)[,1],
    share_female_z      = scale(share_female)[,1],
    share_immigration_z    = scale(share_immigration)[,1],
    physical_prox_z     = scale(physical_prox)[,1],
    face2face_z         = scale(face2face)[,1],
    interpersonal_z     = scale(interpersonal)[,1],
    ext_customers_z     = scale(ext_customers)[,1],
    responsible_z       = scale(responsible_health)[,1],
    share_white_z       = scale(share_white)[,1],
    share_black_z       = scale(share_black)[,1],
    share_asian_z       = scale(share_asian)[,1],
    share_hispanic_z    = scale(share_hispanic)[,1]
  ) %>%
  ungroup()

results <- data.frame(
  group=c(rep("PM",10),rep("Service",10)),
  variable=c(c("Annual income (log)", "Share of overwork (male)",
               "Share of overwork (female)",
             "    Weekly working hours (male)", "    Weekly working hours (female)",
             "Share of full-time working",
             "Share of white", "Share of Black", "Share of Asian",
             "Share of Hispanic"),
  c("Share of white", "Share of Black", "Share of Asian",
    "Share of Hispanic", "Share of immigrants",
    "Physical proximity",
    "Face-to-face discussions", "Interpersonally oriented",
    "Serves external clients", "Care responsibility")),
  feature=c(c("demand","demand","demand","demand","demand","demand",
              "race/ethnicity","race/ethnicity",
              "race/ethnicity","race/ethnicity"),
            c("race/ethnicity","race/ethnicity",
            "race/ethnicity","race/ethnicity","",
            "work content/style","work content/style",
            "work content/style","work content/style",
            "work content/style")),
                      coef=NA,
                      se=NA)

for (i in 1:dim(results[results$group=="PM",])[1]) {
  v <- c("log_incwage_z","overwork_male_z","overwork_female_z",
         "uhrswork_male_z","uhrswork_female_z","fulltime_z",
         "share_white_z","share_black_z","share_asian_z","share_hispanic_z")[i]
  fit <- lm(as.formula(paste("cos_change_z", "~", v)), data = change_occ %>% filter(da == 1))
  s <- summary(fit)
  results$coef[i] <- s$coefficients[2, "Estimate"]
  results$se[i]   <- s$coefficients[2, "Std. Error"]
  rm(v,s)
}

for (i in 1:dim(results[results$group=="Service",])[1]) {
  v <- c("share_white_z","share_black_z","share_asian_z","share_hispanic_z",
         "share_immigration_z","physical_prox_z","face2face_z",
         "interpersonal_z","ext_customers_z","responsible_z")[i]
  fit <- lm(as.formula(paste("cos_change_z", "~", v)), data = change_occ %>% filter(da == 9), weight = weight)
  s <- summary(fit)
  results$coef[i+10] <- s$coefficients[2, "Estimate"]
  results$se[i+10]   <- s$coefficients[2, "Std. Error"]
  rm(v,s)
}

results %>%
  mutate(lo = coef - 1.96 * se, hi = coef + 1.96 * se) %>%
  {
    pm_df <- filter(., group == "PM") %>% mutate(variable = as.character(variable))
    svc_df <- filter(., group == "Service") %>% mutate(variable = as.character(variable))

    x_m_pm  <- max(abs(c(pm_df$lo, pm_df$hi)), na.rm = TRUE)
    x_m_svc <- max(abs(c(svc_df$lo, svc_df$hi)), na.rm = TRUE)
    xlim_pm  <- c(-x_m_pm-0.2,  x_m_pm+0.2)
    xlim_svc <- c(-x_m_svc-1, x_m_svc+1)

    reorder_within_feature <- function(df) {
      df %>%
        group_by(feature) %>%
        arrange(desc(coef)) %>%
        mutate(.order_in_feat = row_number()) %>%
        ungroup() %>%
        group_by(feature) %>%
        mutate(feature_mean = mean(coef, na.rm = TRUE)) %>%
        ungroup() %>%
        arrange(desc(feature_mean), feature, .order_in_feat) %>%
        mutate(variable = factor(variable, levels = unique(variable)))
    }

    pm_df  <- reorder_within_feature(pm_df)
    svc_df <- reorder_within_feature(svc_df)

    sep_pm <- pm_df %>%
      mutate(.y = as.numeric(variable)) %>%
      group_by(feature) %>%
      summarize(max_y = max(.y), .groups = "drop") %>%
      arrange(max_y) %>%
      pull(max_y) %>%
      { if(length(.)>1) head(., -1) + 0.5 else numeric(0) }

    sep_svc <- svc_df %>%
      mutate(.y = as.numeric(variable)) %>%
      group_by(feature) %>%
      summarize(max_y = max(.y), .groups = "drop") %>%
      arrange(max_y) %>%
      pull(max_y) %>%
      { if(length(.)>1) head(., -1) + 0.5 else numeric(0) }

    p_pm <- ggplot(pm_df, aes(x = coef, y = variable)) +
      { if(length(sep_pm)>0) geom_hline(yintercept = sep_pm, color = "grey85", size = 0.35, inherit.aes = FALSE) } +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, color = "grey40") +
      geom_point(size = 3.7, shape = 21, fill = "#82243b", color = "grey40") +
      scale_y_discrete(limits = levels(pm_df$variable)) +
      coord_cartesian(xlim = xlim_pm, expand = T) +
      theme_classic() +
      labs(x = NULL, y = NULL, title = "Cosine Change (1980-2021) and Occupational Features") +
      theme(
        axis.text.x     = element_text(size = 12),
        axis.text.y     = element_text(size = 12, hjust = 1),
        axis.title      = element_text(size = 13),
        legend.title    = element_text(size = 13),
        legend.text     = element_text(size = 12),
        plot.title = element_text(hjust=1.2),
        legend.position = "bottom"
      )

    p_svc <- ggplot(svc_df, aes(x = coef, y = variable)) +
      { if(length(sep_svc)>0) geom_hline(yintercept = sep_svc, color = "grey85", size = 0.35, inherit.aes = FALSE) } +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, color = "grey40") +
      geom_point(size = 3.7, shape = 21, fill = "#ffe4ec", color = "grey40") +
      scale_y_discrete(limits = levels(svc_df$variable)) +
      scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
      coord_cartesian(xlim = xlim_svc, expand = T) +
      theme_classic() +
      geom_point(
        data = data.frame(group = c("Professional-Management", "Personal Service"),
                          coef = NA, variable = NA),
        aes(fill = group), shape = 21, size = 3.7, color = "black", show.legend = TRUE
      ) +
      scale_fill_manual(
        values = c("Personal Service" = "#82243b", "Professional-Management" = "#ffe4ec"),
        labels = c("Professional-Management","Personal Service"),
        name = NULL
      ) +
      labs(x = "coefficient (bivariate)", y = NULL, title = NULL) +
      theme(
        axis.text.x     = element_text(size = 12),
        axis.text.y     = element_text(size = 12, hjust = 1),
        axis.title      = element_text(size = 13),
        legend.title    = element_text(size = 13),
        legend.text     = element_text(size = 12),
        legend.position = "bottom"
      )

    heights_vec <- c(nrow(pm_df), nrow(svc_df))
    p_pm / p_svc + plot_layout(heights = heights_vec) +
      theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.key.height = unit(0.5, "lines")) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE, keyheight = unit(0.4, "lines"))) +
      theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0.5, unit = "cm"))
  }

ggsave(
  file.path(fig_dir, "ACS_DA_occ_ols.png"),
  width = 14, height = 17.5, units = "cm", dpi = 800
)


## ── Section 4: ACS EGP typology (SI S2) ──────────────────────────────────────
## → Output: Analysis/figures/Cosine Similarity/ACS_EGP_PM_base.png
##           (SI Appendix S2, ACS EGP figure)
##
## This section produces the standalone ACS EGP panel.

source(file.path(code_dir, "02 functions.R"))
select <- dplyr::select; filter <- dplyr::filter

## collapse EGP sub-classes in acs_agg (reuses data from Section 1)
acs_agg[acs_agg$egp %in% c("IIIa","IIIb"), "egp"] <- "III"
acs_agg[acs_agg$egp %in% c("VIIa"),         "egp"] <- "VII"
egp_keep <- c("I","II","III","V","VI","VII")

cos_sim <- data.frame(year = 2002:2021)

for (i in 2002:2021){

  occvec <-
    create_occvec(data = acs_agg[which(acs_agg$egp %in% egp_keep), ],
                  nu = 50, crosswalk = occegp, cps = FALSE)
  
  ## calculate centroid and weighted by individual counts
  occvec <- occvec_centroid(data = occvec)
  
  ## calculate cosine similarity
  cosine_similarity(occvec,category=6,standardize=T)
  
  ## monitor progress
  print(paste("year", as.character(i), "is done!"))
  
}

## V3=I_II, V4=I_III, V5=I_V, V6=I_VI, V7=I_VII (category=6, cols 3:7)
egp_level_order <- c("Class VII (Semi+Unskilled)", "Class II (Lower Service)",
                     "Class III (Routine Non-Manual)", "Class VI (Skilled Manual)",
                     "Class V (Manual Supervisors)")
acs_egp_df <- cos_sim[, c(1, 3:7)] %>%
  pivot_longer(-year, names_to="pair", values_to="value") %>%
  mutate(label = recode(pair,
    V3="Class II (Lower Service)",  V4="Class III (Routine Non-Manual)",
    V5="Class V (Manual Supervisors)", V6="Class VI (Skilled Manual)",
    V7="Class VII (Semi+Unskilled)")) %>%
  mutate(label = factor(label, levels=egp_level_order))

acs_x_breaks <- c(2002, 2006, 2010, 2014, 2018, 2021)
acs_x_labels <- c("1980","2006","2010","2014","2018","2021")

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

lvls <- levels(acs_egp_df$label)
p_acs_egp <- ggplot(acs_egp_df,
       aes(x=year, y=value, group=label,
           color=label, fill=label, shape=label, lty=label)) +
  geom_point(size=3.5) +
  geom_line() +
  scale_color_manual(values=egp_colors) +
  scale_fill_manual(values=egp_fills, na.value=NA) +
  scale_shape_manual(values=setNames(rep(21L, length(lvls)), lvls)) +
  scale_linetype_manual(values=egp_ltypes) +
  scale_x_continuous(breaks=acs_x_breaks, labels=acs_x_labels) +
  scale_y_continuous(breaks=seq(-2.5, 1.5, 0.5)) +
  ylab("cosine similarity (z-standardized)") +
  ggtitle("ACS/Census, EGP") +
  theme_classic() +
  theme(
    axis.text.x     = element_text(size=14),
    axis.text.y     = element_text(size=14, angle=90, hjust=0.5),
    axis.title.x    = element_text(size=16, hjust=1, vjust=8, margin=margin(t=-12)),
    axis.title.y    = element_text(size=16, angle=90),
    plot.title      = element_text(size=16, hjust=0.5),
    legend.title    = element_blank(),
    legend.text     = element_text(size=16),
    legend.position = "bottom",
    legend.direction= "vertical",
    legend.key.width= unit(1.2, "cm")
  ) +
  guides(color=guide_legend(ncol=1))

## fixed-height save: 10 cm plot + 6 cm legend grob
H_main <- 10; H_leg <- 6
get_legend_grob <- function(p) {
  gt  <- ggplot_gtable(ggplot_build(p))
  idx <- which(sapply(gt$grobs, function(x) x$name) == "guide-box")
  gt$grobs[[idx]]
}
p_clean <- p_acs_egp + theme(legend.position="none")
leg     <- get_legend_grob(p_acs_egp)
ggsave(file.path(fig_dir, "ACS_EGP_PM_base.png"),
       plot=arrangeGrob(p_clean, leg, nrow=2, heights=unit(c(H_main, H_leg),"cm")),
       width=11.5, height=H_main+H_leg, units="cm", dpi=800)
cat("ACS EGP figure saved.\n")


## ── Section 5: ACS tract analysis (SI S6) ────────────────────────────────────
## → Output: Analysis/figures/Cosine Similarity/acs_M_base_3_tract_pooled_alltracts.png
##           (SI Appendix S6)
##
## Part A downloads tract-level occupation data from the Census API and saves
## acs_tract.RData. Part B loads that file and produces the paper figure.
## Requires a Census API key: census_api_key("YOUR_KEY_HERE", install = TRUE)

library(tidycensus)
library(censusapi)
library(purrr)
library(Matrix)
library(irlba)
census_api_key("8aad65112dfa6a24854bc860e265658ae6c11af9", install = TRUE)

## ── Part A: Download tract data and save ─────────────────────────────────────

states <- unique(fips_codes$state)
states <- states[!is.na(states)][1:51]
pooled <- data.frame()

## ---- 2000 Census ----
vars_male <- c(
  management             = "Total!!Male!!Management, professional, and related occupations!!Management, business, and financial operations occupations!!Management occupations, except farmers and farm managers",
  business               = "Total!!Male!!Management, professional, and related occupations!!Management, business, and financial operations occupations!!Business and financial operations occupations",
  computer               = "Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Computer and mathematical occupations",
  architecture           = "Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Architecture and engineering occupations",
  science                = "Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Life, physical, and social science occupations",
  community_service      = "Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Community and social services occupations",
  legal                  = "Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Legal occupations",
  education              = "Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Education, training, and library occupations",
  arts                   = "Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Arts, design, entertainment, sports, and media occupations",
  healthcare_practitioner= "Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Healthcare practitioners and technical occupations",
  healthcare_support     = "Total!!Male!!Service occupations!!Healthcare support occupations",
  protective_service     = "Total!!Male!!Service occupations!!Protective service occupations",
  food_preparation       = "Total!!Male!!Service occupations!!Food preparation and serving related occupations",
  building_maintenance   = "Total!!Male!!Service occupations!!Building and grounds cleaning and maintenance occupations",
  personal_care          = "Total!!Male!!Service occupations!!Personal care and service occupations",
  sales_related          = "Total!!Male!!Sales and office occupations!!Sales and related occupations",
  office_admin           = "Total!!Male!!Sales and office occupations!!Office and administrative support occupations",
  farming_fishing        = "Total!!Male!!Farming, fishing, and forestry occupations",
  construction           = "Total!!Male!!Construction, extraction, and maintenance occupations!!Construction and extraction occupations",
  installation           = "Total!!Male!!Construction, extraction, and maintenance occupations!!Installation, maintenance, and repair occupations",
  production             = "Total!!Male!!Production, transportation, and material moving occupations!!Production occupations",
  transportation_material= "Total!!Male!!Production, transportation, and material moving occupations!!Transportation and material moving occupations"
)
vars_female <- gsub("!!Male!!", "!!Female!!", vars_male)

census_vars <- listCensusMetadata(name="dec/sf3", vintage=2000, type="variables") %>%
  filter(str_detect(name,"P05")) %>%
  filter(str_detect(concept,"SEX BY OCCUPATION FOR THE EMPLOYED CIVILIAN POPULATION 16 YEARS AND OVER")) %>%
  filter(str_detect(label,"Total!!")) %>%
  filter(label %in% c(vars_male, vars_female)) %>%
  pull(name)
census_vars <- setNames(census_vars,
  c(paste0(names(vars_male),"_male"), paste0(names(vars_female),"_female")))

census <- map_dfr(
  setdiff(sprintf("%02d", c(1:56)), c("03","07","14","43","52")),
  function(st) getCensus(name="dec/sf3", vintage=2000, region="tract:*",
                         vars=census_vars, regionin=paste0("state:",st)),
  .id = "state_fips"
) %>%
  pivot_longer(cols=starts_with("P"), names_to="variable", values_to="estimate") %>%
  mutate(variable = names(census_vars)[match(variable, census_vars)]) %>%
  mutate(state2  = if_else(str_length(state)==2, state,
                           str_pad(state_fips, 2, "left", "0")),
         county3 = str_pad(str_replace_all(county,"\\D",""), 3, "left", "0"),
         tract6  = str_pad(str_replace_all(str_replace_all(tract,"\\D",""),"",""), 6, "right", "0"),
         GEOID   = paste0(state2, county3, tract6)) %>%
  filter(str_length(GEOID)==11) %>%
  select(GEOID, variable, estimate) %>%
  tidyr::separate(variable, into=c("variable","sex"), sep="_(?=[^_]+$)") %>%
  mutate(year=2000, moe=NA)
pooled <- rbind(pooled, census)

## ---- 2005–2009 ACS (5-year ending 2009, middle year = 2007) ----
vars_male <- c(
  management             = "Estimate!!Total!!Male!!Management, professional, and related occupations!!Management, business, and financial occupations!!Management occupations",
  business               = "Estimate!!Total!!Male!!Management, professional, and related occupations!!Management, business, and financial occupations!!Business and financial operations occupations",
  computer               = "Estimate!!Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Computer and mathematical occupations",
  architecture           = "Estimate!!Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Architecture and engineering occupations",
  science                = "Estimate!!Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Life, physical, and social science occupations",
  community_service      = "Estimate!!Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Community and social services occupations",
  legal                  = "Estimate!!Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Legal occupations",
  education              = "Estimate!!Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Education, training, and library occupations",
  arts                   = "Estimate!!Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Arts, design, entertainment, sports, and media occupations",
  healthcare_practitioner= "Estimate!!Total!!Male!!Management, professional, and related occupations!!Professional and related occupations!!Healthcare practitioner and technical occupations",
  healthcare_support     = "Estimate!!Total!!Male!!Service occupations!!Healthcare support occupations",
  protective_service     = "Estimate!!Total!!Male!!Service occupations!!Protective service occupations",
  food_preparation       = "Estimate!!Total!!Male!!Service occupations!!Food preparation and serving related occupations",
  building_maintenance   = "Estimate!!Total!!Male!!Service occupations!!Building and grounds cleaning and maintenance occupations",
  personal_care          = "Estimate!!Total!!Male!!Service occupations!!Personal care and service occupations",
  sales_related          = "Estimate!!Total!!Male!!Sales and office occupations!!Sales and related occupations",
  office_admin           = "Estimate!!Total!!Male!!Sales and office occupations!!Office and administrative support occupations",
  farming_fishing        = "Estimate!!Total!!Male!!Farming, fishing, and forestry occupations",
  construction           = "Estimate!!Total!!Male!!Construction, extraction, maintenance, and repair occupations!!Construction and extraction occupations",
  installation           = "Estimate!!Total!!Male!!Construction, extraction, maintenance, and repair occupations!!Installation, maintenance, and repair occupations",
  production             = "Estimate!!Total!!Male!!Production, transportation, and material moving occupations!!Production occupations",
  transportation_material= "Estimate!!Total!!Male!!Production, transportation, and material moving occupations!!Transportation and material moving occupations"
)
vars_female <- gsub("!!Male!!", "!!Female!!", vars_male)
acs_vars <- load_variables(2009, "acs5", cache=TRUE) %>%
  filter(str_detect(label,"occupation"), str_detect(label,"Total!!"),
         label %in% c(vars_male, vars_female),
         str_detect(name,"C"), !str_detect(concept,"FULL-TIME")) %>%
  pull(name)
acs_vars <- setNames(acs_vars,
  c(paste0(names(vars_male),"_male"), paste0(names(vars_female),"_female")))

acs_09 <- map_dfr(states,
  ~get_acs(geography="tract", variables=acs_vars, survey="acs5",
           state=.x, year=2009, cache_table=TRUE), .id="state_abbr") %>%
  tidyr::separate(variable, into=c("variable","sex"), sep="_(?=[^_]+$)") %>%
  select(GEOID, variable, sex, estimate, moe) %>%
  mutate(year=2007)
pooled <- rbind(pooled, acs_09)

## ---- 2010–2017 ACS ----
vars_male <- c(
  management             = "Estimate!!Total!!Male!!Management, business, science, and arts occupations!!Management, business, and financial occupations!!Management occupations",
  business               = "Estimate!!Total!!Male!!Management, business, science, and arts occupations!!Management, business, and financial occupations!!Business and financial operations occupations",
  computer               = "Estimate!!Total!!Male!!Management, business, science, and arts occupations!!Computer, engineering, and science occupations!!Computer and mathematical occupations",
  architecture           = "Estimate!!Total!!Male!!Management, business, science, and arts occupations!!Computer, engineering, and science occupations!!Architecture and engineering occupations",
  science                = "Estimate!!Total!!Male!!Management, business, science, and arts occupations!!Computer, engineering, and science occupations!!Life, physical, and social science occupations",
  community_service      = "Estimate!!Total!!Male!!Management, business, science, and arts occupations!!Education, legal, community service, arts, and media occupations!!Community and social service occupations",
  legal                  = "Estimate!!Total!!Male!!Management, business, science, and arts occupations!!Education, legal, community service, arts, and media occupations!!Legal occupations",
  education              = "Estimate!!Total!!Male!!Management, business, science, and arts occupations!!Education, legal, community service, arts, and media occupations!!Education, training, and library occupations",
  arts                   = "Estimate!!Total!!Male!!Management, business, science, and arts occupations!!Education, legal, community service, arts, and media occupations!!Arts, design, entertainment, sports, and media occupations",
  healthcare_practitioner= "Estimate!!Total!!Male!!Management, business, science, and arts occupations!!Healthcare practitioners and technical occupations",
  healthcare_support     = "Estimate!!Total!!Male!!Service occupations!!Healthcare support occupations",
  protective_service     = "Estimate!!Total!!Male!!Service occupations!!Protective service occupations",
  food_preparation       = "Estimate!!Total!!Male!!Service occupations!!Food preparation and serving related occupations",
  building_maintenance   = "Estimate!!Total!!Male!!Service occupations!!Building and grounds cleaning and maintenance occupations",
  personal_care          = "Estimate!!Total!!Male!!Service occupations!!Personal care and service occupations",
  sales_related          = "Estimate!!Total!!Male!!Sales and office occupations!!Sales and related occupations",
  office_admin           = "Estimate!!Total!!Male!!Sales and office occupations!!Office and administrative support occupations",
  farming_fishing        = "Estimate!!Total!!Male!!Natural resources, construction, and maintenance occupations!!Farming, fishing, and forestry occupations",
  construction           = "Estimate!!Total!!Male!!Natural resources, construction, and maintenance occupations!!Construction and extraction occupations",
  installation           = "Estimate!!Total!!Male!!Natural resources, construction, and maintenance occupations!!Installation, maintenance, and repair occupations",
  production             = "Estimate!!Total!!Male!!Production, transportation, and material moving occupations!!Production occupations",
  transportation         = "Estimate!!Total!!Male!!Production, transportation, and material moving occupations!!Transportation occupations",
  material_moving        = "Estimate!!Total!!Male!!Production, transportation, and material moving occupations!!Material moving occupations"
)
vars_female <- gsub("!!Male!!", "!!Female!!", vars_male)
acs_vars <- load_variables(2017, "acs5", cache=TRUE) %>%
  filter(str_detect(label,"occupation"), str_detect(label,"Total!!"),
         label %in% c(vars_male, vars_female),
         str_detect(name,"C"), !str_detect(concept,"FULL-TIME")) %>%
  pull(name)
acs_vars <- setNames(acs_vars,
  c(paste0(names(vars_male),"_male"), paste0(names(vars_female),"_female")))

for (i in seq(2010, 2017)) {
  tmp <- map_dfr(states,
    ~get_acs(geography="tract", variables=acs_vars, survey="acs5",
             state=.x, year=i, cache_table=TRUE), .id="state_abbr") %>%
    tidyr::separate(variable, into=c("variable","sex"), sep="_(?=[^_]+$)") %>%
    select(GEOID, variable, sex, estimate, moe) %>%
    mutate(year = i - 2)
  pooled <- rbind(pooled, tmp)
}

## ---- 2018 ACS ----
vars_male["education"] <- "Estimate!!Total!!Male!!Management, business, science, and arts occupations!!Education, legal, community service, arts, and media occupations!!Educational instruction, and library occupations"
vars_female <- gsub("!!Male!!", "!!Female!!", vars_male)
acs_vars <- load_variables(2018, "acs5", cache=TRUE) %>%
  filter(str_detect(label,"occupation"),
         label %in% c(vars_male, vars_female),
         str_detect(name,"C"), !str_detect(concept,"FULL-TIME")) %>%
  pull(name)
acs_vars <- setNames(acs_vars,
  c(paste0(names(vars_male),"_male"), paste0(names(vars_female),"_female")))

tmp <- map_dfr(states,
  ~get_acs(geography="tract", variables=acs_vars, survey="acs5",
           state=.x, year=2018, cache_table=TRUE), .id="state_abbr") %>%
  tidyr::separate(variable, into=c("variable","sex"), sep="_(?=[^_]+$)") %>%
  select(GEOID, variable, sex, estimate, moe) %>%
  mutate(year=2016)
pooled <- rbind(pooled, tmp)

## ---- 2019–2023 ACS ----
add_colons_after_first <- function(x) {
  prefix <- sub("!!.*$","!!",x)
  rest   <- sub("^[^!]*!!","",x)
  paste0(prefix, gsub("!!",":!!",rest,fixed=TRUE))
}
vars_male <- vapply(vars_male, add_colons_after_first, character(1))
vars_male["protective_service"]     <- "Estimate!!Total:!!Male:!!Service occupations:!!Protective service occupations:"
vars_male["healthcare_practitioner"]<- "Estimate!!Total:!!Male:!!Management, business, science, and arts occupations:!!Healthcare practitioners and technical occupations:"
vars_female <- gsub("!!Male","!!Female",vars_male)
acs_vars <- load_variables(2022,"acs5",cache=TRUE) %>%
  filter(str_detect(label,"occupation"),
         label %in% c(vars_male,vars_female),
         str_detect(name,"C"),
         !str_detect(concept,"FULL-TIME"), !str_detect(concept,"Full-Time")) %>%
  pull(name)
acs_vars <- setNames(acs_vars,
  c(paste0(names(vars_male),"_male"), paste0(names(vars_female),"_female")))

for (i in seq(2019, 2023)) {
  tmp <- map_dfr(states,
    ~get_acs(geography="tract", variables=acs_vars, survey="acs5",
             state=.x, year=i, cache_table=TRUE), .id="state_abbr") %>%
    tidyr::separate(variable, into=c("variable","sex"), sep="_(?=[^_]+$)") %>%
    select(GEOID, variable, sex, estimate, moe) %>%
    mutate(year = i - 2)
  pooled <- rbind(pooled, tmp)
}

## combine transportation + material_moving into one category
pooled <- pooled %>%
  bind_rows(
    pooled %>%
      filter(variable %in% c("transportation","material_moving")) %>%
      group_by(GEOID, sex, year) %>%
      summarize(estimate=sum(estimate,na.rm=TRUE),
                moe=if(any(is.na(moe))) NA_real_ else sqrt(sum(moe^2,na.rm=TRUE)),
                .groups="drop") %>%
      mutate(variable="transportation_material")
  ) %>%
  filter(!variable %in% c("transportation","material_moving")) %>%
  arrange(GEOID, year, variable, sex)

rm(tmp, acs_09, census, acs_vars, census_vars, vars_female, vars_male,
   add_colons_after_first, i)
save(pooled, file=file.path(data_dir, "acs_tract.RData"))
cat("Tract data saved to acs_tract.RData\n")

## ── Part B: Pooled cosine similarity figure (paper figure) ───────────────────

if (!exists("pooled")) load(file.path(data_dir, "acs_tract.RData"))

## recode occupation labels to SOC codes
occ_recode <- c(
  management="11", business="13", computer="15", architecture="17",
  science="19", community_service="21", legal="23", education="25",
  arts="27", healthcare_practitioner="29", healthcare_support="31",
  protective_service="33", food_preparation="35", building_maintenance="37",
  personal_care="39", sales_related="41", office_admin="43",
  farming_fishing="45", construction="47", installation="49",
  production="51", transportation_material="53"
)
pooled <- pooled %>%
  mutate(variable = ifelse(variable %in% names(occ_recode),
                           occ_recode[variable], variable))

## pool both sexes; restrict to 2007+
pooled_both <- pooled %>%
  filter(year >= 2007) %>%
  group_by(GEOID, year, variable) %>%
  summarize(estimate=sum(estimate,na.rm=TRUE), .groups="drop")

cos_sim <- data.frame(year=2007:2021)
library(irlba)
for (i in 2007:2021) {
  dat <- pooled_both %>% filter(year==i)
  olm <- with(dat, questionr::wtd.table(variable, GEOID, estimate)) %>%
    prop.table(margin=2)
  olm <- scale(olm, center=TRUE, scale=TRUE); olm[is.na(olm)] <- 0
  sv  <- irlba(Matrix(olm, sparse=TRUE), nu=10, nv=10)
  occvec <- sv$u %*% diag(sv$d[1:10])
  cosine_similarity(occvec, category=22, standardize=TRUE)
  print(paste("year", i, "done"))
}

## management (col 1) vs healthcare support, protective service, building,
## personal care, sales, office admin — positions 12,13,15,16,17,18
compare <- c(12, 13, 15, 16, 17, 18)
cos_sim_tract <- cos_sim[, c(1, compare)] %>%
  pivot_longer(-year, names_to="variable", values_to="value") %>%
  mutate(variable = recode(variable,
    V12="Healthcare Support", V13="Protective Service",
    V15="Building Maintenance", V16="Personal Care and Service",
    V17="Sales and Related", V18="Office and Admin Support")) %>%
  mutate(variable = factor(variable, levels=c(
    "Personal Care and Service","Protective Service","Building Maintenance",
    "Healthcare Support","Sales and Related","Office and Admin Support")))

ggplot(cos_sim_tract,
       aes(x=year, y=value, group=variable, color=variable, fill=variable,
           shape=variable, lty=variable)) +
  geom_point(aes(size=variable)) +
  geom_line() +
  scale_x_continuous(breaks=seq(2007,2021,2)) +
  scale_y_continuous(limits=c(-2.2,2.2), breaks=seq(-2.2,2.2,0.7)) +
  theme_classic() +
  scale_colour_manual(values=c("#82243b","grey30","grey30","#82243b","#ed5278","grey30")) +
  scale_fill_manual(values=c("#82243b","#ed5278","#ffb6c1",NA,NA,NA)) +
  scale_shape_manual(values=rep(21,6)) +
  scale_linetype_manual(values=c(2,3,3,3,3,3)) +
  scale_size_manual(values=rep(3.5,6)) +
  ylab("cosine similarity") +
  ggtitle("Pooled, ACS, 2007-2021 (All Tracts)") +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14,angle=90,hjust=0.5),
        axis.title.x=element_text(size=16,hjust=1,vjust=8,margin=margin(t=-12)),
        axis.title.y=element_text(size=16,angle=90),
        plot.title=element_text(size=18,hjust=0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.position="bottom") +
  guides(color=guide_legend(ncol=1), size=guide_legend(direction="vertical"))

ggsave(
  file.path(fig_dir, "acs_M_base_3_tract_pooled_alltracts.png"),
  width=11.5, height=16, units="cm", dpi=800
)
cat("Tract pooled figure saved.\n")

## ── Male and female separate figures ─────────────────────────────────────────
## Local cosine similarity: same logic as cosine_similarity() but writes to a
## local data frame instead of the global cos_sim.
cosine_sim_local <- function(occvec, cs, i, category=22) {
  coslist <- c()
  for (k in 1:category) for (m in 1:category)
    if (k!=m) coslist <- c(coslist, lsa::cosine(occvec[k,], occvec[m,])[1,1])
  for (f in 1:category) for (s in 1:category) {
    cossim <- (lsa::cosine(occvec[f,], occvec[s,])[1,1] - mean(coslist)) / sd(coslist)
    cs[which(cs$year==i), 1+(f-1)*category+s] <- cossim
  }
  cs
}

tract_plot_sex <- function(dat_sex) {
  cs <- data.frame(year=2007:2021)
  for (i in 2007:2021) {
    dat <- dat_sex %>% filter(year==i)
    if (nrow(dat) == 0) next
    olm <- with(dat, questionr::wtd.table(variable, GEOID, estimate)) %>%
      prop.table(margin=2)
    olm <- scale(olm, center=TRUE, scale=TRUE); olm[is.na(olm)] <- 0
    sv     <- irlba(Matrix(olm, sparse=TRUE), nu=10, nv=10)
    occvec <- sv$u %*% diag(sv$d[1:10])
    cs     <- cosine_sim_local(occvec, cs, i, category=22)
    print(paste("year", i, "done"))
  }
  cs[, c(1, 12, 13, 15, 16, 17, 18)] %>%
    pivot_longer(-year, names_to="variable", values_to="value") %>%
    mutate(variable = recode(variable,
      V12="Healthcare Support", V13="Protective Service",
      V15="Building Maintenance", V16="Personal Care and Service",
      V17="Sales and Related", V18="Office and Admin Support"),
      variable = factor(variable, levels=c(
        "Personal Care and Service","Protective Service","Building Maintenance",
        "Healthcare Support","Sales and Related","Office and Admin Support")))
}

for (sx in c("male","female")) {
  d <- pooled %>% filter(year >= 2007, sex == sx) %>%
    group_by(GEOID, year, variable) %>%
    summarize(estimate=sum(estimate,na.rm=TRUE), .groups="drop")
  ttl  <- if (sx=="male") "Male, ACS, 2007-2021 (All Tracts)" else "Female, ACS, 2007-2021 (All Tracts)"
  fout <- if (sx=="male") "acs_M_base_3_tract_male.png" else "acs_M_base_3_tract_female.png"

  df_plot <- tract_plot_sex(d)

  ggplot(df_plot,
         aes(x=year, y=value, group=variable, color=variable, fill=variable,
             shape=variable, lty=variable)) +
    geom_point(aes(size=variable)) +
    geom_line() +
    scale_x_continuous(breaks=seq(2007,2021,2)) +
    scale_y_continuous(limits=c(-2.2,2.2), breaks=seq(-2.2,2.2,0.7)) +
    theme_classic() +
    scale_colour_manual(values=c("#82243b","grey30","grey30","#82243b","#ed5278","grey30")) +
    scale_fill_manual(values=c("#82243b","#ed5278","#ffb6c1",NA,NA,NA)) +
    scale_shape_manual(values=rep(21,6)) +
    scale_linetype_manual(values=c(2,3,3,3,3,3)) +
    scale_size_manual(values=rep(3.5,6)) +
    ylab("cosine similarity") +
    ggtitle(ttl) +
    theme(axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14,angle=90,hjust=0.5),
          axis.title.x=element_text(size=16,hjust=1,vjust=8,margin=margin(t=-12)),
          axis.title.y=element_text(size=16,angle=90),
          plot.title=element_text(size=18,hjust=0.5),
          legend.title=element_blank(),
          legend.text=element_text(size=16),
          legend.position="bottom") +
    guides(color=guide_legend(ncol=1), size=guide_legend(direction="vertical"))

  ggsave(file.path(fig_dir, fout), width=11.5, height=16, units="cm", dpi=800)
  cat(fout, "saved.\n")
}


## ── Section 6: Table 1 ACS column (Panel A lead-lag + Panel B Bartik) ────────
##
## Build CZ-year panel from acs_agg (already in memory): for each year 2005-2021,
## run global SVD on 5-yr window occ-CZ matrix, then compute CZ-level PM and PS
## centroids and their cosine similarity.

nu          <- 50
vec_cols_ll <- paste0("V", 1:nu)
years_panel <- 2005:2021
da_all      <- c(1L, 3L, 4L, 5L, 6L, 7L, 8L, 9L)
min_pm      <- 30; min_ps <- 10

panel <- list()
for (yr in years_panel) {
  cat("Year", yr, "... ")
  sub <- acs_agg %>% ungroup() %>%
    filter(year >= yr-2, year <= yr+2, da %in% da_all)

  olm <- questionr::wtd.table(sub$occ1990, sub$metarea, sub$wtfinl) %>%
    prop.table(margin = 2)
  olm <- scale(olm, center = TRUE, scale = TRUE); olm[is.na(olm)] <- 0

  occ_names <- rownames(olm)
  sv        <- svd(olm, nu = nu, nv = nu)
  occ_vecs  <- t(apply(sv$u %*% diag(sv$d[1:nu]), 1, l2.norm))
  rownames(occ_vecs) <- occ_names

  occ_df <- as.data.frame(occ_vecs)
  colnames(occ_df) <- vec_cols_ll
  occ_df$occ1990 <- as.integer(occ_names)

  cz_occ <- sub %>%
    group_by(metarea, occ1990, da) %>%
    summarize(emp = sum(wtfinl, na.rm = TRUE), .groups = "drop") %>%
    mutate(occ1990 = as.integer(occ1990)) %>%
    left_join(occ_df, by = "occ1990")

  centroids <- cz_occ %>%
    filter(da %in% c(1L, 9L)) %>%
    group_by(metarea, da) %>%
    summarize(across(all_of(vec_cols_ll), ~weighted.mean(., w = emp, na.rm = TRUE)),
              emp_total = sum(emp), n_occ = n(), .groups = "drop") %>%
    filter(n_occ >= 2)

  pm_ct  <- centroids %>% filter(da == 1L, emp_total >= min_pm)
  ps_ct  <- centroids %>% filter(da == 9L, emp_total >= min_ps)
  common <- as.character(intersect(pm_ct$metarea, ps_ct$metarea))

  pm_ct <- pm_ct %>% filter(as.character(metarea) %in% common) %>% arrange(metarea)
  ps_ct <- ps_ct %>% filter(as.character(metarea) %in% common) %>% arrange(metarea)

  pm_mat   <- as.matrix(pm_ct[, vec_cols_ll])
  ps_mat   <- as.matrix(ps_ct[, vec_cols_ll])
  cos_vals <- rowSums(pm_mat * ps_mat) /
    (sqrt(rowSums(pm_mat^2)) * sqrt(rowSums(ps_mat^2)))

  tot_emp <- cz_occ %>%
    group_by(metarea) %>%
    summarize(total_emp = sum(emp), .groups = "drop") %>%
    filter(metarea %in% common) %>% arrange(metarea)

  panel[[length(panel)+1]] <- data.frame(
    czone     = common,
    year      = yr,
    cos_sim   = cos_vals,
    pm_share  = pm_ct$emp_total / tot_emp$total_emp,
    ps_share  = ps_ct$emp_total / tot_emp$total_emp
  )
  cat(sprintf("%d CZs, mean cos_sim=%.3f\n", length(common), mean(cos_vals)))
}

panel_acs <- bind_rows(panel) %>%
  mutate(czone    = as.character(czone),
         year     = as.integer(year),
         cos_sim  = as.numeric(cos_sim),
         pm_share = as.numeric(pm_share),
         ps_share = as.numeric(ps_share)) %>%
  arrange(czone, year)

## ── Panel A: 4-window stacked lead-lag ───────────────────────────────────────
bal_acs <- panel_acs %>%
  filter(czone %in% (panel_acs %>% dplyr::count(czone) %>%
    filter(n == n_distinct(panel_acs$year)) %>% pull(czone))) %>%
  arrange(czone, year)

all_yrs_a <- sort(unique(bal_acs$year))
cuts_a    <- round(seq(1, length(all_yrs_a)+1, length.out=5))
stacked_acs <- do.call(rbind, lapply(1:4, function(i) {
  bal_acs %>% filter(year %in% all_yrs_a[cuts_a[i]:(cuts_a[i+1]-1)]) %>%
    group_by(czone) %>%
    summarize(cos_sim=mean(cos_sim), pm_share=mean(pm_share),
              ps_share=mean(ps_share), .groups="drop") %>% mutate(window=i)
})) %>% arrange(czone, window) %>% group_by(czone) %>%
  mutate(d_cos=cos_sim-lag(cos_sim), d_pm=pm_share-lag(pm_share),
         d_ps=ps_share-lag(ps_share), l_d_pm=lag(d_pm), l_d_ps=lag(d_ps)) %>%
  ungroup()

acs_ll_pm <- feols(d_cos ~ l_d_pm | czone + window, data=stacked_acs, vcov=~czone)
acs_ll_ps <- feols(d_cos ~ l_d_ps | czone + window, data=stacked_acs, vcov=~czone)
acs_n_ll  <- n_distinct(stacked_acs$czone[!is.na(stacked_acs$l_d_pm)])
cat("\nTable 1 Panel A (ACS):\n")
cat(sprintf("  Lag DeltaPM: %.3f (%.3f) t=%.2f\n",
            coef(acs_ll_pm)["l_d_pm"], se(acs_ll_pm)["l_d_pm"],
            coef(acs_ll_pm)["l_d_pm"]/se(acs_ll_pm)["l_d_pm"]))
cat(sprintf("  Lag DeltaPS: %.3f (%.3f) t=%.2f\n",
            coef(acs_ll_ps)["l_d_ps"], se(acs_ll_ps)["l_d_ps"],
            coef(acs_ll_ps)["l_d_ps"]/se(acs_ll_ps)["l_d_ps"]))
acs_ll <- acs_ll_pm

## ── Panel B: Bartik IV ────────────────────────────────────────────────────────
cat("Building ACS Bartik instrument...\n")
t0 <- 2005:2009; t1 <- 2017:2021

d_pm_acs <- acs_agg %>% ungroup() %>%
  filter(year %in% c(t0, t1), da %in% da_all) %>%
  mutate(period = if_else(year %in% t0, "t0", "t1")) %>%
  group_by(metarea, period) %>%
  summarize(pm  = sum(wtfinl[da == 1L], na.rm = TRUE),
            tot = sum(wtfinl, na.rm = TRUE), .groups = "drop") %>%
  mutate(pm_sh = pm / tot) %>%
  select(czone = metarea, period, pm_sh) %>%
  pivot_wider(names_from = period, values_from = pm_sh) %>%
  filter(!is.na(t0), !is.na(t1)) %>%
  mutate(d_pm = t1 - t0)

pm_occ0 <- acs_agg %>% ungroup() %>%
  filter(da == 1L, year %in% t0) %>%
  group_by(czone = metarea, occ1990) %>%
  summarize(e0 = sum(wtfinl, na.rm = TRUE), .groups = "drop") %>%
  group_by(czone) %>% mutate(sh0 = e0 / sum(e0)) %>% ungroup()
pm_occ1 <- acs_agg %>% ungroup() %>%
  filter(da == 1L, year %in% t1) %>%
  group_by(czone = metarea, occ1990) %>%
  summarize(e1 = sum(wtfinl, na.rm = TRUE), .groups = "drop")
pm_nat  <- pm_occ0 %>% left_join(pm_occ1, by = c("czone","occ1990")) %>%
  mutate(e1 = if_else(is.na(e1), 0, e1)) %>%
  group_by(occ1990) %>%
  summarize(n0 = sum(e0), n1 = sum(e1), .groups = "drop")
bartik_acs <- pm_occ0 %>%
  left_join(pm_occ1, by = c("czone","occ1990")) %>%
  mutate(e1 = if_else(is.na(e1), 0, e1)) %>%
  left_join(pm_nat, by = "occ1990") %>%
  mutate(loo_g = ((n1 - e1) - (n0 - e0)) / (n0 - e0)) %>%
  filter(is.finite(loo_g)) %>%
  group_by(czone) %>%
  summarize(bartik_pm = sum(sh0 * loo_g, na.rm = TRUE), .groups = "drop")

cos_ld_acs <- panel_acs %>%
  filter(year %in% c(t0, t1)) %>%
  mutate(period = if_else(year %in% t0, "t0", "t1")) %>%
  group_by(czone, period) %>%
  summarize(cos = mean(cos_sim), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = cos) %>%
  filter(!is.na(t0), !is.na(t1)) %>%
  mutate(d_cos = t1 - t0)

acs_iv_df <- d_pm_acs %>%
  inner_join(cos_ld_acs %>% select(czone, d_cos), by = "czone") %>%
  inner_join(bartik_acs, by = "czone")
cat(sprintf("ACS IV sample: %d CZs\n", nrow(acs_iv_df)))

acs_fs  <- feols(d_pm ~ bartik_pm,              data = acs_iv_df, vcov = "HC1")
acs_ols <- feols(d_cos ~ d_pm,                  data = acs_iv_df, vcov = "HC1")
acs_iv  <- feols(d_cos ~ 1 | d_pm ~ bartik_pm,  data = acs_iv_df, vcov = "HC1")

cat("\nTable 1 Panel B (ACS):\n")
cat(sprintf("  First-stage F = %.1f\n", fitstat(acs_fs, "f")$f$stat))
etable(acs_ols, acs_iv, se.below = TRUE, fitstat = ~ivf,
       headers = c("OLS: DeltaPM->Deltacos", "2SLS: DeltaPM->Deltacos"))

rm(panel, pm_occ0, pm_occ1, pm_nat, bartik_acs, d_pm_acs, cos_ld_acs, acs_iv_df); gc()

## ── Interactive ──────────────────────────────────────────────────────────────
summary(acs_ll_pm)
summary(acs_iv)
