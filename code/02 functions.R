## 02 functions.R
## Helper functions used by all analysis scripts.
## Source this at the top of each analysis script: source(file.path(code_dir, "02 functions.R"))

## ── Packages ─────────────────────────────────────────────────────────────────
library(dplyr)
library(haven)
library(sjlabelled)
library(ggplot2)
library(zoo)
library(MASS)
library(lsa)
library(sda)
library(ggpubr)
library(gridExtra)
library(data.table)
library(questionr)
library(Rtsne)
library(textshape)
library(akima)
library(tidyverse)
library(emdist)
library(rlang)

## 5-year moving window to read data (original data)
five_rd <- function(x,year){
  return(x[which(x$year>=year-2&
                   x$year<=year+2),])
}

## 5-year moving window to smooth data (occupation weight)
five_mw <- function(x){
  return(c((x-2):(x+2)))
}

## 3-year moving average (smooth plot)
three_ma <- function(x) c(x[1], rollmean(x, 3), x[length(x)])

## normalize vector function
l2.norm <- function(x){x/norm(x,type="2")}

## cosine similarity between occupation vector pairs for each year
cosine_similarity <- function(data=occvec,
                              category=6,
                              standardize=TRUE){

  coslist <- c()
  for (k in 1:category){
    for (m in 1:category){
      if (k!=m){
        coslist <- append(coslist,cosine(occvec[k,],occvec[m,])[1,1])
      } else {}
    }
  }

  ## if standardize
  ## rescale cosine similarity by the mean cosine similarity between all possible class pairs
  if (standardize==TRUE){
    ## calculate cosine similarity for all possible combinations of classes
    for (f in 1:category){
      for (s in 1:category){
        cossim <- (cosine(occvec[f,],occvec[s,])[1,1] - mean(coslist))/sd(coslist) ## center and scale
        cos_sim[which(cos_sim$year==i),1+(f-1)*category+s] <<- cossim
      }
    }
  } else { ## if not standardize
    for (f in 1:category){
      for (s in 1:category){
        cossim <- cosine(occvec[f,],occvec[s,])[1,1]
        cos_sim[which(cos_sim$year==i),1+(f-1)*category+s] <<- cossim
      }
    }
  }
}

## create weighted count by year and occupation group
create_count_group <- function(data=cps,
                          crosswalk="egp"){

  if (crosswalk=="egp"){
    ## calculate group count by year
    classshare <-
      five_rd(data,i) %>%
      filter(egp!="IVa"&egp!="IVb"&egp!="IVc") %>%
      mutate(egp = case_when(egp %in% c("IIIa","IIIb") ~ "III",
                             egp %in% c("VIIa","VIIb") ~ "VII",
                             TRUE ~ as.character(egp))) %>%
      filter(!is.na(wtfinl)) %>%
      group_by(metarea,egp) %>%
      summarize(Nbyegp=sum(wtfinl)) %>% ungroup() %>%
      group_by(metarea) %>%
      mutate(Ntotal=sum(Nbyegp),
             share=Nbyegp/Ntotal,
             count=n())

    share <- classshare %>%
      ungroup() %>%
      dplyr::select(metarea,egp,share)
    share <- dcast(share,metarea~egp)
    return(share)


  } else if (crosswalk=="da") {
    ## calculate group count by year
    classshare <-
      five_rd(data,i) %>%
      mutate(da = case_when(da=="managers/executives" ~ "1",
                            da=="professionals" ~ "2",
                            da=="technicians" ~ "3",
                            da=="sales" ~ "4",
                            da=="administrative/office" ~ "5",
                            da=="production" ~ "6",
                            da=="laborers" ~ "7",
                            da=="clean and protect services" ~ "8",
                            da=="personal services" ~ "9",
                            TRUE ~ as.character(da))) %>%
      mutate(da=as.integer(da)) %>%
      filter(!is.na(wtfinl)) %>%
      group_by(metarea,da) %>%
      summarize(Nbyegp=sum(wtfinl)) %>% ungroup() %>%
      group_by(metarea) %>%
      mutate(Ntotal=sum(Nbyegp),
             share=Nbyegp/Ntotal,
             count=n())

    share <- classshare %>%
      ungroup() %>%
      dplyr::select(metarea,da,share)
    share <- dcast(share,metarea~da)
    return(share)


  } else if (crosswalk=="ind") {

    ## calculate group count by year (industry)
    classshare <-
      five_rd(data,i) %>%
      merge(occind,by="occ_code",all.x=T)%>%
      group_by(metarea,ind) %>%
      summarize(Nbyegp=sum(tot_emp)) %>% ungroup() %>%
      group_by(metarea) %>%
      mutate(Ntotal=sum(Nbyegp),
             share=Nbyegp/Ntotal,
             count=n())

    share <- classshare %>%
      ungroup() %>%
      dplyr::select(metarea,ind,share)
    share <- dcast(share,metarea~ind)
    return(share)

  } else {
    ## calculate group count by year (GWmeso)
    classshare <-
      five_rd(data,i) %>%
    filter(!is.na(wtfinl)) %>%
      group_by(metarea,gw_meso) %>%
      summarize(Nbyegp=sum(wtfinl)) %>% ungroup() %>%
      group_by(metarea) %>%
      mutate(Ntotal=sum(Nbyegp),
             share=Nbyegp/Ntotal,
             count=n())

    share <- classshare %>%
      ungroup() %>%
      dplyr::select(metarea,gw_meso,share)
    share <- dcast(share,metarea~gw_meso)
    return(share)
  }
  }

## create centroid with weights
occvec_centroid <- function(data=occvec){
  occvec <- as.data.table(occvec)
  occvec <- occvec[,lapply(.SD,weighted.mean,w=Nocc),
                   by=eval(colnames(occvec)[length(colnames(occvec))-1])]
  occvec <- occvec %>% arrange(!!sym(colnames(occvec)[1])) %>% dplyr::select(-c(Nocc,colnames(occvec)[1]))
  occvec <- as.matrix(occvec)
  return(occvec)
  }

## create occupation vectors
create_occvec <- function(data=cps,
                          nu=200,
                          crosswalk=occegp,
                          cps=TRUE){

  ## create Occupation-Location Matrix in each year

  ## TF
  olm <- with(five_rd(data,i), questionr::wtd.table(occ1990, metarea, wtfinl)) %>%
    prop.table(margin=2)
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

## create occupation vectors for OEWS
create_occvec_owes <- function(data=oews,
                               nu=200,
                               crosswalk=occind){

  ## calculate number of people (weighted) in each occupation in each year for later centroid weight
   Nocc <-
     data %>% group_by(year,occ_code) %>%
     summarize(Nocc=sum(tot_emp,na.rm=T))

  ## create Occupation-Location Matrix in each year
  olm <- with(five_rd(data,i),
              questionr::wtd.table(occ_code, metarea, weights=tot_emp)) %>%
    prop.table(margin=2)
  olm <- scale(olm,center=TRUE,scale=TRUE) ## Occ-Loc Matrix is centered and scaled by default

  ## SVD decomposition and a vector representation of occupations
  svd <- svd(olm,nu=nu)
  occvec <- svd$u %*% diag(svd$d)[1:nu,1:nu]

  ## create class identifier and add occupation weight
  occvec <- t(apply(occvec,1,l2.norm))
  occvec <- as.data.frame(occvec)
  occvec$occ_code <- rownames(olm)

  occvec <- merge(occvec,crosswalk[,c("occ_code","ind")],by="occ_code",all.x=T)
  occvec <- merge(occvec,Nocc[which(Nocc$year %in% five_mw(i)),]
                  %>% group_by(occ_code) %>%
                    summarize(Nocc=sum(Nocc,na.rm=T)),
                  by="occ_code",all.x=T)
  occvec <- occvec[,-which(names(occvec) %in% c("occ_code"))]

  return(occvec)

}

## scaleFUN
scaleFUN <- function(x) sprintf("%.3f", x)
