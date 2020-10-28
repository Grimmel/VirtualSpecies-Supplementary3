###################################################################################
#
# Author:       Liam Grimmett
# 
# Email:        ligrimmett@csu.edu.au
# 
# Date:         `r paste(Sys.Date())`
#
# Script Name:  Landscape Metrics
#
# Description:  Calculate landscape metrics for high suitability habitat patches
#
# Notes:        None          
# 
# Copyright (c) Liam Grimmett, 2020
###################################################################################

library(raster)
library(landscapemetrics)

landscapeCategories <- c('LLM','LLS','LMM','LMS','LSS')
landscapeReplicates <- 10

df <- data.frame(matrix(ncol=18,nrow=0))
for(lc in landscapeCategories){
  for(i in 1:landscapeReplicates){
    r <- raster(paste('D:\\PHDExperimentOutputs\\SimLandscapes\\suitability\\',lc,i,'_suitability.asc',sep=''))
    rdat <- as.data.frame(r,xy=TRUE)
    names(rdat) <- c('x','y','suit')
    coordinates(rdat) = ~x + y
    # Patch considered as having suitability value >= 75/85
    core_75 <- r
    core_75[core_75<85] <- 0 #SET as 85 or 75 depending on patch definition
    core_75[core_75>0] <- 1
    # Calculate metrics
    area <- lsm_c_area_mn(core_75, directions = 8)
    area_sd <- lsm_c_area_sd(core_75, directions = 8)
    area_cv <- lsm_c_area_cv(core_75, directions = 8)
    core_area <- lsm_c_core_mn(core_75, directions = 8)
    core_area_sd <- lsm_c_core_sd(core_75, directions = 8)
    core_area_cv <- lsm_c_core_cv(core_75, directions = 8)
    total_core <- lsm_c_tca(core_75,directions=8)
    para_mn <- lsm_c_para_mn(core_75, directions = 8)
    para_sd <- lsm_c_para_sd(core_75, directions = 8)
    para_cv <- lsm_c_para_cv(core_75, directions = 8)
    con_mn <- lsm_c_contig_mn(core_75, directions = 8)
    con_sd <- lsm_c_contig_sd(core_75, directions = 8)
    con_cv <- lsm_c_contig_cv(core_75, directions = 8)
    coh <- lsm_c_cohesion(core_75, directions = 8)
    lpi <- lsm_c_lpi(core_75,directions=8)
    lsi <- lsm_c_lsi(core_75,directions=8)
    tot <- cellStats(core_75,stat='sum')
    df[nrow(df)+1,] = c(paste(lc,i,sep=''),area[2,]$value,area_sd[2,]$value,area_cv[2,]$value,core_area[2,]$value,core_area_sd[2,]$value,core_area_cv[2,]$value,total_core[2,]$value,para_mn[2,]$value,para_sd[2,]$value,para_cv[2,]$value,con_mn[2,]$value,con_sd[2,]$value,con_cv[2,]$value,lpi[2,]$value,lsi[2,]$value,tot,coh[2,]$value)
  }
}
names(df) <- c('ls','area_mn','area_sd','area_cv','core_mn','core_sd','core_cv','total_core','para_mn','para_sd','para_cv','con_mn','con_sd','con_cv','lpi','lsi','total_area','cohesion')
write.csv(df,'D:\\PHDExperimentOutputs\\SimLandscapes\\suitability\\landscape_metrics85.csv')