###################################################################################
#
# Author:       Liam Grimmett
# 
# Email:        ligrimmett@csu.edu.au
# 
# Date:         `r paste(Sys.Date())`
#
# Script Name:  Patch Measures
#
# Description:  Calculates landscape metrics for each patch in each landscape. Then
#               measures the prevalence of the species within each patch.
#
# Notes:        None          
# 
# Copyright (c) Liam Grimmett, 2020
###################################################################################
library(raster)
library(landscapemetrics)

species <- c('EO_Species1','EO_Species2','DO_Species1','DO_Species2','DO_Species3','DO_Species4','DP_Species1','DP_Species2','DP_Species3','DP_Species4')
exptFolder <- "D:\\PHDExperimentOutputs\\MainSims\\"
landscapeCategories <- c('LLM','LLS','LMM','LMS','LSS')
landscapeReplicates <- 10
replicates <- 9

results <- 0

for (lc in landscapeCategories){
  print(lc)
  for (lr in 1:landscapeReplicates){
    ls <- paste(lc,lr,sep='')
    suit <- raster(paste('D:\\PHDExperimentOutputs\\SimLandscapes\\suitability\\',lc,lr,'_suitability.asc',sep=''))
    patch <- suit
    patch[patch<85,]<-NA # Change to 75 for high suitability patches
    patch[patch>0,] <- 1
    patches <- get_patches(patch,class=1)
    # show_patches(pa)
    patch_area <- lsm_p_area(patch)
    patch_core <- lsm_p_core(patch)
    patch_contig <- lsm_p_contig(patch)
    patch_perim <-lsm_p_perim(patch)
    patch_enn <- lsm_p_enn(patch)
    patch_frac <- lsm_p_frac(patch)
    patch_mean_suit <- zonal(suit,patches$`1`,fun='mean')
    patch_sd_suit <- zonal(suit,patches$`1`,fun='sd',na.rm=TRUE)
    
    patch_dat <- as.data.frame(cbind(patch_area$id,patch_area$value,patch_core$value,patch_contig$value,patch_perim$value,patch_enn$value,patch_frac$value,patch_mean_suit[,2],patch_sd_suit[,2]))
    names(patch_dat) <- c('patch_id','area','core_area','contig','perim','enn','frac','mean','sd')

    for(spec in species){
      print(spec)
      for(i in 1:2){
        for (j in 0:9){
          if(file.exists(paste(exptFolder,spec,'\\Output_Maps\\abundance\\abundance_s',i,'_',ls,'_r',j,'.tif',sep=''))){
            ras <- raster(paste(exptFolder,spec,'\\Output_Maps\\abundance\\abundance_s',i,'_',ls,'_r',j,'.tif',sep=''))
            pa <- ras
            pa[pa>0,] <- 1
            paCount <- zonal(pa,patches$`1`,fun='sum')
            out <- as.data.frame(cbind(paCount[,1],paCount[,2]))
            names(out) <- c('patch_id','pa_count')
            out <- cbind(out,patch_dat)
            out$ls <- ls
            out$lr <- lr
            out$sim <- i
            out$rep <- j
            out$species <- spec
            if(results==0){
              results <- out
            }else{
              results <- rbind(results,out)
            }
          }
        }
      }
    }
  }
}
write.csv(results,"D:\\PHDExperimentOutputs\\MainSims\\Analysis\\Structure\\patches_measures85.csv")