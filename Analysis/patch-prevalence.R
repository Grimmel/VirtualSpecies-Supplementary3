###################################################################################
#
# Author:       Liam Grimmett
# 
# Email:        ligrimmett@csu.edu.au
# 
# Date:         `r paste(Sys.Date())`
#
# Script Name:  Patch prevalence
#
# Description:  Calculates the number of cells outside high suitability patches
#               that are occupied. This is to determine the prevalence of each
#               species within and outside of habitat patches.
#
# Notes:        None          
#
# Copyright (c) Liam Grimmett, 2020
###################################################################################

library(raster)
species <-c("Species2","Species3","Species4","Species5") #
exptFolder <- "D:\\PHDExperimentOutputs\\MainSims\\"
landscapeCategories <- c('LLM','LLS','LMM','LMS','LSS')
landscapeReplicates <- 10
replicates <- 9
nSims <- 2
# Looping through all my simulation types, and landscapes
for (spec in species){
  # Data frame to store each correlogram
  for (i in 2:nSims){
    for(lc in landscapeCategories){
      for (lr in 1:landscapeReplicates){
        suit <- raster(paste('D:\\PHDExperimentOutputs\\SimLandscapes\\suitability\\',lc,lr,'_suitability.asc',sep=''))
        # Convert landscape so high suitability patches are 1 and low suitabilty are 0. 
        suit[suit<75,]<-0
        suit[suit>0,] <- 1
        ls <- paste(lc,lr,sep='')
        for (j in 0:9){
          if(file.exists(paste(exptFolder,spec,'\\Output_Maps\\abundance\\abundance_s',i,'_',ls,'_r',j,'.tif',sep=''))){
            # Load abundance raster
            pa <- raster(paste(exptFolder,spec,'\\Output_Maps\\abundance\\abundance_s',i,'_',ls,'_r',j,'.tif',sep=''))
            # Convert to occurence
            pa[is.na(pa)]<-0
            pa[pa>0]<-1
            # Sum the number of occurences within each zone (patch/non-patch)
            occ_out <- zonal(pa,suit,fun='sum')
            # patch_sum = sum within patches
            # out_sum = sum outside patches
            results <- data.frame(species=spec,lc=lc,lr=lr,ls=ls,sim=i,rep=j,patch_sum=occ_out[[2,2]],out_sum=occ_out[[1,2]])
            write.table(results,paste(exptFolder,'Analysis\\Structure\\prevalence_75.csv',sep=''),append=TRUE,col.names=!file.exists(paste(exptFolder,'Analysis\\Structure\\prevalence_suit_ibm75.csv',sep='')),sep=',')
          }
        }
      }
    }
  }
}
