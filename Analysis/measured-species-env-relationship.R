###################################################################################
#
# Author:       Liam Grimmett
# 
# Email:        ligrimmett@csu.edu.au
# 
# Date:         `r paste(Sys.Date())`
#
# Script Name:  Measured species-environment relationship
#
# Description:  Measure the mean response of the species along the suitability 
#               gradient.
#
# Notes:        None          
#  
# Copyright (c) Liam Grimmett, 2020
###################################################################################
library(raster)

species <- c('EO_Species1','EO_Species2','DO_Species1','DO_Species2','DO_Species3','DO_Species4','DP_Species1','DP_Species2','DP_Species3','DP_Species4')
exptFolder <- "D:\\PHDExperimentOutputs\\MainSims\\"
landscapeCategories <- c('LLM','LLS','LMM','LMS','LSS')
landscapeReplicates <- 10
replicates <- 9

results <- 0
for(spec in species){
  for (lc in landscapeCategories){
    for (lr in 1:landscapeReplicates){
      # Load suitability raster
      suit <- raster(paste('D:\\PHDExperimentOutputs\\SimLandscapes\\suitability\\',lc,lr,'_suitability.asc',sep=''))
      # Discretise suitability into 20 classes
      c <- cut(suit,breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100))
      # Calculate the number of cells within each suitability class
      freqSuit <- freq(c,useNA="no")
      ls <- paste(lc,lr,sep='')
      for(i in 1:1){
        for (j in 0:9){
          if(file.exists(paste(exptFolder,spec,'\\Output_Maps\\abundance\\abundance_s',i,'_',ls,'_r',j,'.tif',sep=''))){
            # Load abundance raster
            ras <- raster(paste(exptFolder,spec,'\\Output_Maps\\abundance\\abundance_s',i,'_',ls,'_r',j,'.tif',sep=''))
            # Convert to occurence
            pa <- ras
            pa[pa>0,] <- 1
            # Zonal statistics by suitability class
            abCount <- zonal(ras,c,fun='sum')       ## Abundance within each suitability class
            abCountPU <- abCount[,2]/freqSuit[,2]   ## Mean abundance for each suitability class
            paCount <- zonal(pa,c,fun='sum')        ## Number of cells occupied within each suitability class
            paCountPU <- paCount[,2]/freqSuit[,2]   ## Proportion of cells occupied for each suitability class
            # Process outputs
            out <- as.data.frame(cbind(freqSuit,abCount[,2],abCountPU,paCount[,2],paCountPU))
            names(out) <- c('break','suit_freq','abun_sum','abun_prop','pa_count','pa_prop')
            out$ls <- ls
            out$lr <- lr
            out$sim <- i
            out$rep <- j
            out$species <- spec
            if(results==0){
              results = out
            }else{
              results <- rbind(results,out)
           }
          }
        }
      }
    }
  }
}
write.csv(results,"D:\\PHDExperimentOutputs\\MainSims\\Analysis\\Structure\\structure.csv")