###################################################################################
#
# Author:       Liam Grimmett
# 
# Email:        ligrimmett@csu.edu.au
# 
# Date:         `r paste(Sys.Date())`
#
# Script Name:  Occurence - Suitability Relationship
#
# Description:  Perform logistic regression using occurence ~ sutiablity to describe
#               the species responses relationship with suitability. Moran's I 
#               correlograms calculated to describe spatial structure of the species
#               response unexplainable by suitability alone.
#
# Notes:        None          
# 
# Copyright (c) Liam Grimmett, 2020
###################################################################################

library(raster)
library(elsa)

species <- c('EO_Species1','EO_Species2','DO_Species1','DO_Species2','DO_Species3','DO_Species4','DP_Species1','DP_Species2','DP_Species3','DP_Species4')
exptFolder <- "D:\\PHDExperimentOutputs\\MainSims\\"
landscapeCategories <- c('LLM','LLS','LMM','LMS','LSS')
landscapeReplicates <- 10
replicates <- 9
nSims <- 1
stratRaster <- raster('D:\\PHDExperimentOutputs\\landscapesstratified400.asc')

# Looping through all my simulation types, and landscapes
t <- Sys.time()
for (spec in species){
  # Data frame to store each correlogram
  df <- data.frame(matrix(ncol=35,nrow=0))
  setNames(df,c('ls','lsrep','species','sim','simrep','int','int_se','int_z','int_p','suit','suit_se','suit_z','suit_p','null_dev','dev'))
  for (i in 1:nSims){
    for(lc in landscapeCategories){
      for (lr in 1:landscapeReplicates){
        suit <- raster(paste('D:\\PHDExperimentOutputs\\SimLandscapes\\suitability\\',lc,lr,'_suitability.asc',sep=''))
        ls <- paste(lc,lr,sep='')
        for (j in 0:9){
          if(file.exists(paste(exptFolder,spec,'\\Output_Maps\\abundance\\abundance_s',i,'_',ls,'_r',j,'.tif',sep=''))){
            # Load abundance
            ras <- raster(paste(exptFolder,spec,'\\Output_Maps\\abundance\\abundance_s',i,'_',ls,'_r',j,'.tif',sep=''))
            ras[is.na(ras)]<-0
            # Convert to occurence
            pa <- ras
            pa[pa>0,]<-1 
            # Prepare data for regression
            # NOTE: All cells are being converted to points. Subsequent regression then uses ALL available simulated data.
            s <- stack(pa,suit)
            abunSample <- as.data.frame(s,xy=TRUE)
            names(abunSample) <- c('x','y','pa','suitability')
            # Run regression
            m <- glm(pa~suitability,data=abunSample,family=binomial())
            m_summary <- summary(m)
            # Output coefficients
            results <- c(m_summary$coefficients[1,1],m_summary$coefficients[1,2],m_summary$coefficients[1,3],m_summary$coefficients[1,4],m_summary$coefficients[2,1],m_summary$coefficients[2,2],m_summary$coefficients[2,3],m_summary$coefficients[2,4],m_summary$null.deviance,m_summary$deviance)
            # Convert residuals into raster
            res <- cbind(abunSample[1:2],resid(m,type="pearson"))
            res <- rasterFromXYZ(res,res=c(100,100))
            # Correlogram
            c <- correlogram(res,500,10000)
            # Add a new row to my data frame
            df[nrow(df)+1,] <- c(lc,lr,spec,i,j,c@correlogram$moran,results)
          }
        }
      }
    }
  }
  write.csv(df,paste(exptFolder,spec,'_lrpa_output.csv',sep=''))
}

