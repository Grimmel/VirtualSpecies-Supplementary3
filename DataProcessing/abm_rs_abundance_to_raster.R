library(raster)

exptFolder <- "D:\\PHDExperimentOutputs\\MainSims\\Pop_1\\"
landscapeCategories <- c('LLM','LLS','LMM','LMS','LSS')
landscapeReplicates <- 10
replicates <- 9
nSims <- 2
sp.raster <- raster(paste("D:\\PHDExperimentOutputs\\MainSims\\Species5\\Inputs\\LLM1_suitability.txt",sep=""))


for (i in 2:nSims){
  counter <- 1
  for(lc in landscapeCategories){
    for (lr in 1:landscapeReplicates){
      ls <- paste(lc,lr,sep='')
      
      occupancyfileName <- paste(exptFolder,'Outputs\\','Sim',i,'_land',counter,'_Occupancy.txt',sep='')
      occupancy <- read.table(occupancyfileName,h = T, sep = "\t")
      occ <- rasterize(occupancy[, c("x", "y")], y = sp.raster, field = occupancy[, "Year_600"])
      writeRaster(occ,paste(exptFolder,'Output_Maps\\occupancy\\occupancy_s',i,'_',ls,'.tif',sep=''),overwrite=TRUE)
      for (j in 0:replicates){
        populationfileName <- paste(exptFolder,'Outputs\\','Sim',i,'_land',counter,'_Pop.txt',sep='')
        population <- read.table(populationfileName,h = T, sep = "\t")
        population <- population[population$Rep==j & population$Year ==600,]
        if(nrow(population)>0){
          abundance <- rasterize(population[, c("x", "y")], y = sp.raster, field = population[, "NInd"])
          pa <- abundance
          pa[pa>0] <- 1
          writeRaster(abundance,paste(exptFolder,'Output_Maps\\abundance\\abundance_s',i,'_',ls,'_r',j,'.tif',sep=''),overwrite=TRUE)
        }
      }
      counter <- counter + 1
    }
  }
}

