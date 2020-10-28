library(raster)

exptFolder <- "D:\\PHDExperimentOutputs\\MainSims\\NullP\\"
landscapeCategories <- c('LLM','LLS','LMM','LMS','LSS')
landscapeReplicates <- 10
replicates <- 9

toAb <- function(prob.raster)
{
  calc(prob.raster, fun = function (x)
  {
    sapply(x, FUN = function(y)
    {
      if(is.na(y))
      { NA } else
      {
        rpois(1, 10*(y))
      }
    }
    )
  })
}

for(lc in landscapeCategories){
  for (lr in 1:landscapeReplicates){
    r <- NULL
    ls <- paste(lc,lr,sep='')
    suitability <- raster(paste('D:\\PHDExperimentOutputs\\MainSims\\NullLP\\Inputs\\',ls,'_suitability.txt',sep=''))
    # 10 = 5 when nulllp
    suitability <- exp(10*((suitability/100)-1))
    for (j in 0:9){
      abundance <- toAb(suitability)
      writeRaster(abundance,paste('D:\\PHDExperimentOutputs\\MainSims\\NullVLP\\Output_Maps\\abundance\\abundance_s1_',ls,'_r',j,'.tif',sep=''),overwrite=TRUE)
    }
  }
}

