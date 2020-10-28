library(NLMR)
library(raster)
library(landscapemetrics)
library(virtualspecies)
resolution <- 400
nreps <- 10
output <- "D:\\PHDExperimentOutputs\\SimLandscapes\\"

combinations <- list(c(0.25,0.25,0.10,'LLM'),
                  c(0.25,0.25,0.025,'LLS'),
                  c(0.25,0.10,0.10,'LMM'),
                  c(0.25,0.10,0.025,'LMS'),
                  c(0.25,0.025,0.025,'LSS'))

for(i in 1:5){
  vals <- combinations[[i]]
  for(j in 1:nreps){
    landscapeName <- paste(vals[4],j,sep='')
    e1 <- nlm_gaussianfield(resolution, resolution, resolution = 100, autocorr_range = resolution * as.double(vals[1]),
                            mag_var = 2, nug = 0.1, mean = 1, user_seed = NULL,
                            rescale = TRUE)
    e2 <- nlm_gaussianfield(resolution, resolution, resolution = 100, autocorr_range =resolution * as.double(vals[2]),
                            mag_var = 2, nug = 0.1, mean = 1, user_seed = NULL,
                            rescale = TRUE)
    e3 <- nlm_gaussianfield(resolution, resolution, resolution = 100, autocorr_range =resolution * as.double(vals[3]),
                            mag_var = 2, nug = 0.1, mean = 1, user_seed = NULL,
                            rescale = TRUE)
    writeRaster(e1,paste(output,'evar\\',landscapeName,'_e1.asc',sep=''),format='ascii',overwrite=TRUE)
    writeRaster(e2,paste(output,'evar\\',landscapeName,'_e2.asc',sep=''),format='ascii',overwrite=TRUE)
    writeRaster(e3,paste(output,'evar\\',landscapeName,'_e3.asc',sep=''),format='ascii',overwrite=TRUE)
    l1_predictors <- stack(e1,e2,e3)
    names(l1_predictors) <- c('temp','precip','habitat')
    # Formatting of the response functions
    species1.parameters <- formatFunctions(temp = c(fun = 'dnorm', mean = 0.65, sd = 0.3),
                                           precip = c(fun = 'dnorm', mean = 0.35, sd = 0.3),
                                           habitat = c(fun = 'dnorm', mean = 0.65, sd = 0.3))
    # Generation of the virtual species
    species1 <- generateSpFromFun(raster.stack = l1_predictors,parameters = species1.parameters,plot=TRUE,rescale=FALSE)
    writeRaster(species1$suitab.raster*100,paste(output,'suitability\\',landscapeName,'_suitability.asc',sep=''),format='ascii',overwrite=TRUE)
  }
}



                              