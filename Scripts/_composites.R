##  ---------------------------------------------------------------------------------------------------------------------  ##
                           # Composite Maps for Functional Groups
##  ---------------------------------------------------------------------------------------------------------------------  ##
# Written by Nicholas J Lyon
  ## nicholasjlyon@gmail.com

# PURPOSE
  ## This script is in response to the reviewers who have asked for 
  ## a synthetic map of all milkweeds predicted to do
  ## better than some threshold (50% suitability) across the landscape

# Get necessary libraries
library(rgbif); library(rbison); library(ecoengine); library(spocc); library(rgdal);
library(rgeos); library(dismo); library(raster); library(scrubr); library(geosphere);
library(scales); library(rJava); library(mapr); library(ggmap); library(ggplot2);
library(sp)

# Check to make sure maxent is running through R
system.file('java', package = 'dismo'); maxent(); .jinit()

# Set working directory & clear environment
setwd("~/Documents/School/1. Iowa State/_MS Project/_Leopold Project/Lyon_etal_2018_SDM_Project")
rm(list = ls())

# Handy shapefiles denoting country/state/ecotone borders
countries <- rgdal::readOGR('./Borders', 'TM_WORLD_BORDERS-0.3')
states <- rgdal::readOGR('./Borders', 'cb_2015_us_state_20m')
eco <- rgdal::readOGR('./Borders', 'provinces', verbose = F)

# Get predictor variable lists for each species
elyvir.pred <- c('WC08', 'WC09', 'WC10', 'WC16', 'WC17')
koemac.pred <- c('WC09', 'WC10', 'WC16', 'WC17', 'WC18')
stispa.pred <- c('WC08', 'WC09', 'WC16', 'WC17')
boucur.pred <- c('WC08', 'WC09', 'WC10', 'WC16', 'WC17', 'WC19')
schsco.pred <- c('WC08', 'WC09', 'WC10', 'WC16', 'WC17')
sornut.pred <- c('WC08', 'WC09', 'WC10', 'WC16', 'WC17')
ascinc.pred <- c('WC08', 'WC09', 'WC10', 'WC16', 'WC17')
ascsyr.pred <- c('WC08', 'WC09', 'WC10', 'WC16', 'WC17')
asctub.pred <- c('WC08', 'WC09', 'WC10', 'WC16', 'WC17')
dryarg.pred <- c('WC08', 'WC10', 'WC11', 'WC16', 'WC17')
lobsip.pred <- c('WC08', 'WC09', 'WC10', 'WC16', 'WC17')
monfis.pred <- c('WC08', 'WC09', 'WC10', 'WC17')
amocan.pred <- c('WC08', 'WC09', 'WC16', 'WC17')
petcan.pred <- c('WC08', 'WC09', 'WC10', 'WC16', 'WC19')

# Set threshold of suitability here and all following code will call this object
thresh <- 0.5

##  ---------------------------------------------------------------------------------------------------------------------  ##
                                    # Ready Climate Data
##  ---------------------------------------------------------------------------------------------------------------------  ##


# Make an object of all climate variables used in the modeling
clim <- stack(c('./Composite Maps/Clim Data/Elevation.tif',
                    list.files('./Composite Maps/Clim Data/Current/', full.names = T)))

##  ---------------------------------------------------------------------------------------------------------------------  ##
                # Load Finalized "Current" Baseline Maps for all Asclepias Species
##  ---------------------------------------------------------------------------------------------------------------------  ##
# Load in finished baseline model
ASP.intermed <- load("./Modeled Spp/ASP/ASP Model V3 - Baseline.Rdata")

# Get the model into a different object
ASP <- baselineModel

# Pull in predictors for this species
ASP.predictors <- c('ASP_WC08', 'ASP_WC09', 'ASP_WC10', 'ASP_WC11', 'ASP_WC16', 'ASP_WC17')

# Subset the full climate data for just these predictor variables
ASP.predclim <- subset(clim, ASP.predictors)

# We'll want it eventually, so let's go ahead and make the predicted map for each species too
ASP.map <- predict(ASP, ASP.predclim, filename = './Composite Maps/Species Maps/ASP - Baseline Map',
                       format = 'GTiff', progress = "text", overwrite = T)

# Get/crop all the climate data for all the species
  ## Naming conventions are dumb because of species-specific naming conventions in baseline models
for (i in 1:19) {
  print(paste('Clipping current WORLDCLIM raster ', i))
  flush.console()
  
  rast <- raster(paste0('./Raw WORLDCLIM/1960-1990/bio', i, '.bil'))
  rast <- crop(rast, studyExtent)
  rast <- setMinMax(rast)
  projection(rast) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  
  writeRaster(rast, paste0('./Composite Maps/Clim Data/Current/ASP_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/Current/ERO_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/Current/FAS_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/Current/INC_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/Current/OEN_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/Current/SPE_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/Current/SYR_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/Current/TUB_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/Current/VER_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/Current/VIR_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
}

# Want to convert it to threshold value too though for later summation across species
ASP.map[ASP.map < thresh] = 0
ASP.map[ASP.map >= thresh] = 1

# Now do that for every other modeled species (comments excluded for others for space efficiency)
# ERO
ERO.intermed <- load("./Modeled Spp/ERO/ERO Model V3 - Baseline.Rdata")
ERO <- baselineModel
ERO.predictors <- c('ERO_WC08', 'ERO_WC09', 'ERO_WC16', 'ERO_WC18')
ERO.predclim <- subset(clim, ERO.predictors)
ERO.map <- predict(ERO, ERO.predclim, 
                   filename = './Composite Maps/Species Maps/ERO - Baseline Map',
                   format = 'GTiff', progress = "text", overwrite = T)
ERO.map[ERO.map < thresh] = 0
ERO.map[ERO.map >= thresh] = 1

# Make a mega map
present.map <- (ASP.map + ERO.map + FAS.map + INC.map + OEN.map + 
  SPE.map + SYR.map + TUB.map + VER.map + VIR.map)

# Plot it!
jpeg(file = "./Composite Maps/Composite Map Pics/current_compmap.jpg", 
     width = 1080, height = 1138, quality = 100)
plot(present.map, col = rev(bpy.colors(7)))
sp::plot(countries, add = T, border = 'gray45')
sp::plot(states, add = T, border = 'gray45')
dev.off()

##  ---------------------------------------------------------------------------------------------------------------------  ##
                                  # Future Synthetic Maps
##  ---------------------------------------------------------------------------------------------------------------------  ##
# Get and crop CC 4.5 Future
for (i in 1:19) {
  print(paste('Clipping CC 4.5 WORLDCLIM raster ', i))
  flush.console()
  
  rast <- raster(paste0('./Raw WORLDCLIM/CC 4.5/cc45bi50', i, '.tif'))
  rast <- crop(rast, studyExtent)
  rast <- setMinMax(rast)
  projection(rast) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 4.5/ASP_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 4.5/ERO_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 4.5/FAS_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 4.5/INC_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 4.5/OEN_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 4.5/SPE_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 4.5/SYR_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 4.5/TUB_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 4.5/VER_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 4.5/VIR_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
}

# Stack them into a raster layer
cc45.clim <- stack(list.files('./Composite Maps/Clim Data/CC 4.5', full.names = T))

# Ditto for the 8.5 future
for (i in 1:19) {
  print(paste('Clipping CC 8.5 WORLDCLIM raster ', i))
  flush.console()
  
  rast <- raster(paste0('./Raw WORLDCLIM/CC 8.5/cc85bi70', i, '.tif'))
  rast <- crop(rast, studyExtent)
  rast <- setMinMax(rast)
  projection(rast) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 8.5/ASP_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 8.5/ERO_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 8.5/FAS_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 8.5/INC_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 8.5/OEN_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 8.5/SPE_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 8.5/SYR_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 8.5/TUB_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 8.5/VER_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  writeRaster(rast, paste0('./Composite Maps/Clim Data/CC 8.5/VIR_WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
}

# Stack 'em
cc85.clim <- stack(list.files('./Composite Maps/Clim Data/CC 8.5', full.names = T))

##  ------------------------------------  ##
  # Get Futures for Asclepias Spp.
##  ------------------------------------  ##
# Get only predictor climate variables for both futures
ASP.cc45.clim <- subset(cc45.clim, ASP.predictors)
ASP.cc85.clim <- subset(cc85.clim, ASP.predictors)

# Get maps for both futures
ASP.cc45.map <- predict(ASP, ASP.cc45.clim,
                               filename = './Composite Maps/Species Maps/ASP - CC 45 Map',
                               format = 'GTiff',  progress = "text", overwrite = T)
ASP.cc85.map <- predict(ASP, ASP.cc85.clim,
                      filename = './Composite Maps/Species Maps/ASP - CC 85 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)

# And convert to our now beloved thresholds
ASP.cc45.map[ASP.cc45.map < thresh] = 0
ASP.cc45.map[ASP.cc45.map >= thresh] = 1
ASP.cc85.map[ASP.cc85.map < thresh] = 0
ASP.cc85.map[ASP.cc85.map >= thresh] = 1

# Now, do that for all the other species

# ERO
ERO.cc45.clim <- subset(cc45.clim, ERO.predictors)
ERO.cc85.clim <- subset(cc85.clim, ERO.predictors)
ERO.cc45.map <- predict(ERO, ERO.cc45.clim,
                      filename = './Composite Maps/Species Maps/ERO - CC 45 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
ERO.cc85.map <- predict(ERO, ERO.cc85.clim,
                      filename = './Composite Maps/Species Maps/ERO - CC 85 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
ERO.cc45.map[ERO.cc45.map < thresh] = 0
ERO.cc45.map[ERO.cc45.map >= thresh] = 1
ERO.cc85.map[ERO.cc85.map < thresh] = 0
ERO.cc85.map[ERO.cc85.map >= thresh] = 1

# FAS
FAS.cc45.clim <- subset(cc45.clim, FAS.predictors)
FAS.cc85.clim <- subset(cc85.clim, FAS.predictors)
FAS.cc45.map <- predict(FAS, FAS.cc45.clim,
                      filename = './Composite Maps/Species Maps/FAS - CC 45 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
FAS.cc85.map <- predict(FAS, FAS.cc85.clim,
                      filename = './Composite Maps/Species Maps/FAS - CC 85 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
FAS.cc45.map[FAS.cc45.map < thresh] = 0
FAS.cc45.map[FAS.cc45.map >= thresh] = 1
FAS.cc85.map[FAS.cc85.map < thresh] = 0
FAS.cc85.map[FAS.cc85.map >= thresh] = 1

# INC
INC.cc45.clim <- subset(cc45.clim, INC.predictors)
INC.cc85.clim <- subset(cc85.clim, INC.predictors)
INC.cc45.map <- predict(INC, INC.cc45.clim,
                      filename = './Composite Maps/Species Maps/INC - CC 45 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
INC.cc85.map <- predict(INC, INC.cc85.clim,
                      filename = './Composite Maps/Species Maps/INC - CC 85 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
INC.cc45.map[INC.cc45.map < thresh] = 0
INC.cc45.map[INC.cc45.map >= thresh] = 1
INC.cc85.map[INC.cc85.map < thresh] = 0
INC.cc85.map[INC.cc85.map >= thresh] = 1

# OEN
OEN.cc45.clim <- subset(cc45.clim, OEN.predictors)
OEN.cc85.clim <- subset(cc85.clim, OEN.predictors)
OEN.cc45.map <- predict(OEN, OEN.cc45.clim,
                      filename = './Composite Maps/Species Maps/OEN - CC 45 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
OEN.cc85.map <- predict(OEN, OEN.cc85.clim,
                      filename = './Composite Maps/Species Maps/OEN - CC 85 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
OEN.cc45.map[OEN.cc45.map < thresh] = 0
OEN.cc45.map[OEN.cc45.map >= thresh] = 1
OEN.cc85.map[OEN.cc85.map < thresh] = 0
OEN.cc85.map[OEN.cc85.map >= thresh] = 1

# SPE
SPE.cc45.clim <- subset(cc45.clim, SPE.predictors)
SPE.cc85.clim <- subset(cc85.clim, SPE.predictors)
SPE.cc45.map <- predict(SPE, SPE.cc45.clim,
                      filename = './Composite Maps/Species Maps/SPE - CC 45 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
SPE.cc85.map <- predict(SPE, SPE.cc85.clim,
                      filename = './Composite Maps/Species Maps/SPE - CC 85 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
SPE.cc45.map[SPE.cc45.map < thresh] = 0
SPE.cc45.map[SPE.cc45.map >= thresh] = 1
SPE.cc85.map[SPE.cc85.map < thresh] = 0
SPE.cc85.map[SPE.cc85.map >= thresh] = 1

# SYR
SYR.cc45.clim <- subset(cc45.clim, SYR.predictors)
SYR.cc85.clim <- subset(cc85.clim, SYR.predictors)
SYR.cc45.map <- predict(SYR, SYR.cc45.clim,
                      filename = './Composite Maps/Species Maps/SYR - CC 45 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
SYR.cc85.map <- predict(SYR, SYR.cc85.clim,
                      filename = './Composite Maps/Species Maps/SYR - CC 85 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
SYR.cc45.map[SYR.cc45.map < thresh] = 0
SYR.cc45.map[SYR.cc45.map >= thresh] = 1
SYR.cc85.map[SYR.cc85.map < thresh] = 0
SYR.cc85.map[SYR.cc85.map >= thresh] = 1

# TUB
TUB.cc45.clim <- subset(cc45.clim, TUB.predictors)
TUB.cc85.clim <- subset(cc85.clim, TUB.predictors)
TUB.cc45.map <- predict(TUB, TUB.cc45.clim,
                      filename = './Composite Maps/Species Maps/TUB - CC 45 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
TUB.cc85.map <- predict(TUB, TUB.cc85.clim,
                      filename = './Composite Maps/Species Maps/TUB - CC 85 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
TUB.cc45.map[TUB.cc45.map < thresh] = 0
TUB.cc45.map[TUB.cc45.map >= thresh] = 1
TUB.cc85.map[TUB.cc85.map < thresh] = 0
TUB.cc85.map[TUB.cc85.map >= thresh] = 1

# VER
VER.cc45.clim <- subset(cc45.clim, VER.predictors)
VER.cc85.clim <- subset(cc85.clim, VER.predictors)
VER.cc45.map <- predict(VER, VER.cc45.clim,
                      filename = './Composite Maps/Species Maps/VER - CC 45 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
VER.cc85.map <- predict(VER, VER.cc85.clim,
                      filename = './Composite Maps/Species Maps/VER - CC 85 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
VER.cc45.map[VER.cc45.map < thresh] = 0
VER.cc45.map[VER.cc45.map >= thresh] = 1
VER.cc85.map[VER.cc85.map < thresh] = 0
VER.cc85.map[VER.cc85.map >= thresh] = 1

# VIR
VIR.cc45.clim <- subset(cc45.clim, VIR.predictors)
VIR.cc85.clim <- subset(cc85.clim, VIR.predictors)
VIR.cc45.map <- predict(VIR, VIR.cc45.clim,
                      filename = './Composite Maps/Species Maps/VIR - CC 45 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
VIR.cc85.map <- predict(VIR, VIR.cc85.clim,
                      filename = './Composite Maps/Species Maps/VIR - CC 85 Map',
                      format = 'GTiff',  progress = "text", overwrite = T)
VIR.cc45.map[VIR.cc45.map < thresh] = 0
VIR.cc45.map[VIR.cc45.map >= thresh] = 1
VIR.cc85.map[VIR.cc85.map < thresh] = 0
VIR.cc85.map[VIR.cc85.map >= thresh] = 1

# Make synthetic maps
future.cc45.map <- (ASP.cc45.map + ERO.cc45.map + FAS.cc45.map + INC.cc45.map + OEN.cc45.map +
                    SPE.cc45.map + SYR.cc45.map + TUB.cc45.map + VER.cc45.map + VIR.cc45.map)

jpeg(file = "./Composite Maps/Composite Map Pics/cc45_compmap.jpg", 
     width = 1080, height = 1138, quality = 100)
plot(future.cc45.map, main = "CC 4.5", col = rev(bpy.colors(7)))
sp::plot(countries, add = T, border = 'gray45')
sp::plot(states, add = T, border = 'gray45')
dev.off()

future.cc85.map <- (ASP.cc85.map + ERO.cc85.map + FAS.cc85.map + INC.cc85.map + OEN.cc85.map +
                    SPE.cc85.map + SYR.cc85.map + TUB.cc85.map + VER.cc85.map + VIR.cc85.map)

jpeg(file = "./Composite Maps/Composite Map Pics/cc85_compmap.jpg", 
     width = 1080, height = 1138, quality = 100)
plot(future.cc85.map, main = "CC 8.5", col = rev(bpy.colors(7)))
sp::plot(countries, add = T, border = 'gray45')
sp::plot(states, add = T, border = 'gray45')
dev.off()




