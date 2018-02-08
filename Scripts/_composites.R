##  ---------------------------------------------------------------------------------------------------------------------  ##
                           # Composite Maps for Functional Groups
##  ---------------------------------------------------------------------------------------------------------------------  ##
# Written by Nicholas J Lyon
  ## nicholasjlyon@gmail.com

# PURPOSE
  ## This script will generate composite maps to show full functional group response
  ## If areas of high co-occurrence (within functional group) are more abundant/obvious than
  ## areas of low co-occurence, then functional group is a good predictor of climate response
  ## I'll also plan on making 2 plots per future per functional group
  ## one showing areas of decreasing suitability and one showing increasing suitability
  ## It's more than possible that what is "good" might vary while what is "bad" might be universal
  ## or vice versa, but we'll cross that bridge when we come to it

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

# Handy shapefiles denoting country/state borders
countries <- rgdal::readOGR('./Borders', 'TM_WORLD_BORDERS-0.3')
states <- rgdal::readOGR('./Borders', 'cb_2015_us_state_20m')

##  ---------------------------------------------------------------------------------------------------------------------  ##
                                # Get the Climate Data Ready
##  ---------------------------------------------------------------------------------------------------------------------  ##
# Stack the current climate data into a single RasterStack
clim <- stack(c('./WORLDCLIM Data/Cropped Elevation/elevation.tif',
                    list.files('./WORLDCLIM Data/Cropped 1960-1990/', full.names = T)))

# Now do the same for both RCP 4.5 and 8.5 futures
clim.45 <- stack(c('./WORLDCLIM Data/Cropped Elevation/elevation.tif',
                   list.files('./WORLDCLIM Data/Cropped CC 45/', full.names = T)))

clim.85 <- stack(c('./WORLDCLIM Data/Cropped Elevation/elevation.tif',
                   list.files('./WORLDCLIM Data/Cropped CC 85/', full.names = T)))

##  ---------------------------------------------------------------------------------------------------------------------  ##
                          # Prep Maps for all Species Individually
##  ---------------------------------------------------------------------------------------------------------------------  ##
# Load in finished baseline model
elyvir.intermed <- load("./Preliminary Models/ELYVIR Model V3 - Baseline.Rdata")

# Get the model into a different object
elyvir <- baselineModel

# Subset the full climate data for just these predictor variables
elyvir.predclim <- subset(clim, elyvir.pred)

# We'll want it eventually, so let's go ahead and make the predicted map for each species too
elyvir.map <- predict(elyvir, elyvir.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/ELYVIR - Baseline Map',
                       format = 'GTiff', progress = "text", overwrite = T)

# Get a map that will show bad and a map that will show good!
elyvir.neg <- elyvir.map
elyvir.pos <- elyvir.map

# Want to convert the values in both maps to thresholds for later summation across species (w/in fxnl group)
elyvir.neg[elyvir.neg < thresh] = 1
elyvir.neg[elyvir.neg >= thresh] = 0
elyvir.pos[elyvir.pos < thresh] = 0
elyvir.pos[elyvir.pos >= thresh] = 1

# Once you have that for current conditions, let's go ahead and do it for both futures too!
  ## RCP 4.5
elyvir.predclim.45 <- subset(clim.45, elyvir.pred)
elyvir.45map <- predict(elyvir, elyvir.predclim.45, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/ELYVIR - CC 4.5 Map',
                      format = 'GTiff', progress = "text", overwrite = T)
elyvir.45neg <- elyvir.45map
elyvir.45pos <- elyvir.45map
elyvir.45neg[elyvir.45neg < thresh] = 1
elyvir.45neg[elyvir.45neg >= thresh] = 0
elyvir.45pos[elyvir.45pos < thresh] = 0
elyvir.45pos[elyvir.45pos >= thresh] = 1

  ## RCP 8.5
elyvir.predclim.85 <- subset(clim.85, elyvir.pred)
elyvir.85map <- predict(elyvir, elyvir.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/ELYVIR - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
elyvir.85neg <- elyvir.85map
elyvir.85pos <- elyvir.85map
elyvir.85neg[elyvir.85neg < thresh] = 1
elyvir.85neg[elyvir.85neg >= thresh] = 0
elyvir.85pos[elyvir.85pos < thresh] = 0
elyvir.85pos[elyvir.85pos >= thresh] = 1

# Now do that for every other modeled species (comments excluded for others for space efficiency)
  ## KOEMAC (Koeleria macrantha)
koemac.intermed <- load("./Preliminary Models/KOEMAC Model V3 - Baseline.Rdata")
koemac <- baselineModel
koemac.predclim <- subset(clim, koemac.pred)
koemac.map <- predict(koemac, koemac.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/KOEMAC - Baseline Map',
                      format = 'GTiff', progress = "text", overwrite = T)
koemac.neg <- koemac.map
koemac.pos <- koemac.map
koemac.neg[koemac.neg < thresh] = 1
koemac.neg[koemac.neg >= thresh] = 0
koemac.pos[koemac.pos < thresh] = 0
koemac.pos[koemac.pos >= thresh] = 1
koemac.predclim.45 <- subset(clim.45, koemac.pred)
koemac.45map <- predict(koemac, koemac.predclim.45, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/KOEMAC - CC 4.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
koemac.45neg <- koemac.45map
koemac.45pos <- koemac.45map
koemac.45neg[koemac.45neg < thresh] = 1
koemac.45neg[koemac.45neg >= thresh] = 0
koemac.45pos[koemac.45pos < thresh] = 0
koemac.45pos[koemac.45pos >= thresh] = 1
koemac.predclim.85 <- subset(clim.85, koemac.pred)
koemac.85map <- predict(koemac, koemac.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/KOEMAC - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
koemac.85neg <- koemac.85map
koemac.85pos <- koemac.85map
koemac.85neg[koemac.85neg < thresh] = 1
koemac.85neg[koemac.85neg >= thresh] = 0
koemac.85pos[koemac.85pos < thresh] = 0
koemac.85pos[koemac.85pos >= thresh] = 1

  ## STISPA (Stipa spartea)
stispa.intermed <- load("./Preliminary Models/STISPA Model V3 - Baseline.Rdata")
stispa <- baselineModel
stispa.predclim <- subset(clim, stispa.pred)
stispa.map <- predict(stispa, stispa.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/STISPA - Baseline Map',
                      format = 'GTiff', progress = "text", overwrite = T)
stispa.neg <- stispa.map
stispa.pos <- stispa.map
stispa.neg[stispa.neg < thresh] = 1
stispa.neg[stispa.neg >= thresh] = 0
stispa.pos[stispa.pos < thresh] = 0
stispa.pos[stispa.pos >= thresh] = 1
stispa.predclim.45 <- subset(clim.45, stispa.pred)
stispa.45map <- predict(stispa, stispa.predclim.45, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/STISPA - CC 4.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
stispa.45neg <- stispa.45map
stispa.45pos <- stispa.45map
stispa.45neg[stispa.45neg < thresh] = 1
stispa.45neg[stispa.45neg >= thresh] = 0
stispa.45pos[stispa.45pos < thresh] = 0
stispa.45pos[stispa.45pos >= thresh] = 1
stispa.predclim.85 <- subset(clim.85, stispa.pred)
stispa.85map <- predict(stispa, stispa.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/STISPA - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
stispa.85neg <- stispa.85map
stispa.85pos <- stispa.85map
stispa.85neg[stispa.85neg < thresh] = 1
stispa.85neg[stispa.85neg >= thresh] = 0
stispa.85pos[stispa.85pos < thresh] = 0
stispa.85pos[stispa.85pos >= thresh] = 1

  ## BOUCUR (Bouteloua curtipendula)
boucur.intermed <- load("./Preliminary Models/BOUCUR Model V3 - Baseline.Rdata")
boucur <- baselineModel
boucur.predclim <- subset(clim, boucur.pred)
boucur.map <- predict(boucur, boucur.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/BOUCUR - Baseline Map',
                      format = 'GTiff', progress = "text", overwrite = T)
boucur.neg <- boucur.map
boucur.pos <- boucur.map
boucur.neg[boucur.neg < thresh] = 1
boucur.neg[boucur.neg >= thresh] = 0
boucur.pos[boucur.pos < thresh] = 0
boucur.pos[boucur.pos >= thresh] = 1
boucur.predclim.45 <- subset(clim.45, boucur.pred)
boucur.45map <- predict(boucur, boucur.predclim.45, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/BOUCUR - CC 4.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
boucur.45neg <- boucur.45map
boucur.45pos <- boucur.45map
boucur.45neg[boucur.45neg < thresh] = 1
boucur.45neg[boucur.45neg >= thresh] = 0
boucur.45pos[boucur.45pos < thresh] = 0
boucur.45pos[boucur.45pos >= thresh] = 1
boucur.predclim.85 <- subset(clim.85, boucur.pred)
boucur.85map <- predict(boucur, boucur.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/BOUCUR - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
boucur.85neg <- boucur.85map
boucur.85pos <- boucur.85map
boucur.85neg[boucur.85neg < thresh] = 1
boucur.85neg[boucur.85neg >= thresh] = 0
boucur.85pos[boucur.85pos < thresh] = 0
boucur.85pos[boucur.85pos >= thresh] = 1

  ## SCHSCO (Schizachyrium scoparium)
schsco.intermed <- load("./Preliminary Models/SCHSCO Model V3 - Baseline.Rdata")
schsco <- baselineModel
schsco.predclim <- subset(clim, schsco.pred)
schsco.map <- predict(schsco, schsco.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/SCHSCO - Baseline Map',
                      format = 'GTiff', progress = "text", overwrite = T)
schsco.neg <- schsco.map
schsco.pos <- schsco.map
schsco.neg[schsco.neg < thresh] = 1
schsco.neg[schsco.neg >= thresh] = 0
schsco.pos[schsco.pos < thresh] = 0
schsco.pos[schsco.pos >= thresh] = 1
schsco.predclim.45 <- subset(clim.45, schsco.pred)
schsco.45map <- predict(schsco, schsco.predclim.45, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/SCHSCO - CC 4.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
schsco.45neg <- schsco.45map
schsco.45pos <- schsco.45map
schsco.45neg[schsco.45neg < thresh] = 1
schsco.45neg[schsco.45neg >= thresh] = 0
schsco.45pos[schsco.45pos < thresh] = 0
schsco.45pos[schsco.45pos >= thresh] = 1
schsco.predclim.85 <- subset(clim.85, schsco.pred)
schsco.85map <- predict(schsco, schsco.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/SCHSCO - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
schsco.85neg <- schsco.85map
schsco.85pos <- schsco.85map
schsco.85neg[schsco.85neg < thresh] = 1
schsco.85neg[schsco.85neg >= thresh] = 0
schsco.85pos[schsco.85pos < thresh] = 0
schsco.85pos[schsco.85pos >= thresh] = 1

  ## SORNUT (Sorghastrum nutans)
sornut.intermed <- load("./Preliminary Models/SORNUT Model V3 - Baseline.Rdata")
sornut <- baselineModel
sornut.predclim <- subset(clim, sornut.pred)
sornut.map <- predict(sornut, sornut.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/SORNUT - Baseline Map',
                      format = 'GTiff', progress = "text", overwrite = T)
sornut.neg <- sornut.map
sornut.pos <- sornut.map
sornut.neg[sornut.neg < thresh] = 1
sornut.neg[sornut.neg >= thresh] = 0
sornut.pos[sornut.pos < thresh] = 0
sornut.pos[sornut.pos >= thresh] = 1
sornut.predclim.45 <- subset(clim.45, sornut.pred)
sornut.45map <- predict(sornut, sornut.predclim.45, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/SORNUT - CC 4.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
sornut.45neg <- sornut.45map
sornut.45pos <- sornut.45map
sornut.45neg[sornut.45neg < thresh] = 1
sornut.45neg[sornut.45neg >= thresh] = 0
sornut.45pos[sornut.45pos < thresh] = 0
sornut.45pos[sornut.45pos >= thresh] = 1
sornut.predclim.85 <- subset(clim.85, sornut.pred)
sornut.85map <- predict(sornut, sornut.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/SORNUT - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
sornut.85neg <- sornut.85map
sornut.85pos <- sornut.85map
sornut.85neg[sornut.85neg < thresh] = 1
sornut.85neg[sornut.85neg >= thresh] = 0
sornut.85pos[sornut.85pos < thresh] = 0
sornut.85pos[sornut.85pos >= thresh] = 1

  ## ASCINC (Asclepias incarnata)
ascinc.intermed <- load("./Preliminary Models/ASCINC Model V3 - Baseline.Rdata")
ascinc <- baselineModel
ascinc.predclim <- subset(clim, ascinc.pred)
ascinc.map <- predict(ascinc, ascinc.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/ASCINC - Baseline Map',
                      format = 'GTiff', progress = "text", overwrite = T)
ascinc.neg <- ascinc.map
ascinc.pos <- ascinc.map
ascinc.neg[ascinc.neg < thresh] = 1
ascinc.neg[ascinc.neg >= thresh] = 0
ascinc.pos[ascinc.pos < thresh] = 0
ascinc.pos[ascinc.pos >= thresh] = 1
ascinc.predclim.45 <- subset(clim.45, ascinc.pred)
ascinc.45map <- predict(ascinc, ascinc.predclim.45, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/ASCINC - CC 4.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
ascinc.45neg <- ascinc.45map
ascinc.45pos <- ascinc.45map
ascinc.45neg[ascinc.45neg < thresh] = 1
ascinc.45neg[ascinc.45neg >= thresh] = 0
ascinc.45pos[ascinc.45pos < thresh] = 0
ascinc.45pos[ascinc.45pos >= thresh] = 1
ascinc.predclim.85 <- subset(clim.85, ascinc.pred)
ascinc.85map <- predict(ascinc, ascinc.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/ASCINC - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
ascinc.85neg <- ascinc.85map
ascinc.85pos <- ascinc.85map
ascinc.85neg[ascinc.85neg < thresh] = 1
ascinc.85neg[ascinc.85neg >= thresh] = 0
ascinc.85pos[ascinc.85pos < thresh] = 0
ascinc.85pos[ascinc.85pos >= thresh] = 1

  ## ASCSYR (Asclepias syriaca)
ascsyr.intermed <- load("./Preliminary Models/ASCSYR Model V3 - Baseline.Rdata")
ascsyr <- baselineModel
ascsyr.predclim <- subset(clim, ascsyr.pred)
ascsyr.map <- predict(ascsyr, ascsyr.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/ASCSYR - Baseline Map',
                      format = 'GTiff', progress = "text", overwrite = T)
ascsyr.neg <- ascsyr.map
ascsyr.pos <- ascsyr.map
ascsyr.neg[ascsyr.neg < thresh] = 1
ascsyr.neg[ascsyr.neg >= thresh] = 0
ascsyr.pos[ascsyr.pos < thresh] = 0
ascsyr.pos[ascsyr.pos >= thresh] = 1
ascsyr.predclim.45 <- subset(clim.45, ascsyr.pred)
ascsyr.45map <- predict(ascsyr, ascsyr.predclim.45, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/ASCSYR - CC 4.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
ascsyr.45neg <- ascsyr.45map
ascsyr.45pos <- ascsyr.45map
ascsyr.45neg[ascsyr.45neg < thresh] = 1
ascsyr.45neg[ascsyr.45neg >= thresh] = 0
ascsyr.45pos[ascsyr.45pos < thresh] = 0
ascsyr.45pos[ascsyr.45pos >= thresh] = 1
ascsyr.predclim.85 <- subset(clim.85, ascsyr.pred)
ascsyr.85map <- predict(ascsyr, ascsyr.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/ASCSYR - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
ascsyr.85neg <- ascsyr.85map
ascsyr.85pos <- ascsyr.85map
ascsyr.85neg[ascsyr.85neg < thresh] = 1
ascsyr.85neg[ascsyr.85neg >= thresh] = 0
ascsyr.85pos[ascsyr.85pos < thresh] = 0
ascsyr.85pos[ascsyr.85pos >= thresh] = 1

  ## ASCTUB (Asclepias tuberosa)
asctub.intermed <- load("./Preliminary Models/ASCTUB Model V3 - Baseline.Rdata")
asctub <- baselineModel
asctub.predclim <- subset(clim, asctub.pred)
asctub.map <- predict(asctub, asctub.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/ASCTUB - Baseline Map',
                      format = 'GTiff', progress = "text", overwrite = T)
asctub.neg <- asctub.map
asctub.pos <- asctub.map
asctub.neg[asctub.neg < thresh] = 1
asctub.neg[asctub.neg >= thresh] = 0
asctub.pos[asctub.pos < thresh] = 0
asctub.pos[asctub.pos >= thresh] = 1
asctub.predclim.45 <- subset(clim.45, asctub.pred)
asctub.45map <- predict(asctub, asctub.predclim.45, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/ASCTUB - CC 4.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
asctub.45neg <- asctub.45map
asctub.45pos <- asctub.45map
asctub.45neg[asctub.45neg < thresh] = 1
asctub.45neg[asctub.45neg >= thresh] = 0
asctub.45pos[asctub.45pos < thresh] = 0
asctub.45pos[asctub.45pos >= thresh] = 1
asctub.predclim.85 <- subset(clim.85, asctub.pred)
asctub.85map <- predict(asctub, asctub.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/ASCTUB - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
asctub.85neg <- asctub.85map
asctub.85pos <- asctub.85map
asctub.85neg[asctub.85neg < thresh] = 1
asctub.85neg[asctub.85neg >= thresh] = 0
asctub.85pos[asctub.85pos < thresh] = 0
asctub.85pos[asctub.85pos >= thresh] = 1

  ## DRYARG (Drymocallis arguta)
dryarg.intermed <- load("./Preliminary Models/DRYARG Model V3 - Baseline.Rdata")
dryarg <- baselineModel
dryarg.predclim <- subset(clim, dryarg.pred)
dryarg.map <- predict(dryarg, dryarg.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/DRYARG - Baseline Map',
                      format = 'GTiff', progress = "text", overwrite = T)
dryarg.neg <- dryarg.map
dryarg.pos <- dryarg.map
dryarg.neg[dryarg.neg < thresh] = 1
dryarg.neg[dryarg.neg >= thresh] = 0
dryarg.pos[dryarg.pos < thresh] = 0
dryarg.pos[dryarg.pos >= thresh] = 1
dryarg.predclim.45 <- subset(clim.45, dryarg.pred)
dryarg.45map <- predict(dryarg, dryarg.predclim.45, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/DRYARG - CC 4.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
dryarg.45neg <- dryarg.45map
dryarg.45pos <- dryarg.45map
dryarg.45neg[dryarg.45neg < thresh] = 1
dryarg.45neg[dryarg.45neg >= thresh] = 0
dryarg.45pos[dryarg.45pos < thresh] = 0
dryarg.45pos[dryarg.45pos >= thresh] = 1
dryarg.predclim.85 <- subset(clim.85, dryarg.pred)
dryarg.85map <- predict(dryarg, dryarg.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/DRYARG - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
dryarg.85neg <- dryarg.85map
dryarg.85pos <- dryarg.85map
dryarg.85neg[dryarg.85neg < thresh] = 1
dryarg.85neg[dryarg.85neg >= thresh] = 0
dryarg.85pos[dryarg.85pos < thresh] = 0
dryarg.85pos[dryarg.85pos >= thresh] = 1

  ## LOBSIP (Lobelia siphilitica)
lobsip.intermed <- load("./Preliminary Models/LOBSIP Model V3 - Baseline.Rdata")
lobsip <- baselineModel
lobsip.predclim <- subset(clim, lobsip.pred)
lobsip.map <- predict(lobsip, lobsip.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/LOBSIP - Baseline Map',
                      format = 'GTiff', progress = "text", overwrite = T)
lobsip.neg <- lobsip.map
lobsip.pos <- lobsip.map
lobsip.neg[lobsip.neg < thresh] = 1
lobsip.neg[lobsip.neg >= thresh] = 0
lobsip.pos[lobsip.pos < thresh] = 0
lobsip.pos[lobsip.pos >= thresh] = 1
lobsip.predclim.45 <- subset(clim.45, lobsip.pred)
lobsip.45map <- predict(lobsip, lobsip.predclim.45, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/LOBSIP - CC 4.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
lobsip.45neg <- lobsip.45map
lobsip.45pos <- lobsip.45map
lobsip.45neg[lobsip.45neg < thresh] = 1
lobsip.45neg[lobsip.45neg >= thresh] = 0
lobsip.45pos[lobsip.45pos < thresh] = 0
lobsip.45pos[lobsip.45pos >= thresh] = 1
lobsip.predclim.85 <- subset(clim.85, lobsip.pred)
lobsip.85map <- predict(lobsip, lobsip.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/LOBSIP - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
lobsip.85neg <- lobsip.85map
lobsip.85pos <- lobsip.85map
lobsip.85neg[lobsip.85neg < thresh] = 1
lobsip.85neg[lobsip.85neg >= thresh] = 0
lobsip.85pos[lobsip.85pos < thresh] = 0
lobsip.85pos[lobsip.85pos >= thresh] = 1

  ## MONFIS (Monarda fistulosa)
monfis.intermed <- load("./Preliminary Models/MONFIS Model V3 - Baseline.Rdata")
monfis <- baselineModel
monfis.predclim <- subset(clim, monfis.pred)
monfis.map <- predict(monfis, monfis.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/MONFIS - Baseline Map',
                      format = 'GTiff', progress = "text", overwrite = T)
monfis.neg <- monfis.map
monfis.pos <- monfis.map
monfis.neg[monfis.neg < thresh] = 1
monfis.neg[monfis.neg >= thresh] = 0
monfis.pos[monfis.pos < thresh] = 0
monfis.pos[monfis.pos >= thresh] = 1
monfis.predclim.45 <- subset(clim.45, monfis.pred)
monfis.45map <- predict(monfis, monfis.predclim.45, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/MONFIS - CC 4.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
monfis.45neg <- monfis.45map
monfis.45pos <- monfis.45map
monfis.45neg[monfis.45neg < thresh] = 1
monfis.45neg[monfis.45neg >= thresh] = 0
monfis.45pos[monfis.45pos < thresh] = 0
monfis.45pos[monfis.45pos >= thresh] = 1
monfis.predclim.85 <- subset(clim.85, monfis.pred)
monfis.85map <- predict(monfis, monfis.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/MONFIS - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
monfis.85neg <- monfis.85map
monfis.85pos <- monfis.85map
monfis.85neg[monfis.85neg < thresh] = 1
monfis.85neg[monfis.85neg >= thresh] = 0
monfis.85pos[monfis.85pos < thresh] = 0
monfis.85pos[monfis.85pos >= thresh] = 1

  ## AMOCAN (Amorpha canscens)
amocan.intermed <- load("./Preliminary Models/AMOCAN Model V3 - Baseline.Rdata")
amocan <- baselineModel
amocan.predclim <- subset(clim, amocan.pred)
amocan.map <- predict(amocan, amocan.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/AMOCAN - Baseline Map',
                      format = 'GTiff', progress = "text", overwrite = T)
amocan.neg <- amocan.map
amocan.pos <- amocan.map
amocan.neg[amocan.neg < thresh] = 1
amocan.neg[amocan.neg >= thresh] = 0
amocan.pos[amocan.pos < thresh] = 0
amocan.pos[amocan.pos >= thresh] = 1
amocan.predclim.45 <- subset(clim.45, amocan.pred)
amocan.45map <- predict(amocan, amocan.predclim.45, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/AMOCAN - CC 4.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
amocan.45neg <- amocan.45map
amocan.45pos <- amocan.45map
amocan.45neg[amocan.45neg < thresh] = 1
amocan.45neg[amocan.45neg >= thresh] = 0
amocan.45pos[amocan.45pos < thresh] = 0
amocan.45pos[amocan.45pos >= thresh] = 1
amocan.predclim.85 <- subset(clim.85, amocan.pred)
amocan.85map <- predict(amocan, amocan.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/AMOCAN - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
amocan.85neg <- amocan.85map
amocan.85pos <- amocan.85map
amocan.85neg[amocan.85neg < thresh] = 1
amocan.85neg[amocan.85neg >= thresh] = 0
amocan.85pos[amocan.85pos < thresh] = 0
amocan.85pos[amocan.85pos >= thresh] = 1

  ## PETCAN (Petalostemon candida)
petcan.intermed <- load("./Preliminary Models/PETCAN Model V3 - Baseline.Rdata")
petcan <- baselineModel
petcan.predclim <- subset(clim, petcan.pred)
petcan.map <- predict(petcan, petcan.predclim, 
                      filename = './Maps/Composite Maps/Necessary Pre-Products/PETCAN - Baseline Map',
                      format = 'GTiff', progress = "text", overwrite = T)
petcan.neg <- petcan.map
petcan.pos <- petcan.map
petcan.neg[petcan.neg < thresh] = 1
petcan.neg[petcan.neg >= thresh] = 0
petcan.pos[petcan.pos < thresh] = 0
petcan.pos[petcan.pos >= thresh] = 1
petcan.predclim.45 <- subset(clim.45, petcan.pred)
petcan.45map <- predict(petcan, petcan.predclim.45, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/PETCAN - CC 4.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
petcan.45neg <- petcan.45map
petcan.45pos <- petcan.45map
petcan.45neg[petcan.45neg < thresh] = 1
petcan.45neg[petcan.45neg >= thresh] = 0
petcan.45pos[petcan.45pos < thresh] = 0
petcan.45pos[petcan.45pos >= thresh] = 1
petcan.predclim.85 <- subset(clim.85, petcan.pred)
petcan.85map <- predict(petcan, petcan.predclim.85, 
                        filename = './Maps/Composite Maps/Necessary Pre-Products/PETCAN - CC 8.5 Map',
                        format = 'GTiff', progress = "text", overwrite = T)
petcan.85neg <- petcan.85map
petcan.85pos <- petcan.85map
petcan.85neg[petcan.85neg < thresh] = 1
petcan.85neg[petcan.85neg >= thresh] = 0
petcan.85pos[petcan.85pos < thresh] = 0
petcan.85pos[petcan.85pos >= thresh] = 1

##  ---------------------------------------------------------------------------------------------------------------------  ##
                                    # Construct Composite Maps!
##  ---------------------------------------------------------------------------------------------------------------------  ##
# Make a composite map for the positive and negative versions of present conditions
comp.pos <- (elyvir.pos + koemac.pos + stispa.pos + 
               boucur.pos + schsco.pos + sornut.pos +
               ascinc.pos + ascinc.pos + ascinc.pos +
               dryarg.pos + lobsip.pos + monfis.pos + 
               amocan.pos + petcan.pos)

comp.neg <- (elyvir.neg + koemac.neg + stispa.neg + 
               boucur.neg + schsco.neg + sornut.neg +
               ascinc.neg + ascinc.neg + ascinc.neg +
               dryarg.neg + lobsip.neg + monfis.neg + 
               amocan.neg + petcan.neg)

# Do it for the future conditions too!
  ## RCP 4.5
comp.45pos <- (elyvir.45pos + koemac.45pos + stispa.45pos + 
               boucur.45pos + schsco.45pos + sornut.45pos +
               ascinc.45pos + ascinc.45pos + ascinc.45pos +
               dryarg.45pos + lobsip.45pos + monfis.45pos + 
               amocan.45pos + petcan.45pos)

comp.45neg <- (elyvir.45neg + koemac.45neg + stispa.45neg + 
               boucur.45neg + schsco.45neg + sornut.45neg +
               ascinc.45neg + ascinc.45neg + ascinc.45neg +
               dryarg.45neg + lobsip.45neg + monfis.45neg + 
               amocan.45neg + petcan.45neg)

  ## RCP 8.5
comp.85pos <- (elyvir.85pos + koemac.85pos + stispa.85pos + 
                 boucur.85pos + schsco.85pos + sornut.85pos +
                 ascinc.85pos + ascinc.85pos + ascinc.85pos +
                 dryarg.85pos + lobsip.85pos + monfis.85pos + 
                 amocan.85pos + petcan.85pos)

comp.85neg <- (elyvir.85neg + koemac.85neg + stispa.85neg + 
                 boucur.85neg + schsco.85neg + sornut.85neg +
                 ascinc.85neg + ascinc.85neg + ascinc.85neg +
                 dryarg.85neg + lobsip.85neg + monfis.85neg + 
                 amocan.85neg + petcan.85neg)

# Plot them!
jpeg(file = "./Maps/Composite Maps/comp_positive.jpg", 
     width = 1080, height = 1138, quality = 100)
plot(present.map, col = rev(topo.colors(7)))
sp::plot(countries, add = T, border = 'gray45')
sp::plot(states, add = T, border = 'gray45')
dev.off()



