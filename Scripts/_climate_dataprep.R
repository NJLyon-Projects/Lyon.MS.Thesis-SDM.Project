##  ---------------------------------------------------------------------------------------------------------------------  ##
                                  # Climate Dataprep Script
##  ---------------------------------------------------------------------------------------------------------------------  ##
# Written by Nicholas J Lyon
  ## nicholasjlyon@gmail.com

# SCRIPT PURPOSE
  ## To generate these models one must first obtain and crop climate data through the chosen extent
  ## This is computationally time-consuming, so this script was created to offer a short cut
  ## By getting an aggregate pull of all modeled species, climate can be obtained over the shared extent
  ## of all the species, thus necessitating only a single data acquisition step and climate raster cropping step
  ## Per future condition of course, but this is still a dramatic improvement over the species-specific way
  ## of doing things. Cheers!

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

# You'll need a shapefile of ecotones if you want to define the study extent for later climate cropping
eco <- rgdal::readOGR('./Borders', 'provinces', verbose = F)

##  ---------------------------------------------------------------------------------------------------------------------  ##
                           # Aggregate Occurrence Record Pull
##  ---------------------------------------------------------------------------------------------------------------------  ##
# Get a vector of all modeled species' synonyms for easy later reference
spnames <- c('Elymus virginicus', 'Elymus carolinianus', 'Hordeum virginicum', 
             'Elymus virginicanus', 'Elymus striatus', 'Koeleria macrantha', 'Koeleria gracilis', 
             'Koeleria albescens', 'Koeleria cristata', 'Stipa spartea', 'Hesperostipa spartea',
             'Stipa robusta', 'Hesperostipa spartea', 'Bouteloua curtipendula', 'Bouteloua racemosa',
             'Atheropogon apludoides', 'Schizachyrium scoparium', 'Andropogon spadiceus',
             'Andropogon divergens', 'Sorghastrum nutans', 'Andropogon nutans', 'Sorghastrum avenaceum', 
             'Asclepias incarnata', 'Asclepias syriaca', 'Asclepias tuberosa', 'Drymocallis arguta',
             'Bootia sericea', 'Drymocallis agrimonioides', 'Potentilla arguta', 'Lobelia siphilitica',
             'Lobelia reflexa', 'Lobelia antisyphilitica', 'Rapuntium siphiliticum', 'Monarda fistulosa',
             'Monarda menthifolia', 'Amorpha canescens', 'Amorpha brachycarpa', 'Dalea candida',
             'Psoralea candida', 'Kuhnistera candida', 'Petalostemon candidum')

# Actually get these occurrence records
  ## This'll take awhile, so be sure this says exactly what you want it to before you run it
df_mult <- occ(query = spnames, from = c('gbif', 'ecoengine', 'bison'), limit = 10000,
               geometry = c(-140, 22, -58, 55), has_coords = T)

# "Fix" the species names to match the entries provided in the "spnames" vector
df_nam <- fixnames(df_mult, how = "query")

# And get the tibbles into a single dataframe
df_comb <- occ2df(df_nam)

# Get rid of records with impossible geo-references or ones where the reference is missing
records_geoclean <- df_comb %>% coord_impossible(drop = T) %>% coord_incomplete(drop = T) %>% coord_unlikely(drop = T)

# Get rid of records that are not appropriately date referenced
records_dateclean <- records_geoclean
records_dateclean <- date_standardize(records_dateclean, "%Y")
records_dateclean <- date_missing(records_dateclean, format = "%Y", drop = T)
records_dateclean$date <- as.numeric(records_dateclean$date)

# And ditch records that are from outside of the training climate range
records <- subset(records_dateclean, records_dateclean$date >= 1960 & records_dateclean$date <= 1990)

# How many does that leave us with?
nrow(records)

# Get a spatial dataframe for climate cropping
recordsSpatial <- SpatialPointsDataFrame(coords = cbind(records$longitude, records$latitude), data = records,
                                proj4string = CRS('+proj=longlat +datumWGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

# Save out these species records
saveRDS(records, './Species Records/All Spp - Aggregate Pull.rds')

##  ----------------------------------------------------------------------------------------------------------------------  ##
                                     # Define Study Extent
##  ----------------------------------------------------------------------------------------------------------------------  ##
# Ditch ocean data
eco <- eco[-which(eco$DOM_DESC == 'outside polygon'), ]

# Get that spatial dataframe
recordsSpatial <- SpatialPointsDataFrame(coords = cbind(records$longitude, records$latitude), data = records,
                          proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

# Get a touch matrix for the ecoregions where either an occurrence record was or adjacent to one that had one
ecoContain <- eco[recordsSpatial, ]
touchMatrix <- gTouches(eco, ecoContain, byid = T)
touchVector <- colSums(touchMatrix)
ecoStudyRegion <- eco[touchVector > 0, ]
ecoStudyRegion <- rbind(ecoStudyRegion, ecoContain, makeUniqueIDs = T)
studyExtent <- extent(ecoStudyRegion)

save(studyExtent, file = './WORLDCLIM Data/All Spp - Aggregate Study Region.Rdata', compress = T)

##  ----------------------------------------------------------------------------------------------------------------------  ##
                               # Clipping Current Climate Rasters
##  ----------------------------------------------------------------------------------------------------------------------  ##
# Get the global (i.e. uncropped) climate data
elevation <- raster('./WORLDCLIM Data/Global Elevation/alt.bil')

# Crop it to the aggregate study extent & set the min/max
elevation <- crop(elevation, studyExtent)
elevation <- setMinMax(elevation)

# Now save that out for later use
projection(elevation) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
writeRaster(elevation, './WORLDCLIM Data/Cropped Elevation/elevation',
            format = 'GTiff', datatype = 'INT2S', overwrite = T)

# Do the same iteratively for all the BIOCLIM data that you are interested in
for (i in 1:19) {
  print(paste('Clipping current WORLDCLIM raster ', i))
  flush.console()
  
  rast <- raster(paste0('./WORLDCLIM Data/Global 1960-1990/bio', i, '.bil'))
  rast <- crop(rast, studyExtent)
  rast <- setMinMax(rast)
  projection(rast) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  
  writeRaster(rast, paste0('./WORLDCLIM Data/Cropped 1960-1990/WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  }

##  ---------------------------------------------------------------------------------------------------------------------  ##
                                # Clipping Future Climate Rasters
##  ---------------------------------------------------------------------------------------------------------------------  ##
# CCSM4 RCP4.5 @ 2050 (Moderate Hot + No Precip ∆)
for (i in 1:19) {
  print(paste('Clipping CCSM4 4.5 future WORLDCLIM raster ', i))
  flush.console()
  
  future1 <- raster(paste0('./WORLDCLIM Data/Global CC 45/cc45bi50', i, '.tif'))
  future1 <- crop(future1, studyExtent)
  future1 <- setMinMax(future1)
  projection(future1) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  
  writeRaster(future1, paste0('./WORLDCLIM Data/Cropped CC 45/WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
  }

# CCSM4 RCP8.5 @ 2050 (Hotter with relatively drier conditions
for (i in 1:19) {
  print(paste('Clipping CCSM4 8.5 future WORLDCLIM raster ', i))
  flush.console()
  
  future2 <- raster(paste0('./WORLDCLIM Data/Global CC 85/cc85bi50', i, '.tif'))
  future2 <- crop(future2, studyExtent)
  future2 <- setMinMax(future2)
  projection(future2) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  
  writeRaster(future2, paste0('./WORLDCLIM Data/Cropped CC 85/WC', ifelse(i < 10, '0', ''), i),
              format = 'GTiff', datatype = 'INT2S', overwrite = T)
}




