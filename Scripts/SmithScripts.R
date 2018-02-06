# Functions written by Adam Smith
# lab webpage is here:
  ## http://www.earthskysea.org/

elimCellDups <- function(pointsFrame, baseRaster, longLatFields=c('LONG_WGS84', 'LAT_WGS84'), priority=NULL) {
  # elimCellDups For a set of geographic points, this script eliminates all but one point in each cell. An optional priority vector can be passed to the function so that points with higher priority are kept.
  #
  # ARGUMENTS
  #
  # pointsFrame
  # data frame of points with at least the fields named in longLatFields
  #
  # baseRaster
  # Raster* object against which to compare points data frame
  #
  # longLatFields
  # character list of the names of the variables that have longitude and latitude... longitude must be listed first
  #
  # priority
  # optional, numeric or character vector same length as pointsFrame has rows, with values corresponding to priority in which records in the same cell will be kept (though only one point per cell will be kept)... these MUST be specified in alphanumeric order, so that, for example, a record with a priority value of '1' will be preferred over one with a value of '2', which will be preferred over a value of 'alpha' which is preferred over 'beta'
  #
  # VALUES
  #
  # cleanedData
  # data frame with just one point per cell in baseRaster
  #
  # BAUHAUS
  # - for all points get cell numbers in which they fall in baseRaster
  # - adjoin pointsFrame and cell number vector
  # - sort by priority if extant
  # - eliminate all but one point with duplicate cell numbers
  #
  # EXAMPLE
  # theseRecords <- elimCellDups(pointsFrame=data.frame(LAT_WGS84=c(38, 38.001, 39, 40, 40.01, 38), LONG_WGS84=c(-90, -90, -90, -90, -90, -90), Extra=letters[1:6] ), longLatFields=c('LONG_WGS84', 'LAT_WGS84'), baseRaster=raster('C:/ecology/DropBox/r/areaRasters/areaRaster_easternUsgsProvinces_terrestrialCellsEqual1.tif'), Priority=c(2, 2, 2, 2, 2, 1) )
  #
  # SOURCE	source('C:/ecology/Drive/Workshops/SDM from Start to Finish (KSU, 2016-02)/Scripts/Eliminate Points in Same Cell of a Raster.r')
  #
  # LICENSE
  # This document is copyright ?2011 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/
  #
  # AUTHOR	Adam B. Smith | Missouri Botanical Garden | adamDOTsmithATmobotDOTorg
  # DATE		2012-01 Adapted from "Utility - Eliminate Species Records in Same Cell of a Raster.r"
  # REVISIONS 2012-01 Now returns empty data frame if input frame had no rows
  
  #############################
  ## libraries and functions ##
  #############################
  
  require(raster)
  
  ##########################
  ## initialize variables ##
  ##########################
  
  ##########
  ## MAIN ##
  ##########
  
  if (nrow(pointsFrame) > 0) { # frame has at least one row
    
    # get cell numbers for each point and adjoin with data frame
    pointsFrame$cellNoTEMP <- cellFromXY(object=baseRaster, xy=cbind(pointsFrame[ , longLatFields[1] ], pointsFrame[ , longLatFields[2] ]) )
    
    # remember original row names
    pointsFrame$origRowNamesTEMP <- rownames(pointsFrame)
    
    # sort by priority
    if (!is.null(priority)) pointsFrame <- pointsFrame[ order(priority), ]
    
    # sort by cell number
    pointsFrame <- pointsFrame[ order(pointsFrame$cellNoTEMP), ]
    
    # reassign rownames so that the "unique" function keeps top-most row
    rownames(pointsFrame) <- 1:nrow(pointsFrame)
    
    # get top-most point in each cell
    cleanFrame <- pointsFrame[ rownames( unique( as.data.frame(pointsFrame$cellNoTEMP) ) ), ]
    
    # re-assign original row names
    rownames(cleanFrame) <- cleanFrame$origRowNamesTEMP
    
    # remove column with original row names
    cleanFrame$origRowNamesTEMP <- NULL
    
    # remove cell number column
    cleanFrame$cellNoTEMP <- NULL
    
    # original order
    cleanFrame <- cleanFrame[ order( rownames(cleanFrame) ), ]
    
  } else { # frame was empty!
    
    cleanFrame <- data.frame()
    
  }
  
  return(cleanFrame)
  
}



biovarsArbitrary <- function(
  targetStack,
  refStack,
  targetFunctSummary,
  refFunctSummary,
  refFunctCriterion,
  offset,
  targetDuration,
  refDuration,
  roundTarget=FALSE,
  outDir,
  outName,
  verbose=TRUE,
  ...
) {
  # biovarsArbitrary Calculates a BIOCLIM or BIOCLIM-like raster.  For example, can be used to calculate the mean of the mean temperature of the three hottest consecutive months (BIO10) or the mean of the minimum temperature of the two months after the three coldest months (a new variable).
  #
  # ARGUMENTS
  # targetStack
  # raster stack with 12 rasters, one per month, of the variable for which the new raster if created
  # 
  # refStack
  # raster stack with 12 rasters, one per month, which is used as a reference for calculating the target raster... e.g., if the mean of the maximum temperature of the wettest 6 months is desired, then targetStack is maximum temperature and refStack is precipitation
  #
  # targetFunctSummary
  # function used to aggregate across months of the target stack... can be a base function (e.g., mean, max, etc.) or a user-defined function...  should be able to handle a vector and produce a single output value... do not put this in quotes... ignored if targetDuration = 1
  #
  # refFunctSummary
  # functions used to summarize values across months in the reference stack... e.g., if wanting minimum of minimum temperature (targetStack) of the three months following the hottest three months (refStack), then this would be "mean" because the temperature of the reference stack is averaged across three month intervals... can be a base function (e.g., mean, max, etc.) or a user-defined function...  should be able to handle a vector and produce a single output value... do not put this in quotes... ignored if refDuration = 1
  #
  # refFunctCriterion
  # functions used to identify a month of period of continuois months in the reference stack... e.g., if wanting minimum of minimum temperature (targetStack) of the three months following the hottest three months (refStack), then this would be "max" because the maximum of the the temperature of the reference stack quarters is desired... can be a base function (e.g., mean, max, etc.) or a user-defined function...  should be able to handle a vector and produce a single output value... do not put this in quotes... ignored if refDuration = 1
  #
  # offset
  # number of months *after* (positive values) or *before* (negative values) the *first month* of the specified window of the reference time frame used to calculate values in the target stack... e.g., if wanting the mean of the maximum temperature (targetStack) of the three months following the three coldest months (refStack), then this is 3 (the first month for which values are to be calculated equals the third month after the start of the three coldest months)... if wanting the maximum temperature of the three months *before* the three coldest months, then this would be -3 (starts 3 months before the starting month of the reference period), and if wanting the maximum temperature of the month that starts on the first month of the three coldest months, this would be 0
  #
  # targetDuration
  # number of months of the target variable across which targetFunctSummary is calculated (must be >=1)
  #
  # refDuration
  # number of months of the reference variable across which refFunctSummary is calcualted (must be >=1)
  #
  # roundTarget
  # logical, will round output values if TRUE... good for reducing file size if extra precision is not needed or unwarranted
  #
  # outDir, outName
  # file path and name of file to which to save the raster
  #
  # verbose
  # logical, indicates progress
  #
  # ...
  # arguments to pass to writeRaster... see expecially format, overwrite, and dataType
  #
  # VALUES
  # raster with arbitrary BIOCLIM/WORLDCLIM values
  #
  # BAUHAUS
  #
  # EXAMPLE
  #
  # SOURCE	source('C:/ecology/Drive/R/Geography/Rasters - Calculate Arbitrary BIOCLIM Variable.r')
  #
  # LICENSE
  # This document is copyright ?2011 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/
  #
  # AUTHOR	Adam B. Smith | Missouri Botanical Garden, St. Louis, Missouri | adamDOTsmithATmobotDOTorg
  # DATE		2011-10 (forked from generic BIOCLIM function)
  # REVISIONS 2015-03-05 Can now handle rows of stacks that are completely NA (all missing values)
  
  if (verbose) cat('Calculating arbitrary BIOCLIM raster...\n'); flush.console()
  
  #############################
  ## libraries and functions ##
  #############################
  
  ##########################
  ## initialize variables ##
  ##########################
  
  targetFunctSummary <- match.fun(targetFunctSummary)
  refFunctSummary <- match.fun(refFunctSummary)
  refFunctCriterion <- match.fun(refFunctCriterion)
  
  
  ####################
  ## pre-processing ##
  ####################
  
  # make directory for output raster
  dir.create(path=outDir, recursive=TRUE, showWarnings=FALSE)
  
  # initiate raster
  targetRast <- raster(extent(targetStack), nrow=nrow(targetStack), ncol=ncol(targetStack))
  targetRast <- writeStart(x=targetRast, filename=paste(outDir, '/', outName, sep='' ), ... )
  
  # get block size
  blocks <- blockSize(targetStack)
  
  # cycle offset by 12 to make it positive
  if (offset < 0) offset <- offset + 12
  
  ##########
  ## MAIN ##
  ##########
  
  # for each block, get temperatures and precipitations and calculate target factor BIOCLIM variables
  for (countBlock in 1:blocks$n) {
    
    if (verbose) cat( paste( 'Calculating target factor BIOCLIM variable for block ', countBlock, ' of ', blocks$n, '...\n', sep='' ) ); flush.console()
    
    # get predictor values, transpose, then make into data frame... each row is a month, each column a different cell
    refVal <- t( getValues(x=refStack, row=blocks$row[countBlock], nrows=blocks$nrows[countBlock] ) )
    targetVal <- t( getValues(x=targetStack, row=blocks$row[countBlock], nrows=blocks$nrows[countBlock] ) )
    
    #######################################
    ### calculate target factor BIOCLIM ###
    #######################################
    
    # if target offset <0 or >0, add rows before/after to target and reference matrices equal to offset value plus/minus 1/0
    if (refDuration > 1) refVal <- rbind(refVal, refVal[1:(refDuration - 1), ])
    if (targetDuration > 1 | offset > 0) targetVal <- rbind(targetVal, targetVal[1:(offset - 1 + targetDuration), ])
    
    # calculate summary statistic for reference values
    if (refDuration != 1) {
      
      refValSummary <- t(
        sapply(
          X=1:12,
          FUN=function(X) apply(refVal[X:(X + refDuration - 1), ], 2, refFunctSummary)
        )
      )
      
    }
    
    # calculate summary statistic for target values
    if (targetDuration != 1) {
      
      targetValSummary <- t(
        sapply(
          X=1:(12 + offset),
          FUN=function(X) apply(targetVal[X:(X + targetDuration - 1), ], 2, targetFunctSummary)
        )
      )
      
    }
    
    # get values of target variable that match criterion WRT reference variable
    nonNa <- which(!is.na(targetVal[1, ]))
    
    theseTargetVal <- rep(NA, ncol(targetVal))
    
    if (length(nonNa) > 0) {
      
      targetValSummaryNonNa <- targetValSummary[ , nonNa]
      refValSummaryNonNa <- refValSummary[ , nonNa]
      
      theseTargetVal[nonNa] <- unlist(
        sapply(X=1:ncol(targetValSummaryNonNa), function(X) targetValSummaryNonNa[which(refValSummaryNonNa[ , X]==refFunctCriterion(refValSummaryNonNa[ , X]))[1] + offset, X])
      )
      
      # round values
      if (roundTarget) theseTargetVal <- round(theseTargetVal)
      
    }
    
    # save to raster
    targetRast <- writeValues(x=targetRast, v=theseTargetVal, start=blocks$row[countBlock])
    
  } # for each block
  
  # stop writing
  targetRast <- writeStop(targetRast)
  
  if (verbose) cat('Done!\n'); flush.console()
  
  return(targetRast)
  
}


predictMax <- function(
  x,
  data,
  outFormat='logistic',
  perm=NULL,
  permLinear=NULL,
  permQuad=NULL,
  permHinge=NULL,
  permThresh=NULL,
  permProd=NULL,
  permProdRule=NULL,
  ...
) {
  # predictMax  Takes a Maxent lambda object or a Maxent object and returns raw or logistic predictions.  Its output is not different from the function predict(maxentModelObject, dataFrame) or predict(maxentModelObject, dataFrame, args='outputformat=raw') (to within rounding error), and in fact those functions should be faster.  However, this function does allow one to perform custom manipulations that those functions do not allow (e.g., permuting product features while leaving other features with the same variables intact).  This function is based on Peter D. Wilson's document "Guidelines for computing MaxEnt model output values from a lambdas file".
  #
  # ARGUMENTS
  # x				Either a Maxent lambda object or a Maxent model object
  # data			Data frame with data to which to make predictions
  # outFormat		Output format: 'raw' ==> produce Maxent raw values; 'logistic' ==> Maxent logistic values
  # perm			Name(s) of variable to permute before calculating predictions. This permutes the variables for *all* features in which they occur.  If a variable is named here, it overrides permutation settings for each feature type.  Note that for product features the variable is permuted before the product is taken.  This permutation is performed before any subsequent permutations (i.e., so if both variables in a product feature are included in perms, then this is equavalent to using the 'beforeProd' rule for permProdRule).  Ignored if NULL.
  # permLinear	Names(s) of variables to permute in linear features before calculating predictions.  Ignored if NULL.
  # permQuad	Names(s) of variables to permute in quadratic features before calculating predictions.  Ignored if NULL.
  # permHinge	Names(s) of variables to permute in forward/reverse hinge features before calculating predictions.  Ignored if NULL.
  # permThresh	Names(s) of variables to permute in threshold features before calculating predictions.  Ignored if NULL.
  # permProd	A list object of n elements, each of which has two character elements listing the variables to permute if they occur in a product feature.  Depending on the value of permProdRule, the function will either permute the individual variables then calculate their product or calculate their product, then permute the product across observations.  Any other features containing the variables will produce values as normal.  Example: permProd=list(c('precipWinter', 'tempWinter'), c('tempSummer', 'precipFall')).  The order of the variables in each element of permProd doesn't matter, so permProd=list(c('temp', 'precip')) is the same as permProd=list(c('precip', 'temp')).  Ignored if NULL.
  # permProdRule Rule for how permutation of product features is applied. 'beforeProd' ==> Permute individual variable values then calculate product. 'afterProd' ==> Calculate product then permute across these values. Ignored if permProd is NULL.
  # ...			Extra arguments... not used.
  #
  # VALUES
  # numeric values representing Maxent predictions
  # 
  # REQUIRED DEPENDANCIES
  # dismo, rJava
  #
  # OPTIONAL DEPENDANCIES
  #
  #
  # BAUHAUS
  # 1. Construct data frame deconstructing lambdas information
  # 2. For each datum (site), predict each feature type
  #
  # EXAMPLE
  # FUNCTION()
  #
  # SOURCE	source('C:/ecology/Drive/R/SDM/SDM - Predict Maxent from Lambdas Object.r')
  #
  # TESTING	Compared to predict() function for many different model forms with permuted data
  #
  #
  # LICENSE
  # This document is copyright ?2014 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/.
  #
  # AUTHOR	Adam B. Smith | Missouri Botanical Garden, St. Louis, Missouri | adamDOTsmithATmobotDOTorg
  # DATE		2015-07-17
  # REVISIONS 
  
  ############################
  ## FUNCTIONS AND PACKAGES ##
  ############################
  
  require(dismo)
  
  ####################
  ## PRE-PROCESSING ##
  ####################
  
  # get lambdas is object is a model object
  if (class(x)=='MaxEnt') x <- x@lambdas
  
  ##########
  ## MAIN ##
  ##########
  
  ### create meta value data frame
  ################################
  metaLines <- x[which(grepl(x, pattern='linearPredictorNormalizer')):length(x)]
  meta <- data.frame(
    item=substr(metaLines, 1, unlist(gregexpr(metaLines, pattern=',')) - 1),
    value=as.numeric(substr(metaLines, unlist(gregexpr(metaLines, pattern=',')) + 2, nchar(metaLines)))
  )
  
  ### create data frame of feature functions
  ##########################################
  feats <- data.frame(
    feature=x[1:(which(grepl(x, pattern='linearPredictorNormalizer')) - 1)],
    type='linear',
    var1=NA,
    var2=NA,
    crit=NA,
    lambda=NA,
    beta1=NA,
    beta2=NA
  )
  
  feats$feature <- as.character(feats$feature)
  feats$type <- as.character(feats$type)
  feats$crit <- as.numeric(feats$crit)
  feats$lambda <- as.numeric(feats$lambda)
  feats$beta1 <- as.numeric(feats$beta1)
  feats$beta2 <- as.numeric(feats$beta2)
  
  ## identify feature types
  feats$type[which(grepl(feats$feature, pattern='[(]') & grepl(feats$feature, pattern='[)]'))] <- 'threshold'
  feats$type[which(grepl(feats$feature, pattern='\\^'))] <- 'quadratic'
  feats$type[which(grepl(feats$feature, pattern='[*]'))] <- 'product'
  feats$type[which(grepl(feats$feature, pattern='\''))] <- 'forward hinge'
  feats$type[which(grepl(feats$feature, pattern='`'))] <- 'reverse hinge'
  
  ## identify (first) variable used
  # see which features have this variable string... note that there is a risk that the parts of a "feature" that have the *substring* including the variable name will falsely imply the variable is in the feature... sorting by length should obviate this
  names <- sort(names(data))
  names <- names[order(nchar(names))] # sort variable names by length
  for (i in 1:length(names)) {
    
    varIn <- grepl(feats$feature, pattern=names[i]) & is.na(feats$var1)
    feats$var1[which(varIn)] <- names[i]
    
  }
  
  ## identify second variable used (if any)... same risks as above
  if ('product' %in% feats$type) {
    indexProd <- which(feats$type=='product')
    for (i in indexProd) {
      for (j in (which(feats$var1[i] == names) + 1):length(names)) {
        secondVarIn <- grepl(feats$feature[i], pattern=names[j])
        if (secondVarIn) feats$var2[i] <- names[j]
      }
    }
  }
  
  ## extract "betas"
  commaPos <- unlist(gregexpr(pattern=',', feats$feature))
  feats$beta1 <- as.numeric(substr(feats$feature, start=commaPos[seq(2, 3 * nrow(feats), by=3)] + 2, stop=commaPos[seq(3, 3 * nrow(feats), by=3)] - 1))
  feats$beta2 <- as.numeric(substr(feats$feature, start=commaPos[seq(3, 3 * nrow(feats), by=3)] + 2, stop=nchar(feats$feature)))
  
  ## extract critical values (hinge/threshold)
  feats$crit[which(feats$type=='threshold')] <- as.numeric(substr(feats$feature[which(feats$type=='threshold')], start=2, stop=unlist(gregexpr(feats$feature[which(feats$type=='threshold')], pattern='<')) - 1))
  
  feats$crit[which(feats$type=='forward hinge')] <- feats$beta1[which(feats$type=='forward hinge')]
  feats$crit[which(feats$type=='reverse hinge')] <- feats$beta2[which(feats$type=='reverse hinge')]
  
  
  ## extract lambdas
  feats$lambda <- as.numeric(substr(feats$feature, start=commaPos[seq(1, 3 * nrow(feats), by=3)] + 2, stop=commaPos[seq(2, 3 * nrow(feats), by=3)] - 1))
  
  ## get values
  value1 <- t(data[ , feats$var1])
  rownames(value1) <- paste0('feat', 1:nrow(feats))
  value2 <- value1 * NA
  
  if ('product' %in% feats$type) {
    value2Compact <- t(data[ , feats$var2[which(feats$type=='product')]])
    value2[which(feats$type=='product'), ] <- value2Compact
  }
  
  # value1 <<- value1
  # value2 <<- value2
  # feats <<- feats
  # print(feats)
  
  ### permutations
  ################
  
  ## permute values REGARDLESS OF FEATURE TYPE
  if (!is.null(perm)) {
    
    for (i in 1:nrow(feats)) {
      if (any(perm %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
      if (any(perm %in% feats$var2[i])) value2[i, ] <- sample(value2[i, ], ncol(value2))
    }
    
  }
  
  ## permute values in LINEAR features
  if (!is.null(permLinear) & 'linear' %in% feats$type) {
    for (i in which(feats$type %in% 'linear')) {
      if (any(permLinear %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
    }
  }
  
  ## permute values in QUADRATIC features
  if (!is.null(permQuad) & 'quadratic' %in% feats$type) {
    for (i in which(feats$type %in% 'quadratic')) {
      if (any(permQuad %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
    }
  }
  
  ## permute values in HINGE features
  if (!is.null(permHinge) & ('forward hinge' %in% feats$type | 'reverse hinge' %in% feats$type)) {
    for (i in sort(unique(c(which(feats$type %in% 'forward hinge'), which(feats$type %in% 'reverse hinge'))))) {
      if (any(permHinge %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
    }
  }
  
  ## permute values in THRESHOLD features
  if (!is.null(permThresh) & 'threshold' %in% feats$type) {
    for (i in which(feats$type %in% 'threshold')) {
      if (any(permThresh %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
    }
  }
  
  ## permute values in PRODUCT features
  if (!is.null(permProd) & 'product' %in% feats$type) {
    
    # permute variables then calculate product
    if (permProdRule=='beforeProd') {
      
      value1Perm <- value1
      value2Perm <- value2
      
      for (i in which(feats$type %in% 'product')) {
        
        # see if permutation is wanted for this particular var1 and var2 of this product feature
        wantPairPerm <- c(sapply(X=permProd, FUN=function(x, y) {x %in% y}, c(feats$var1[i], feats$var2[i])))
        
        if (sum(wantPairPerm)==2) {
          value1Perm[i, ] <- sample(value1Perm[i, ], ncol(value1Perm))
          value2Perm[i, ] <- sample(value2Perm[i, ], ncol(value2Perm))
        }
        
      }
      
      # calculate product
      prodMat <- value1Perm * value2Perm
      
      # calculate product then permute product
    } else if (permProdRule=='afterProd') {
      
      # calculate product
      prodMat <- value1 * value2
      
      for (i in which(feats$type %in% 'product')) {
        
        # see if permutation is wanted for this particular var1 and var2 of this product feature
        wantPairPerm <- c(sapply(X=permProd, FUN=function(x, y) {x %in% y}, c(feats$var1[i], feats$var2[i])))
        if (sum(wantPairPerm)==2) prodMat[i, ] <- sample(prodMat[i, ], ncol(prodMat))
        
      }
      
    }
    
  } else {
    
    # calculate unpermuted product
    prodMat <- value1 * value2
    
  }
  
  ### calculate predictions
  #########################
  
  # initialize exponent
  S <- numeric(nrow(data))
  
  # for each datum
  for (i in 1:nrow(data)) {
    
    # for each type of feature
    for (type in c('linear', 'quadratic', 'product', 'forward hinge', 'reverse hinge', 'threshold')) {
      
      # if feature type occurs in lambdas, evaluate
      if (type %in% feats$type) {
        
        thisValue1 <- value1[which(feats$type==type), i]
        thisValue2 <- value2[which(feats$type==type), i]
        if ('product' %in% feats$type) thisProd <- prodMat[which(feats$type==type), i]
        
        thisLambda <- feats$lambda[which(feats$type==type)]
        thisCrit <- feats$crit[which(feats$type==type)]
        
        thisBeta1 <- feats$beta1[which(feats$type==type)]
        thisBeta2 <- feats$beta2[which(feats$type==type)]
        
        thisRange <- thisBeta2 - thisBeta1
        
        if (type=='linear') S[i] <- S[i] + sum(thisLambda * (thisValue1 - thisBeta1) / thisRange, na.rm=T)
        if (type=='quadratic') S[i] <- S[i] + sum(thisLambda * (thisValue1^2 - thisBeta1) / thisRange, na.rm=T)
        if (type=='product') S[i] <- S[i] + sum(thisLambda * (thisProd - thisBeta1) / thisRange, na.rm=T)
        # if (type=='forward hinge') S[i] <- S[i] + sum(thisLambda[thisValue1 >= thisCrit] * (thisValue1[thisValue1 >= thisCrit] - thisCrit[thisValue1 >= thisCrit]) / thisRange[thisValue1 >= thisCrit], na.rm=T)
        
        
        if (type=='forward hinge') S[i] <- S[i] + sum(thisLambda[thisValue1 >= thisCrit] * (thisValue1[thisValue1 >= thisCrit] - thisCrit[thisValue1 >= thisCrit]) / thisRange[thisValue1 >= thisCrit], na.rm=T)
        if (type=='reverse hinge') S[i] <- S[i] + sum(thisLambda[thisValue1 < thisCrit] * (thisCrit[thisValue1 < thisCrit] - thisValue1[thisValue1 < thisCrit]) / thisRange[thisValue1 < thisCrit], na.rm=T)
        # if (type=='reverse hinge') S[i] <- S[i] + sum(thisLambda[thisValue1 < thisCrit] * (thisCrit[thisValue1 < thisCrit] - thisValue1[thisValue1 < thisCrit]) / thisRange[thisValue1 < thisCrit], na.rm=T)
        
        
        if (type=='threshold') S[i] <- S[i] + sum(thisLambda[which(thisValue1 > thisCrit)], na.rm=T)
        
      }
      
    }
    
  }
  
  raw <- exp(S - meta$value[match('linearPredictorNormalizer', meta$item)]) / meta$value[match('densityNormalizer', meta$item)]
  
  if (outFormat=='raw') {
    return(raw)
  } else {
    scaled <- (raw * exp(meta$value[match('entropy', meta$item)])) / (1 + raw * exp(meta$value[match('entropy', meta$item)]))
    return(scaled)
  }
  
  #####################
  ## POST-PROCESSING ##
  #####################
  
}

# # print(NON)

# set.seed=(57)
# data <- data.frame(
# x=-1 * 1:200,
# y=(1:200)^2,
# z=1:200*runif(200)
# )
# data <- rbind(data, data * -1)

# x <- maxent(data, c(rep(1, 200), rep(0, 200)), args=c('betamultiplier=0.001', 'product=true', 'threshold=false', 'hinge=false', 'quadratic=false', 'jackknife=true'))
# pred <- predict(x, data)
# predMy <- predictMax(x, data)
# plot(pred, predMy)


# predMyAll <- predictMax(x=x, data=data, perm=c('x', 'y', 'z'), permLinear=NULL, permQuad=NULL, permHinge=NULL, permThresh=NULL, permProd=NULL, permProdRule=NULL)
# predMyLinear <- predictMax(x=x, data=data, perm=NULL, permLinear=c('z'), permQuad=NULL, permHinge=NULL, permThresh=NULL, permProd=NULL, permProdRule=NULL)
# predMyQuad <- predictMax(x=x, data=data, perm=NULL, permLinear=NULL, permQuad='x', permHinge=NULL, permThresh=NULL, permProd=NULL, permProdRule=NULL)
# predMyHinge <- predictMax(x=x, data=data, perm=NULL, permLinear=NULL, permQuad=NULL, permHinge='x', permThresh=NULL, permProd=NULL, permProdRule=NULL)
# predMyThresh <- predictMax(x=x, data=data, perm=NULL, permLinear=NULL, permQuad=NULL, permHinge=NULL, permThresh='x', permProd=NULL, permProdRule=NULL)
# predMyProdBefore <- predictMax(x=x, data=data, perm=NULL, permLinear=NULL, permQuad=NULL, permHinge=NULL, permThresh=NULL, permProd=list(c('x', 'y')), permProdRule='beforeProd')
# predMyProdAfter <- predictMax(x=x, data=data, perm=NULL, permLinear=NULL, permQuad=NULL, permHinge=NULL, permThresh=NULL, permProd=list(c('x', 'y')), permProdRule='afterProd')

spokePlot <- function(
  pos,
  neg=NULL,
  ontop='pos',
  labels=NULL,
  shrink=1,
  labelOffset=1.02,
  nudge=1,
  pch=16,
  cexPoints=1,
  cexLabel=1,
  lwdPos=1,
  lwdNeg=1,
  ltyPos='solid',
  ltyNeg='dashed',
  colPos='black',
  colNeg='black',
  colPoints='black',
  colLabel='black',
  ...
) {
  # spokePlot Creates a "spoke" plot to help visualize networks or correlation matrices. Factors are arranged in a circle, and lines connect them if they are linked in some manner (e.g., highly correlated). Two types of spokes can be drawn, one for "positive" associations (denoted with arguments that have "X") and "negative" associations ("Y").
  #
  # ARGUMENTS
  # pos			Binary matrix with 1s indicating row (label) is associated with column (label)
  # neg			As pos, but indicating "negative" associations (however defined)
  # ontop			'pos' ==> plot positive association spokes first; 'neg' ==> plot negative associations first
  # labels		Character vector of names to add to plot.  If NULL then column names of pos will be used.
  # shrink		Numeric, relative size of non-label part of plot... useful is labels are too long to fit onto a plot.  Default = 1.
  # labelOffset	Value (usually > 1) indicating how far from points labels are placed.  If 1, then labels are placed on points and <1 inside points.
  # nudge			Factor by which to multiple y-coordinates of labels. Default is 1. Useful if there are many labels and they tend to overlap one another.
  # pch			Integer, point style (leave as NA to no plot points)
  # cexPoints		Integer, size of points
  # cexLabel	 	Integer, size of labels
  # lwdPos, lwdNeg Integer, line width of spokes
  # ltyPos, ltyNeg Integer or character, line style of spokes (see ?lines)
  # colPos, colNeg Integer or character, color of spokes
  # colPoints		Integer or vector, color of points
  # colLabel		Integer or vector, color of labels
  # ...			Furtehr arguments to pass to plot(), points(), lines(), or text()
  #
  # VALUES
  # None.  By-product is a spoke plot.
  # 
  # REQUIRED DEPENDANCIES
  #
  #
  # OPTIONAL DEPENDANCIES
  #
  #
  # BAUHAUS
  # 
  #
  # EXAMPLE
  # FUNCTION()
  #
  # SOURCE	source('C:/ecology/Drive/Workshops/SDM from Start to Finish (KSU, 2016-02)/Scripts/Spoke Plot.r')
  #
  # TESTING
  #
  #
  # LICENSE
  # This document is copyright ?2014 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/.
  #
  # AUTHOR	Adam B. Smith | Missouri Botanical Garden, St. Louis, Missouri | adamDOTsmithATmobotDOTorg
  # DATE		
  # REVISIONS 
  
  ############################
  ## FUNCTIONS AND PACKAGES ##
  ############################
  
  offset <- 0.2 # amount by which to slide entire plot to left and lower to ensure it's more toward the center of the plotting region
  
  ####################
  ## PRE-PROCESSING ##
  ####################
  
  par(pty='s')
  # plot(x=0, y=0, col=par('bg'), type='n', axes=FALSE, xlab=NA, ylab=NA, xlim=c(-1.1 * shrink * labelOffset, shrink * labelOffset, major * labelOffset, minor * labelOffset), ylim=c(-1.1 * shrink * labelOffset, shrink * labelOffset, major * labelOffset, minor * labelOffset), ...)
  # plot(x=0, y=0, col=par('bg'), type='n', axes=FALSE, xlab=NA, ylab=NA, xlim=c(-1.1 * max(shrink, major) * labelOffset, max(shrink, major) * labelOffset), ylim=c(-1.1 * max(shrink, major) * labelOffset, max(shrink, major) * labelOffset))
  plot(x=0, y=0, col=par('bg'), type='n', axes=FALSE, xlab=NA, ylab=NA, xlim=c(-1.1 * shrink * labelOffset, shrink * labelOffset), ylim=c(-1.1 * shrink * labelOffset, shrink * labelOffset))
  
  ##########
  ## MAIN ##
  ##########
  
  ## get coordinates for connection points
  xLink <- rev(shrink * cos(seq(pi / 2, 2 * pi + pi / 2, length.out=ncol(pos) + 1)))
  yLink <- rev(shrink * sin(seq(pi / 2, 2 * pi + pi / 2, length.out=ncol(pos) + 1)))
  
  ## add points
  if (!is.na(pch)) points(xLink - offset * shrink, yLink, pch=pch, cex=cexPoints, col=colPoints, xpd=NA, ...)
  
  if (ontop=='pos') {
    
    ## add negative spokes
    if (!is.null(neg)) {
      for (i in 1:nrow(neg)) {
        for (j in 1:ncol(neg)) {
          if (!is.na(neg[i, j]) && neg[i, j]==1) lines(x=c(xLink[i] - offset * shrink, xLink[j] - offset * shrink), y=c(yLink[i], yLink[j]), col=colNeg, lwd=lwdNeg, lty=ltyNeg, xpd=NA, ...)
        }
      }
    }
    
    ## add positive spokes
    for (i in 1:nrow(pos)) {
      for (j in 1:ncol(pos)) {
        if (!is.na(pos[i, j]) && pos[i, j]==1) lines(x=c(xLink[i] - offset * shrink, xLink[j] - offset * shrink), y=c(yLink[i], yLink[j]), col=colPos, lwd=lwdPos, lty=ltyPos, xpd=NA, ...)
      }
    }
    
  } else {
    
    ## add positive spokes
    for (i in 1:nrow(pos)) {
      for (j in 1:ncol(pos)) {
        if (!is.na(pos[i, j]) && pos[i, j]==1) lines(x=c(xLink[i] - offset * shrink, xLink[j] - offset * shrink), y=c(yLink[i], yLink[j]), col=colPos, lwd=lwdPos, lty=ltyPos, xpd=NA, ...)
      }
    }
    
    ## add negative spokes
    if (!is.null(neg)) {
      for (i in 1:nrow(neg)) {
        for (j in 1:ncol(neg)) {
          if (!is.na(neg[i, j]) && neg[i, j]==1) lines(x=c(xLink[i] - offset * shrink, xLink[j] - offset * shrink), y=c(yLink[i], yLink[j]), col=colNeg, lwd=lwdNeg, lty=ltyNeg, xpd=NA, ...)
        }
      }
    }
    
  }
  
  ## add labels
  if (is.null(labels)) if (class(pos)=='matrix') { labels <- colnames(pos) } else { labels <- names(pos) }
  
  xLabel <- rev(shrink * labelOffset * cos(seq(pi / 2, 2 * pi + pi / 2, length.out=ncol(pos) + 1)))
  yLabel <- rev(shrink * labelOffset * sin(seq(pi / 2, 2 * pi + pi / 2, length.out=ncol(pos) + 1)))
  yLabel <- yLabel * nudge
  
  position <- rep(1, ncol(pos))
  position[xLabel > 0] <- 4
  position[xLabel < 0] <- 2
  position[xLabel < 10^-5 & xLabel > -10^-5 & yLabel > 0] <- 3
  position[xLabel < 10^-5 & xLabel > -10^-5 & yLabel < 0] <- 1
  
  text(x=xLabel - offset * shrink, y=yLabel, labels=labels, pos=position, cex=cexLabel, col=colLabel, xpd=NA, ...)
  
  #####################
  ## POST-PROCESSING ##
  #####################
  
  
}

maxentAic <- function(
  trainData,
  presentBg,
  betaToTest=betaToTest,
  predStack=NULL,
  scratchDir=NULL,
  params=c('linear=true', 'quadratic=true', 'product=true', 'threshold=true', 'hinge=true', 'addsamplestobackground=false', 'jackknife=true', 'responsecurves=false'),
  predictFx='default',
  threads=1,
  returnModel=TRUE,
  returnAicFrame=TRUE,
  anyway=TRUE,
  verbose=FALSE
) {
  # maxentAic Calibrates master beta in MAXENT model using method described in Warren & Seifert 2011 Ecological niche modeling in Maxent: The importance of model complexity and the performance of model selection criteria.  Ecological Applications 21:335-342.  This script calculates AICc for a series of MAXENT calls across a series of beta values, selects the best beta, then returns a data frame sorted from lowest to highest AICc (then by beta, if there are ties).
  #     The script is able to utilize multiple threads on a computer, speeding things up--however, using a high number of threads on a computer with insufficient memory can really slow things down.  It also uses the disk a lot.  For example, my computer has 16 GB RAM and 8 threads (but a slow hard drive), and using threads >6 causes everything (the mouse, typing, etc.) to be "jumpy".
  #
  # ARGUMENTS
  # trainData
  # data frames with environmental predictors (and no other fields) of training data... training data has presences and background sites
  #
  # presentBg
  # numeric vector of 1/0 to indicate whether line in trainData is a presence or background site
  #
  # predStack
  # raster stack of predictors... leave as NULL to calculate AICc across presence sites and just background sites (which may be the better choice... see Wright et al. 2014 Multiple sources of uncertainty affect metrics for ranking conservation risk under climate change.  Diversity & Distributions 21:111-122.)
  #
  # betaToTest
  # numeric vector of beta values to test (default=1)
  #
  # scratchDir
  # directory to which to write temporary files... leave as NULL to use default directory
  #
  # verbose
  # logical, display extra information
  #
  # params
  # character list of options to pass to maxent() (don't include "betamultiplier"!)
  #
  # predictFx
  #		'default' ==> use standard predict function in raster package
  #		'custom' ==> use custom predictMax() function... this function produces the exact same output as the raster::predict() function except that it dos not handle factors and it extrapolates the Maxent function into non-training areas (i.e., it does not clamp).
  #
  # threads
  # numeric, number of threads for multithreading when writing rasters... can significantly increase calculation speed for large rasters... leave as NULL to use ALL threads available (can slow other processes a lot!)... ignored if not writing rasters!
  #
  # returnModel
  # logical, if TRUE returns model trained using regularization parameter with best AICc... if returnAicFrame is FALSE then returned object will be just the model... if returnAicFrame is TRUE the returned object will be a list with items $model and $aicFrame
  #
  # returnAicFrame
  # logical, if TRUE returns data frame with best beta values... if returnModel is FALSE then returned object will just the this data frame... if returnModel if TRUE then returned object will be a list with items $model and $aicFrame
  #
  # anyway
  # Logical, if no model has fewer coefficients than predictors, return a model anyway with beta=1
  #
  # verbose
  # logical, if TRUE report progress and prints AICc table
  #
  # VALUES
  # AIC-optimized of beta to use for training a MAXENT model
  #
  # BAUHAUS
  # 
  #
  # EXAMPLE
  # require(dismo)
  #
  # ## load maxentAic function (you will need to change the path to the place where you put the file, or just copy and paste that file's contents into R)
  # source('C:/ecology/Drive/r/SDM/SDM - Calibrate MAXENT Model with AIC.r')
  #
  # ## make two rasters to represent predictors
  # pred1 <- pred2 <- raster()
  # pred1[] <- 1:ncell(pred1) # values in first raster go from 1 to total number of cells
  # pred2[] <- runif(ncell(pred2)) # values in second raster ~ normal distribution
  #
  # predStack <- stack(pred1, pred2)
  # names(predStack) <- c('pred1', 'pred2')
  #
  # plot(predStack) # plot rasters
  #
  # ## simulate range of species
  # rangeRast <- pred1 >= 50000 & pred2 > 0.99 # species present in cells with numbers >50000 (pred1) and >0.99 (pred2)
  #
  # plot(rangeRast)
  #
  # ## get cell numbers of all presence sites
  # allPres <- (1:ncell(rangeRast))[c(as.matrix(rangeRast))==1]
  #
  # ## get environment at all presence sites
  # allPresEnv <- data.frame(
  #	pred1=c(as.matrix(pred1))[allPres],
  #	pred2=c(as.matrix(pred2))[allPres]
  # )
  #
  # ## split into training (80%) and testing (20%) presences (note: not really necessary in this example, but typically you would keep aside some data for testing, so we're doing that... if you do supply maxentAic with test data, it should ONLY contain test presences, not presences and absences)
  # trainPres <- allPresEnv[sample(1:nrow(allPresEnv), round(nrow(allPresEnv) * 0.8)), ]
  #
  # ## get environment at background sites
  # bgCoords <- randomPoints(predStack, 1000) # randomly located sites
  # bgEnv <- extract(predStack, bgCoords) # get environment at those sites
  #
  # ## add background to training presences
  # trainData <- rbind(trainPres, bgEnv)
  # presentBg <- c(rep(1, nrow(trainPres)), rep(0, nrow(bgCoords)))
  #
  # ## test regularization parameter values of 1 and 5... note if you're using real rasters which typically have millions of cells, this can take a *very* long time because it must write one raster per beta
  # aicTable <- maxentAic(
  #	trainData=trainData,
  #	presentBg=presentBg,
  #	predStack=predStack,
  #	betaToTest=c(1, 5),
  #	returnModel=TRUE,
  #	returnAicFrame=TRUE,
  #	threads=1,
  #	verbose=TRUE
  # )
  #
  #
  # SOURCE	source('C:/ecology/Drive/R/SDM/SDM - Calibrate MAXENT Model with AIC.r')
  #
  # LICENSE
  # This document is copyright ?2012 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/
  #
  # AUTHOR	Adam B. Smith | Missouri Botanical Garden, St. Louis, Missouri | adamDOTsmithATmobotDOTorg
  # DATE		2012-01
  # REVISIONS 2012-02-22  Test presence data frame can be left as NULL
  #			2012-02-24  Fixed bug associated with data input as vectors, not data frames
  #			2012-06-28  Output is now optimal beta, not model with optimal beta
  #			2012-09-21	Use can state number of cores to use in training
  #			2012-09-27	Returns data frame of results for all betas anstatt just best beta
  #			2012-10-15  Fixed calculations of AICc so they're calculated across all training and test sites, not background sites
  # 			2012-11-30	User can now specify which feature functions to use
  #			2013-07-25  Corrected calculation for log likelihood (thanks to Peter Galante, City College of New York)
  #			2014-02-18  Removed multithreading (not useful except for jackknifing, which isn't used here anyway).
  #			2014-02-26  Added ability to use/not use threshold features
  #			2014-04-07  Added optional multi-core capacity
  #			2014-05-09  User can have script return best model OR data frame with AICc for each model tested
  #			2014-05-16  User can return best model and/or AIC frame
  #			2014-10-15  Fixed bug in handling training data with just one variable (this needs to be handled as a data frame, but was converting to a vector)
  #			2014-11-10	User can now speed calculations by using approximate log likelihood
  #			2015-02-05  User can now calculates AICc using only training/test presences and training background sites (see Wright et al. 2014 Multiple sources of uncertainty affect metrics for ranking conservation risk under climate change.  Diversity & Distributions 21:111-122.)  Also removed testpres argument and assuming trainData has all presences (training and test) in it.
  #			2015-07-09  Changed "params" argument from list to string.
  #			2015-11-05  Added capability to return model with beta=1 if number of predictors > number of coefficients for all models
  
  if (!returnModel & !returnAicFrame) warning('maxentAic function will return empty object because both returnModel and returnAicFrame are FALSE')
  
  # switch off aspects not needed for initial run of all models
  if (any(grepl(params, pattern='jackknife=true'))) {
    jack <- TRUE
    params[which(params=='jackknife=true')] <- 'jackknife=false'
  } else { jack <- FALSE }
  
  if (any(grepl(params, pattern='responsecurves=true'))) {
    resp <- TRUE
    params[which(params=='responsecurves=true')] <- 'responsecurves=false'
  } else { resp <- FALSE }
  
  #############################
  ## libraries and functions ##
  #############################
  
  require(rJava)
  require(dismo)
  # source('C:/ecology/Drive/r/SDM/SDM - Write Raster from Model Basic Function.r')
  # source('C:/ecology/Drive/R/SDM/SDM - Predict Maxent from Lambdas Object.r')
  
  if (is.null(threads)) {
    
    require(parallel)
    threads <- min(detectCores(), length(betaToTest))
    
  } else if (threads > 1) {
    
    require(parallel)
    threads <- min(threads, detectCores(), length(betaToTest))
    
  }
  
  ###########
  ## setup ##
  ###########
  
  # create scratch directory
  scratchDir <- if (is.null(scratchDir)) { paste0(getwd(), '/_scratchDir_dontDelete') } else { scratchDir }
  dir.create(scratchDir, showWarnings=FALSE, recursive=TRUE)
  
  ## collate all presences
  allPresences <- trainData[presentBg==1, ]
  if (class(allPresences)=='factor' | class(allPresences)=='numeric' | class(allPresences)=='integer') allPresences <- data.frame(allPresences)
  names(allPresences) <- names(trainData) # names of data to which to predict
  
  ## collate all background sites
  allBg <- trainData[presentBg==0, ]
  if (class(allBg)=='numeric') allBg <- as.data.frame(allBg)
  names(allBg) <- names(trainData)
  
  ##########
  ## MAIN ##
  ##########
  
  ### SINGLE CORE
  if (threads==1) {
    
    if (verbose) cat('Calculating AICc for beta: ')
    
    for (thisBeta in betaToTest) { # for each beta
      
      if (verbose) { cat(thisBeta, ' '); flush.console() }
      
      # generate random number for temp folder
      thisRand <- round(runif(1) * 10^12)
      dir.create(paste0(scratchDir, '/', thisRand), showWarnings=TRUE, recursive=TRUE)
      
      # train model
      trialModel <- maxent(
        x=trainData,
        p=as.vector(presentBg),
        removeDuplicates=FALSE,
        path=paste0(scratchDir, '/', thisRand),
        args=c(
          paste0('betamultiplier=', thisBeta),
          params
        )
      )
      
      ## predict to training (and maybe test presences)
      predPres <- if (predictFx=='default') {
        
        predict(
          object=trialModel,
          x=allPresences,
          na.rm=TRUE,
          args='outputformat=raw'
        )
        
      } else if (predictFx=='custom') {
        
        predictMax(
          x=trialModel,
          data=allPresences,
          outFormat='raw'
        )			
        
      }
      
      
      ## predict to background
      if (is.null(predStack)) {
        
        predBg <- if (predictFx=='default') {
          
          predict(
            object=trialModel,
            x=allBg,
            na.rm=TRUE,
            args='outputformat=raw'
          )
          
        } else if (predictFx=='custom') {
          
          predictMax(
            x=trialModel,
            data=allBg,
            outFormat='raw'
          )			
          
        }
        
        bgSum <- sum(predBg)
        
        # predict to raster
      } else {
        
        thisRaster <- genericWriteRaster(
          predictorStack=predStack,
          theModel=trialModel,
          fileName=paste0(scratchDir, '/', thisRand, '/mxTuning_beta_', round(100 * thisBeta)),
          predArgs=list(
            maxent=list(
              outputFormat='raw',
              predictFx=predictFx
            )
          ),
          na.rm=T,
          format='raster',
          overwrite=TRUE
        )
        
        bgSum <- cellStats(thisRaster, 'sum')
        
        rm(thisRaster); gc()
        
      }
      
      # waste time... averts a write error, somehow
      Sys.sleep(1)
      
      # delete temp dir
      unlink(paste0(scratchDir, '/', thisRand), recursive=TRUE)
      
      ## calculate log likelihood
      logLik <- sum(log(predPres / bgSum), na.rm=TRUE)
      
      ## calculate number of parameters
      K <- 0
      
      for (thisLambda in Lambdas <- trialModel@lambdas) { # for each line in lambda object
        
        commaPos <- gregexpr( text=thisLambda, pattern=',') # get location of commas
        
        if (length(commaPos[[1]]) > 1) { # if there is >1 comma in this line (this is not a parameter line)
          
          paramValue <- as.numeric( substr(x=thisLambda, start=commaPos[[1]][1]+1, stop=commaPos[[1]][2]-1 ) ) # convert string between first two commas to numeric
          
          if (paramValue !=0) K <- K + 1 # increment number of parameters
          
        } # if there is >1 comma in this line
        
      }
      
      # AICc
      AICc <- -2 * logLik + 2 * K + (2 * K * (K + 1)) / ( sum(presentBg) - K - 1)
      
      # remember
      thisAicFrame <- data.frame(
        beta=thisBeta,
        n=sum(presentBg),
        logLik=logLik,
        K=K,
        AICc=AICc
      )
      
      if (exists('aicFrame', inherits=FALSE)) aicFrame <- rbind(aicFrame, thisAicFrame)
      if (!exists('aicFrame', inherits=FALSE)) aicFrame <- thisAicFrame
      
    } # next beta	
    
    ######################
    ### MULTI-THREADED ###
    ######################
    
  } else {
    
    if (verbose) { cat('Using multiple threads...\n'); flush.console() }
    
    # function to be passed to threads
    workerFunction <- function(trainData, presentBg, thisBeta, params, allPresences, predStack, scratchDir, thisRand=thisRand) {
      
      # generate random number for temp folder
      dir.create(paste0(scratchDir, '/', thisRand), showWarnings=TRUE, recursive=TRUE)
      
      # train model
      trialModel <- maxent(
        x=trainData,
        p=as.vector(presentBg),
        removeDuplicates=FALSE,
        path=paste0(scratchDir, '/', thisRand),
        args=c(
          paste0('betamultiplier=', thisBeta),
          params
        )
      )
      
      ## predict to training (and test presences)
      predPres <- if (predictFx=='default') {
        
        predict(
          object=trialModel,
          x=allPresences,
          na.rm=TRUE,
          args='outputformat=raw'
        )
        
      } else if (predictFx=='custom') {
        
        predictMax(
          x=trialModel,
          data=allPresences,
          outFormat='raw'
        )			
        
      }
      ## predict to background
      if (is.null(predStack)) {
        
        predBg <- if (predictFx=='default') {
          
          predict(
            object=trialModel,
            x=allBg,
            na.rm=TRUE,
            args='outputformat=raw'
          )
          
        } else if (predictFx=='custom') {
          
          predictMax(
            x=trialModel,
            data=allBg,
            outFormat='raw'
          )			
          
        }
        
        bgSum <- sum(predBg)
        
        # predict to raster
      } else {
        
        thisRaster <- genericWriteRaster(
          predictorStack=predStack,
          theModel=trialModel,
          fileName=paste0(scratchDir, '/', thisRand, '/mxTuning_beta_', round(100 * thisBeta)),
          predArgs=list(
            maxent=list(
              outputFormat='raw',
              predictFx=predictFx
            )
          ),
          na.rm=T,
          format='raster',
          overwrite=TRUE
        )
        
        bgSum <- cellStats(thisRaster, 'sum')
        
        rm(thisRaster); gc()
        
      }
      
      
      # waste time... averts a write error, somehow
      Sys.sleep(1)
      
      # delete temp dir
      unlink(paste0(scratchDir, '/', thisRand), recursive=TRUE)
      
      ## calculate log likelihood
      
      logLik <- sum(log(predPres / mapSum), na.rm=TRUE)
      
      ## calculate number of parameters
      K <- 0
      
      for (thisLambda in Lambdas <- trialModel@lambdas) { # for each line in lambda object
        
        commaPos <- gregexpr( text=thisLambda, pattern=',') # get location of commas
        
        if (length(commaPos[[1]]) > 1) { # if there is >1 comma in this line (this is not a parameter line)
          
          paramValue <- as.numeric( substr(x=thisLambda, start=commaPos[[1]][1]+1, stop=commaPos[[1]][2]-1 ) ) # convert string between first two commas to numeric
          
          if (paramValue !=0) K <- K + 1 # increment number of parameters
          
        } # if there is >1 comma in this line
        
      }
      
      # AICc
      AICc <- -2 * logLik + 2 * K + (2 * K * (K + 1)) / ( sum(presentBg) - K - 1)
      
      # remember
      thisAicFrame <- data.frame(
        beta=thisBeta,
        n=sum(presentBg),
        logLik=logLik,
        K=K,
        AICc=AICc
      )
      
      return(thisAicFrame)
      
    }
    
    # initiate cluster
    beginCluster(threads)
    thisCluster <- getCluster()
    # on.exit(returnCluster())
    
    # get all nodes going
    for (i in 1:threads) {
      
      # generate random number for raster name
      thisRand <- round(runif(1) * 10^12)
      
      sendCall(thisCluster[[i]], workerFunction, args=list(trainData=trainData, presentBg=presentBg, thisBeta=betaToTest[i], params=params, allPresences=allPresences, predStack=predStack, scratchDir=scratchDir, thisRand=thisRand, predictMax=predictMax))
      
    }
    
    for (i in 1:length(betaToTest)) {
      
      # receive results from a node
      receivedData <- recvOneData(thisCluster)
      
      # cat('Received i: ', i, '\n'); flush.console()
      # cat('Received beta: ', betaToTest[i], '\n\n'); flush.console()
      
      # error?
      if (!receivedData$value$success) stop('\nCluster error in maxentAic function.')
      
      if (exists('aicFrame', inherits=F)) aicFrame <- rbind(aicFrame, receivedData$value$value)
      if (!exists('aicFrame', inherits=F)) aicFrame <- receivedData$value$value
      
      # need to send more data?
      if ((threads + i) <= length(betaToTest)) {
        
        # generate random number for raster name
        thisRand <- round(runif(1) * 10^12)
        
        sendCall(thisCluster[[receivedData$node]], workerFunction, args=list(trainData=trainData, presentBg=presentBg, thisBeta=betaToTest[threads + i], params=params, allPresences=allPresences, predStack=predStack, scratchDir=scratchDir, thisRand=thisRand, predictMax=predictMax))
        
      }
      
    }
    
    endCluster()	
    
  } # end multi-thread
  
  # remove models with more parameters than data points that have more than 0 parameters
  aicFrame <- subset(aicFrame, aicFrame$n >= aicFrame$K & aicFrame$K > 0) 
  aicFrame <- aicFrame[order(aicFrame$K), ] # do this to ensure if any AICc are the same across K that one with the smaller K gets chosen
  aicFrame <- aicFrame[order(aicFrame$beta, decreasing=TRUE), ] # do this to ensure if any AICc are the same across beta that one with the larger beta gets chosen
  aicFrame <- aicFrame[order(aicFrame$AICc), ]
  
  aicFrame$deltaAIcc <- aicFrame$AICc - min(aicFrame$AICc)
  aicFrame$relLike <- exp(-0.5 * aicFrame$deltaAIcc)
  aicFrame$aicWeight <- aicFrame$relLike / sum(aicFrame$relLike)
  
  if (verbose) {
    
    cat('\n')
    print(aicFrame)
    flush.console()
    
  }
  
  # if user wants best model returned
  if (returnModel) {
    
    if (jack) params[which(params=='jackknife=false')] <- 'jackknife=true'
    if (resp) params[which(params=='responsecurves=false')] <- 'responsecurves=true'
    
    # train model
    if (nrow(aicFrame) > 0) {
      
      model <- maxent(
        x=trainData,
        p=as.vector(presentBg),
        removeDuplicates=FALSE,
        path=paste0(scratchDir, '/', thisRand),
        args=c(
          paste0('betamultiplier=', ifelse(all(is.infinite(aicFrame$AICc)), aicFrame$beta[which.max(aicFrame$logLik)], aicFrame$beta[which.min(aicFrame$AICc)])),
          params
        )
      )
      
      if (!resp) unlink(paste0(scratchDir, '/', thisRand), recursive=TRUE)
      
    } else if (anyway) {
      
      warning('Returning model with beta = 1 even though no model had fewer coefficients than predictors.', immediate.=TRUE)
      
      model <- maxent(
        x=trainData,
        p=as.vector(presentBg),
        removeDuplicates=FALSE,
        path=paste0(scratchDir, '/', thisRand),
        args=c('betamultiplier=1', params)
      )
      
      if (!resp) unlink(paste0(scratchDir, '/', thisRand), recursive=TRUE)
      
    } else {
      
      warning('No models had fewer coefficients than predictors. No model returned.', immediate.=TRUE)
      model <- 'No MAXENT model had number of parameters < number of training presences.'
      
    }
    
  }
  
  # return stuff
  if (returnModel & !returnAicFrame) {
    return(model)
  } else if (!returnModel & returnAicFrame) {
    return(aicFrame)
  } else {
    return(list(aicFrame=aicFrame, model=model))
  }
  
}


## function to sample raster with/out replacement with probabilities proportional to raster values
sampleRast <- function(mask, n, replace, prob=TRUE) {
  
  # replacement function for randomPoints() allowing user to specify if sampling by replacement is allowed
  #
  # mask		raster
  # n			number of points
  # replace	TRUE ==> sample with replacement; FALSE ==> sample without replacement
  # prob		TRUE ==> sample with regards to cell values
  
  val <- as.vector(mask)
  cellNum <- 1:ncell(mask)
  cellNum <- cellNum[!is.na(val)]
  if (prob) val <- val[!is.na(val)]
  
  s <- if (prob) {
    sample(cellNum, size=n, replace=replace, prob=val)
  } else {
    sample(cellNum, size=n, replace=replace)
  }
  xy <- xyFromCell(mask, s)
  
  rm(val, cellNum)
  
  return(xy)
  
}

geogThinApprox <- function(
  x,
  minDist,
  longLatFields=c('longitude', 'latitude'),
  distFunct=NULL,
  verbose=0,
  ...
) {
  # geogThinApprox  Thins a set of geographic points so of the remainder, none are closer than a given distance (< minDist).  The function is nearlt the same as thin.algorithm() in the spThin package, except that it accepts a data frame or SpatialPoints or SpatialPointsDataFrame as a main argument and the user can specify the distance function to be used.  Its key advantage over thin.algorithm() is that 1) it returns the points plus any associated data, whereas that function only returns points; and 2) it is faster, especially for large data sets.
  #
  # ARGUMENTS
  # x				A data frame or a SpatialPoints or SpatialPointsDataFrame object with spatial points to be thinned.
  # minDist		Minimum distance (usually in m) returned points must be from one another.
  # longLatFields If x is a data frame, this is a 2-element charcater vector with the name of the longitude and latitude fields.  Ignored if x is not a data frame.
  # distFunct		If NULL then distCosine() in the geoshere package is used to calculate distances.  More accurate distances can be obtained by using other functions (see ?distCosine).
  # verbose		Numeric, if 0 then silent, if >0 then show progress.
  # ...			Extra arguments to pass to distFunct
  #
  # VALUES
  # A data frame or SpatialPoints* object with geographically thinned points.
  # 
  # 
  # REQUIRED DEPENDANCIES
  # 
  #
  # OPTIONAL DEPENDANCIES
  # geosphere (if distFunct is one of the dist* formulas in that package)
  #
  # BAUHAUS
  # 
  #
  # EXAMPLE
  # FUNCTION()
  #
  # SOURCE	source('C:/ecology/Drive/R/Geography/Points - Thin Geographic Points by Given Minimum Distance (Approximate).r')
  #
  # TESTING
  #
  #
  # LICENSE
  # This document is copyright ?2014 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/.
  #
  # AUTHOR	Adam B. Smith | Missouri Botanical Garden, St. Louis, Missouri | adamDOTsmithATmobotDOTorg
  # DATE		2015-10-15
  # REVISIONS 
  
  ############################
  ## FUNCTIONS AND PACKAGES ##
  ############################
  
  if (is.null(distFunct)) {
    library(geosphere)
    distFunct <- distCosine
  }
  
  ####################
  ## PRE-PROCESSING ##
  ####################
  
  if (verbose > 0) cat('Thinning', nrow(x), 'points geographically by distance of', minDist, '(probably meters)\n'); flush.console()
  
  # get coorindates
  coords <- if (class(x)=='data.frame') { x[ , longLatFields] } else { coordinates(x) }
  index <- 1:nrow(coords)
  
  # calculate distances between points
  dists <- matrix(rep(NA, nrow(coords)^2), nrow=nrow(coords))
  
  if (verbose > 0) say('Calculating pairwise distances from point:', post=0)
  for (i in 1:nrow(coords)) {
    if (verbose > 0) if ((10 * i / nrow(coords)) %% 1 == 0) say('=', post=ifelse(i == nrow(coords), 1, 0), say='')
    dists[i, i:ncol(dists)] <- distFunct(coords[i, ], coords[i:ncol(dists), ])
    dists[i, 1:i] <- dists[1:i, i]
  }
  if (verbose > 0) say('')
  diag(dists) <- NA
  dists <- dists < minDist
  
  # vector to store number of points too close to each particular point
  tooClose <- colSums(dists, na.rm=T)
  
  ##########
  ## MAIN ##
  ##########
  
  Sys.time()
  while (length(index) > 1 & sum(tooClose) > 0) {
    
    indexMaxNeigh <- index[tooClose==max(tooClose)]
    
    # if just one point too close to others, sample among them
    if (length(indexMaxNeigh) > 1) indexMaxNeigh <- indexMaxNeigh[sample(1:length(indexMaxNeigh), 1)]
    
    # remove point
    index <- index[-which(index==indexMaxNeigh)]
    tooClose <- colSums(dists[index, index], na.rm=T)
    
  }
  Sys.time()
  
  
  #####################
  ## POST-PROCESSING ##
  #####################
  
  x <- x[index, ]
  
  return(x)
  
}

maxentAic <- function(
  trainData,
  presentBg,
  betaToTest=betaToTest,
  predStack=NULL,
  scratchDir=NULL,
  params=c('linear=true', 'quadratic=true', 'product=true', 'threshold=true', 'hinge=true', 'addsamplestobackground=false', 'jackknife=true', 'responsecurves=false'),
  predictFx='default',
  threads=1,
  returnModel=TRUE,
  returnAicFrame=TRUE,
  anyway=TRUE,
  verbose=FALSE
) {
  # maxentAic Calibrates master beta in MAXENT model using method described in Warren & Seifert 2011 Ecological niche modeling in Maxent: The importance of model complexity and the performance of model selection criteria.  Ecological Applications 21:335-342.  This script calculates AICc for a series of MAXENT calls across a series of beta values, selects the best beta, then returns a data frame sorted from lowest to highest AICc (then by beta, if there are ties).
  #     The script is able to utilize multiple threads on a computer, speeding things up--however, using a high number of threads on a computer with insufficient memory can really slow things down.  It also uses the disk a lot.  For example, my computer has 16 GB RAM and 8 threads (but a slow hard drive), and using threads >6 causes everything (the mouse, typing, etc.) to be "jumpy".
  #
  # ARGUMENTS
  # trainData
  # data frames with environmental predictors (and no other fields) of training data... training data has presences and background sites
  #
  # presentBg
  # numeric vector of 1/0 to indicate whether line in trainData is a presence or background site
  #
  # predStack
  # raster stack of predictors... leave as NULL to calculate AICc across presence sites and just background sites (which may be the better choice... see Wright et al. 2014 Multiple sources of uncertainty affect metrics for ranking conservation risk under climate change.  Diversity & Distributions 21:111-122.)
  #
  # betaToTest
  # numeric vector of beta values to test (default=1)
  #
  # scratchDir
  # directory to which to write temporary files... leave as NULL to use default directory
  #
  # verbose
  # logical, display extra information
  #
  # params
  # character list of options to pass to maxent() (don't include "betamultiplier"!)
  #
  # predictFx
  #		'default' ==> use standard predict function in raster package
  #		'custom' ==> use custom predictMax() function... this function produces the exact same output as the raster::predict() function except that it dos not handle factors and it extrapolates the Maxent function into non-training areas (i.e., it does not clamp).
  #
  # threads
  # numeric, number of threads for multithreading when writing rasters... can significantly increase calculation speed for large rasters... leave as NULL to use ALL threads available (can slow other processes a lot!)... ignored if not writing rasters!
  #
  # returnModel
  # logical, if TRUE returns model trained using regularization parameter with best AICc... if returnAicFrame is FALSE then returned object will be just the model... if returnAicFrame is TRUE the returned object will be a list with items $model and $aicFrame
  #
  # returnAicFrame
  # logical, if TRUE returns data frame with best beta values... if returnModel is FALSE then returned object will just the this data frame... if returnModel if TRUE then returned object will be a list with items $model and $aicFrame
  #
  # anyway
  # Logical, if no model has fewer coefficients than predictors, return a model anyway with beta=1
  #
  # verbose
  # logical, if TRUE report progress and prints AICc table
  #
  # VALUES
  # AIC-optimized of beta to use for training a MAXENT model
  #
  # BAUHAUS
  # 
  #
  # EXAMPLE
  # require(dismo)
  #
  # ## load maxentAic function (you will need to change the path to the place where you put the file, or just copy and paste that file's contents into R)
  # source('C:/ecology/Drive/r/SDM/SDM - Calibrate MAXENT Model with AIC.r')
  #
  # ## make two rasters to represent predictors
  # pred1 <- pred2 <- raster()
  # pred1[] <- 1:ncell(pred1) # values in first raster go from 1 to total number of cells
  # pred2[] <- runif(ncell(pred2)) # values in second raster ~ normal distribution
  #
  # predStack <- stack(pred1, pred2)
  # names(predStack) <- c('pred1', 'pred2')
  #
  # plot(predStack) # plot rasters
  #
  # ## simulate range of species
  # rangeRast <- pred1 >= 50000 & pred2 > 0.99 # species present in cells with numbers >50000 (pred1) and >0.99 (pred2)
  #
  # plot(rangeRast)
  #
  # ## get cell numbers of all presence sites
  # allPres <- (1:ncell(rangeRast))[c(as.matrix(rangeRast))==1]
  #
  # ## get environment at all presence sites
  # allPresEnv <- data.frame(
  #	pred1=c(as.matrix(pred1))[allPres],
  #	pred2=c(as.matrix(pred2))[allPres]
  # )
  #
  # ## split into training (80%) and testing (20%) presences (note: not really necessary in this example, but typically you would keep aside some data for testing, so we're doing that... if you do supply maxentAic with test data, it should ONLY contain test presences, not presences and absences)
  # trainPres <- allPresEnv[sample(1:nrow(allPresEnv), round(nrow(allPresEnv) * 0.8)), ]
  #
  # ## get environment at background sites
  # bgCoords <- randomPoints(predStack, 1000) # randomly located sites
  # bgEnv <- extract(predStack, bgCoords) # get environment at those sites
  #
  # ## add background to training presences
  # trainData <- rbind(trainPres, bgEnv)
  # presentBg <- c(rep(1, nrow(trainPres)), rep(0, nrow(bgCoords)))
  #
  # ## test regularization parameter values of 1 and 5... note if you're using real rasters which typically have millions of cells, this can take a *very* long time because it must write one raster per beta
  # aicTable <- maxentAic(
  #	trainData=trainData,
  #	presentBg=presentBg,
  #	predStack=predStack,
  #	betaToTest=c(1, 5),
  #	returnModel=TRUE,
  #	returnAicFrame=TRUE,
  #	threads=1,
  #	verbose=TRUE
  # )
  #
  #
  # SOURCE	source('C:/ecology/Drive/R/SDM/SDM - Calibrate MAXENT Model with AIC.r')
  #
  # LICENSE
  # This document is copyright ?2012 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/
  #
  # AUTHOR	Adam B. Smith | Missouri Botanical Garden, St. Louis, Missouri | adamDOTsmithATmobotDOTorg
  # DATE		2012-01
  # REVISIONS 2012-02-22  Test presence data frame can be left as NULL
  #			2012-02-24  Fixed bug associated with data input as vectors, not data frames
  #			2012-06-28  Output is now optimal beta, not model with optimal beta
  #			2012-09-21	Use can state number of cores to use in training
  #			2012-09-27	Returns data frame of results for all betas anstatt just best beta
  #			2012-10-15  Fixed calculations of AICc so they're calculated across all training and test sites, not background sites
  # 			2012-11-30	User can now specify which feature functions to use
  #			2013-07-25  Corrected calculation for log likelihood (thanks to Peter Galante, City College of New York)
  #			2014-02-18  Removed multithreading (not useful except for jackknifing, which isn't used here anyway).
  #			2014-02-26  Added ability to use/not use threshold features
  #			2014-04-07  Added optional multi-core capacity
  #			2014-05-09  User can have script return best model OR data frame with AICc for each model tested
  #			2014-05-16  User can return best model and/or AIC frame
  #			2014-10-15  Fixed bug in handling training data with just one variable (this needs to be handled as a data frame, but was converting to a vector)
  #			2014-11-10	User can now speed calculations by using approximate log likelihood
  #			2015-02-05  User can now calculates AICc using only training/test presences and training background sites (see Wright et al. 2014 Multiple sources of uncertainty affect metrics for ranking conservation risk under climate change.  Diversity & Distributions 21:111-122.)  Also removed testpres argument and assuming trainData has all presences (training and test) in it.
  #			2015-07-09  Changed "params" argument from list to string.
  #			2015-11-05  Added capability to return model with beta=1 if number of predictors > number of coefficients for all models
  
  if (!returnModel & !returnAicFrame) warning('maxentAic function will return empty object because both returnModel and returnAicFrame are FALSE')
  
  # switch off aspects not needed for initial run of all models
  if (any(grepl(params, pattern='jackknife=true'))) {
    jack <- TRUE
    params[which(params=='jackknife=true')] <- 'jackknife=false'
  } else { jack <- FALSE }
  
  if (any(grepl(params, pattern='responsecurves=true'))) {
    resp <- TRUE
    params[which(params=='responsecurves=true')] <- 'responsecurves=false'
  } else { resp <- FALSE }
  
  #############################
  ## libraries and functions ##
  #############################
  
  require(rJava)
  require(dismo)
  # source('C:/ecology/Drive/r/SDM/SDM - Write Raster from Model Basic Function.r')
  # source('C:/ecology/Drive/R/SDM/SDM - Predict Maxent from Lambdas Object.r')
  
  if (is.null(threads)) {
    
    require(parallel)
    threads <- min(detectCores(), length(betaToTest))
    
  } else if (threads > 1) {
    
    require(parallel)
    threads <- min(threads, detectCores(), length(betaToTest))
    
  }
  
  ###########
  ## setup ##
  ###########
  
  # create scratch directory
  scratchDir <- if (is.null(scratchDir)) { paste0(getwd(), '/_scratchDir_dontDelete') } else { scratchDir }
  dir.create(scratchDir, showWarnings=FALSE, recursive=TRUE)
  
  ## collate all presences
  allPresences <- trainData[presentBg==1, ]
  if (class(allPresences)=='factor' | class(allPresences)=='numeric' | class(allPresences)=='integer') allPresences <- data.frame(allPresences)
  names(allPresences) <- names(trainData) # names of data to which to predict
  
  ## collate all background sites
  allBg <- trainData[presentBg==0, ]
  if (class(allBg)=='numeric') allBg <- as.data.frame(allBg)
  names(allBg) <- names(trainData)
  
  ##########
  ## MAIN ##
  ##########
  
  ### SINGLE CORE
  if (threads==1) {
    
    if (verbose) cat('Calculating AICc for beta: ')
    
    for (thisBeta in betaToTest) { # for each beta
      
      if (verbose) { cat(thisBeta, ' '); flush.console() }
      
      # generate random number for temp folder
      thisRand <- round(runif(1) * 10^12)
      dir.create(paste0(scratchDir, '/', thisRand), showWarnings=TRUE, recursive=TRUE)
      
      # train model
      trialModel <- maxent(
        x=trainData,
        p=as.vector(presentBg),
        removeDuplicates=FALSE,
        path=paste0(scratchDir, '/', thisRand),
        args=c(
          paste0('betamultiplier=', thisBeta),
          params
        )
      )
      
      ## predict to training (and maybe test presences)
      predPres <- if (predictFx=='default') {
        
        predict(
          object=trialModel,
          x=allPresences,
          na.rm=TRUE,
          args='outputformat=raw'
        )
        
      } else if (predictFx=='custom') {
        
        predictMax(
          x=trialModel,
          data=allPresences,
          outFormat='raw'
        )			
        
      }
      
      
      ## predict to background
      if (is.null(predStack)) {
        
        predBg <- if (predictFx=='default') {
          
          predict(
            object=trialModel,
            x=allBg,
            na.rm=TRUE,
            args='outputformat=raw'
          )
          
        } else if (predictFx=='custom') {
          
          predictMax(
            x=trialModel,
            data=allBg,
            outFormat='raw'
          )			
          
        }
        
        bgSum <- sum(predBg)
        
        # predict to raster
      } else {
        
        thisRaster <- genericWriteRaster(
          predictorStack=predStack,
          theModel=trialModel,
          fileName=paste0(scratchDir, '/', thisRand, '/mxTuning_beta_', round(100 * thisBeta)),
          predArgs=list(
            maxent=list(
              outputFormat='raw',
              predictFx=predictFx
            )
          ),
          na.rm=T,
          format='raster',
          overwrite=TRUE
        )
        
        bgSum <- cellStats(thisRaster, 'sum')
        
        rm(thisRaster); gc()
        
      }
      
      # waste time... averts a write error, somehow
      Sys.sleep(1)
      
      # delete temp dir
      unlink(paste0(scratchDir, '/', thisRand), recursive=TRUE)
      
      ## calculate log likelihood
      logLik <- sum(log(predPres / bgSum), na.rm=TRUE)
      
      ## calculate number of parameters
      K <- 0
      
      for (thisLambda in Lambdas <- trialModel@lambdas) { # for each line in lambda object
        
        commaPos <- gregexpr( text=thisLambda, pattern=',') # get location of commas
        
        if (length(commaPos[[1]]) > 1) { # if there is >1 comma in this line (this is not a parameter line)
          
          paramValue <- as.numeric( substr(x=thisLambda, start=commaPos[[1]][1]+1, stop=commaPos[[1]][2]-1 ) ) # convert string between first two commas to numeric
          
          if (paramValue !=0) K <- K + 1 # increment number of parameters
          
        } # if there is >1 comma in this line
        
      }
      
      # AICc
      AICc <- -2 * logLik + 2 * K + (2 * K * (K + 1)) / ( sum(presentBg) - K - 1)
      
      # remember
      thisAicFrame <- data.frame(
        beta=thisBeta,
        n=sum(presentBg),
        logLik=logLik,
        K=K,
        AICc=AICc
      )
      
      if (exists('aicFrame', inherits=FALSE)) aicFrame <- rbind(aicFrame, thisAicFrame)
      if (!exists('aicFrame', inherits=FALSE)) aicFrame <- thisAicFrame
      
    } # next beta	
    
    ######################
    ### MULTI-THREADED ###
    ######################
    
  } else {
    
    if (verbose) { cat('Using multiple threads...\n'); flush.console() }
    
    # function to be passed to threads
    workerFunction <- function(trainData, presentBg, thisBeta, params, allPresences, predStack, scratchDir, thisRand=thisRand) {
      
      # generate random number for temp folder
      dir.create(paste0(scratchDir, '/', thisRand), showWarnings=TRUE, recursive=TRUE)
      
      # train model
      trialModel <- maxent(
        x=trainData,
        p=as.vector(presentBg),
        removeDuplicates=FALSE,
        path=paste0(scratchDir, '/', thisRand),
        args=c(
          paste0('betamultiplier=', thisBeta),
          params
        )
      )
      
      ## predict to training (and test presences)
      predPres <- if (predictFx=='default') {
        
        predict(
          object=trialModel,
          x=allPresences,
          na.rm=TRUE,
          args='outputformat=raw'
        )
        
      } else if (predictFx=='custom') {
        
        predictMax(
          x=trialModel,
          data=allPresences,
          outFormat='raw'
        )			
        
      }
      ## predict to background
      if (is.null(predStack)) {
        
        predBg <- if (predictFx=='default') {
          
          predict(
            object=trialModel,
            x=allBg,
            na.rm=TRUE,
            args='outputformat=raw'
          )
          
        } else if (predictFx=='custom') {
          
          predictMax(
            x=trialModel,
            data=allBg,
            outFormat='raw'
          )			
          
        }
        
        bgSum <- sum(predBg)
        
        # predict to raster
      } else {
        
        thisRaster <- genericWriteRaster(
          predictorStack=predStack,
          theModel=trialModel,
          fileName=paste0(scratchDir, '/', thisRand, '/mxTuning_beta_', round(100 * thisBeta)),
          predArgs=list(
            maxent=list(
              outputFormat='raw',
              predictFx=predictFx
            )
          ),
          na.rm=T,
          format='raster',
          overwrite=TRUE
        )
        
        bgSum <- cellStats(thisRaster, 'sum')
        
        rm(thisRaster); gc()
        
      }
      
      
      # waste time... averts a write error, somehow
      Sys.sleep(1)
      
      # delete temp dir
      unlink(paste0(scratchDir, '/', thisRand), recursive=TRUE)
      
      ## calculate log likelihood
      
      logLik <- sum(log(predPres / mapSum), na.rm=TRUE)
      
      ## calculate number of parameters
      K <- 0
      
      for (thisLambda in Lambdas <- trialModel@lambdas) { # for each line in lambda object
        
        commaPos <- gregexpr( text=thisLambda, pattern=',') # get location of commas
        
        if (length(commaPos[[1]]) > 1) { # if there is >1 comma in this line (this is not a parameter line)
          
          paramValue <- as.numeric( substr(x=thisLambda, start=commaPos[[1]][1]+1, stop=commaPos[[1]][2]-1 ) ) # convert string between first two commas to numeric
          
          if (paramValue !=0) K <- K + 1 # increment number of parameters
          
        } # if there is >1 comma in this line
        
      }
      
      # AICc
      AICc <- -2 * logLik + 2 * K + (2 * K * (K + 1)) / ( sum(presentBg) - K - 1)
      
      # remember
      thisAicFrame <- data.frame(
        beta=thisBeta,
        n=sum(presentBg),
        logLik=logLik,
        K=K,
        AICc=AICc
      )
      
      return(thisAicFrame)
      
    }
    
    # initiate cluster
    beginCluster(threads)
    thisCluster <- getCluster()
    # on.exit(returnCluster())
    
    # get all nodes going
    for (i in 1:threads) {
      
      # generate random number for raster name
      thisRand <- round(runif(1) * 10^12)
      
      sendCall(thisCluster[[i]], workerFunction, args=list(trainData=trainData, presentBg=presentBg, thisBeta=betaToTest[i], params=params, allPresences=allPresences, predStack=predStack, scratchDir=scratchDir, thisRand=thisRand, predictMax=predictMax))
      
    }
    
    for (i in 1:length(betaToTest)) {
      
      # receive results from a node
      receivedData <- recvOneData(thisCluster)
      
      # cat('Received i: ', i, '\n'); flush.console()
      # cat('Received beta: ', betaToTest[i], '\n\n'); flush.console()
      
      # error?
      if (!receivedData$value$success) stop('\nCluster error in maxentAic function.')
      
      if (exists('aicFrame', inherits=F)) aicFrame <- rbind(aicFrame, receivedData$value$value)
      if (!exists('aicFrame', inherits=F)) aicFrame <- receivedData$value$value
      
      # need to send more data?
      if ((threads + i) <= length(betaToTest)) {
        
        # generate random number for raster name
        thisRand <- round(runif(1) * 10^12)
        
        sendCall(thisCluster[[receivedData$node]], workerFunction, args=list(trainData=trainData, presentBg=presentBg, thisBeta=betaToTest[threads + i], params=params, allPresences=allPresences, predStack=predStack, scratchDir=scratchDir, thisRand=thisRand, predictMax=predictMax))
        
      }
      
    }
    
    endCluster()	
    
  } # end multi-thread
  
  # remove models with more parameters than data points that have more than 0 parameters
  aicFrame <- subset(aicFrame, aicFrame$n >= aicFrame$K & aicFrame$K > 0) 
  aicFrame <- aicFrame[order(aicFrame$K), ] # do this to ensure if any AICc are the same across K that one with the smaller K gets chosen
  aicFrame <- aicFrame[order(aicFrame$beta, decreasing=TRUE), ] # do this to ensure if any AICc are the same across beta that one with the larger beta gets chosen
  aicFrame <- aicFrame[order(aicFrame$AICc), ]
  
  aicFrame$deltaAIcc <- aicFrame$AICc - min(aicFrame$AICc)
  aicFrame$relLike <- exp(-0.5 * aicFrame$deltaAIcc)
  aicFrame$aicWeight <- aicFrame$relLike / sum(aicFrame$relLike)
  
  if (verbose) {
    
    cat('\n')
    print(aicFrame)
    flush.console()
    
  }
  
  # if user wants best model returned
  if (returnModel) {
    
    if (jack) params[which(params=='jackknife=false')] <- 'jackknife=true'
    if (resp) params[which(params=='responsecurves=false')] <- 'responsecurves=true'
    
    # train model
    if (nrow(aicFrame) > 0) {
      
      model <- maxent(
        x=trainData,
        p=as.vector(presentBg),
        removeDuplicates=FALSE,
        path=paste0(scratchDir, '/', thisRand),
        args=c(
          paste0('betamultiplier=', ifelse(all(is.infinite(aicFrame$AICc)), aicFrame$beta[which.max(aicFrame$logLik)], aicFrame$beta[which.min(aicFrame$AICc)])),
          params
        )
      )
      
      if (!resp) unlink(paste0(scratchDir, '/', thisRand), recursive=TRUE)
      
    } else if (anyway) {
      
      warning('Returning model with beta = 1 even though no model had fewer coefficients than predictors.', immediate.=TRUE)
      
      model <- maxent(
        x=trainData,
        p=as.vector(presentBg),
        removeDuplicates=FALSE,
        path=paste0(scratchDir, '/', thisRand),
        args=c('betamultiplier=1', params)
      )
      
      if (!resp) unlink(paste0(scratchDir, '/', thisRand), recursive=TRUE)
      
    } else {
      
      warning('No models had fewer coefficients than predictors. No model returned.', immediate.=TRUE)
      model <- 'No MAXENT model had number of parameters < number of training presences.'
      
    }
    
  }
  
  # return stuff
  if (returnModel & !returnAicFrame) {
    return(model)
  } else if (!returnModel & returnAicFrame) {
    return(aicFrame)
  } else {
    return(list(aicFrame=aicFrame, model=model))
  }
  
}