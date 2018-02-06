

rm(list = ls())

# spokePlot function
spokePlot <- function(pos, neg=NULL, ontop='pos', labels=NULL, shrink=1, labelOffset=1.02, nudge=1, pch=16,
  cexPoints=1, cexLabel=1, lwdPos=1, lwdNeg=1, ltyPos='solid', ltyNeg='dashed',
  colPos='black', colNeg='black', colPoints='black', colLabel='black', ...) {
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

# maxentAic function
maxentAic <- function(trainData, presentBg, betaToTest=betaToTest, predStack=NULL, scratchDir=NULL,
  params=c('linear=true', 'quadratic=true', 'product=true', 'threshold=true', 'hinge=true', 'addsamplestobackground=false', 'jackknife=true', 'responsecurves=false'),
  predictFx='default',threads=1, returnModel=TRUE, returnAicFrame=TRUE, anyway=TRUE, verbose=FALSE) {
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




