require(foreach)
require(ref)
require(multicore)
require(doMC)


#' Sort a ratios matrix based on the times given in the colMaps
#' 
#' @param ratios  The ratios matrix
#' @param colMaps  The colMap object used by inferelator
#' @export
#' @usage ratios <- sortRatiosByTime(ratios, colMap)
sortRatiosByTime <- function(ratios, colMap) {
	expList<-colMap$time
	names(expList)<-rownames(colMap)
	expList<-sort(expList)
	sortIdx<-match(names(expList),colnames(ratios))
	sortIdx<-sortIdx[!is.na(sortIdx)]
	return(ratios[,sortIdx])
}


#' OBSOLETE Given the genes in a BiCluster and the ratios matrix, calculate the conditions that belong
#' 
#' @param geneNames  The names of the relevant genes
#' @param ratios  The ratios matrix
#' @param numSamples  The number of samples (DEFAULT: 10,000)
#' @param pVal  The pVal cutoff (DEFAULT: 0.05)
#' @export
#' @usage varCutoff <- getCondMembers(geneNames, ratios, numSamples = 10000, pVal = 0.05)
getCondMembers <- function(geneNames, ratios, numSamples = 10000, pVal = 0.05) {
	cutOffs <- getVarianceCutoff(ratios, length(geneNames),numSamples = 10000, pVal = 0.05)
	condVars <- apply(ratios[rownames(ratios) %in% geneNames,],2,var,na.rm=T)
	
	exp2include <- names(condVars)[condVars < cutOffs]
	return(exp2include)
}


#' OBSOLETE Use the distribution of variances to find the variance below a pValue cutoff
#' 
#' @param ratios  The ratios matrix
#' @param n  The number of genes
#' @param numSamples  The number of samples (DEFAULT: 10,000)
#' @param pVal  The pVal cutoff (DEFAULT: 0.05)
#' @export
#' @usage varCutoff <- getVarianceCutoff(ratios, n,numSamples = 10000, pVal = 0.05)
getVarianceCutoff <- function(ratios, n, numSamples = 10000, pVal = 0.05) {
	varDist <- getVarianceDist(ratios, n,numSamples)
	cutOffs <- apply(varDist,2, function(x) qnorm(pVal,mean(x),sd(x)) )
	return(cutOffs)
}

#' Given a ratios matrix and a number of genes, figure out the expected distribution of variances
#'   Old version that will calculate number of samples using 31 samples
#' 
#' @param ratios  A refdata pointing to the ratios matrix
#' @param n  The number of genes
#' @param numSamples  The number of samples (DEFAULT: -1, i.e. auto calculate 99% confidence within 1/30s of maximum sampled standard deviation)
#' @export
#' @usage varDist <- getVarianceDist(ratios, n,numSamples = 10000)
getVarianceDist.bk <- function(ratios, n, numSamples = -1) {
	ratios <- derefdata(ratios)

	#Optionally, determine the number of samples you should take with boot strapping
	if (numSamples <= 0) {
		sampleN <- 31 #Based on Statistics by Trioloa, this should be good assuming a normal distribution
		vars <- matrix(0,nrow=sampleN,ncol=ncol(ratios))
		means <- matrix(0,nrow=sampleN,ncol=ncol(ratios))
		for (i in 1:sampleN) {
			vars[i,] <- apply(ratios[sample.int(nrow(ratios),n),],2,var,na.rm=T)
			#means[i,] <- apply(ratios[sample.int(nrow(ratios),n),],2,mean,na.rm=T)
		}

		alpha <- 0.01 #99% confidence interval
		
		#n = [(Za/2 * sd) / E]
		#E = margin of error, set to 1/30 * variance
		#  Why 30?  Because in test case that would set the error range to .1 (i.e. sd ~= 3)
		#  Why not just use E = .1?  Data can be scaled differently so .1 may not always be best
		#n = [(Za/2 * sd) / (sd^2 / 30)]
		#n = [(30*Za/2)/sd]
		numSamples <- round(( (30*qnorm(alpha/2))/sqrt(max(vars)) )^2)
	}

	vars <- matrix(0,nrow=numSamples,ncol=ncol(ratios))
	for (i in 1:numSamples) {
		vars[i,] <- apply(ratios[sample.int(nrow(ratios),n),],2,var,na.rm=T)
	}
	colnames(vars)<-colnames(ratios)
	return(vars)
}

#' Given a ratios matrix and a number of genes, figure out the expected distribution of variances
#'   Will sample background until the mean and sd converge or the operation times out
#'   Non-parallel due to speed
#' 
#' @param ratios  A refdata pointing to the ratios matrix
#' @param n  The number of genes
#' @param tolerance  The fraction tolance to use as a stopping condition (DEFAULT: 0.01)
#' @param maxTime  The approximate maximum time to run in seconds (DEFAULT: 600)
#' @param chunkSize  The number of samples to add between test (DEFAULT: 200)
#' @param verbose  Set to false to supress output (DEFAULT: F)
#' @param numSamples  Does nothing.  Temp included for backwards compatability (DEFAULT: NULL)
#' @export
#' @usage varDist <- getVarianceMeanSD(ratios, n, tolerance = 0.05 ,maxTime=600, chunkSize=200, verbose=F, numSamples=NULL)
getVarianceMeanSD <- function(ratios, n, tolerance = 0.01 ,maxTime=600, chunkSize=200, verbose=F, numSamples=NULL) {
	ratios <- derefdata(ratios)

	keepRunning <- rep(T,ncol(ratios))
	prevMeans <- rep(0,ncol(ratios))
	prevSds <- rep(0.001,ncol(ratios))
	curMeans <- rep(Inf,ncol(ratios))
	curSds <- rep(Inf,ncol(ratios))
	keepRunning <- rep(T,ncol(ratios))
	names(keepRunning) <- names(curMeans) <- names(curSds) <- names(prevMeans) <- names(prevSds) <- colnames(ratios)
	
	rawScores <- list()
	
	startTime <- proc.time()[3]
	curRep<-1
	
	if (verbose) { cat('\nSample', n, "genes from", ncol(ratios), "conditions") }
	while(any(keepRunning) & ((proc.time()[3]-startTime) < maxTime)) {
		if (verbose) { cat('\nRep = ', curRep," ( ",round(proc.time()[3]-startTime),"s ) ",sep="") }
		prevMeans <- curMeans
		prevSds <- curSds
		newScores <- foreach (i = which(keepRunning)) %do% {
			#if (verbose) { cat(i,",",sep="") }
			curReg <- names(keepRunning)[i]
			curScores <- rep(0,chunkSize)
			for (j in 1:chunkSize) {
				curScores[j] <-var(ratios[sample.int(nrow(ratios),n),curReg],na.rm=T)
			}
			curScores
		}

		#browser()

		for (ctr in 1:sum(keepRunning)) {
			i<-which(keepRunning)[ctr]
			curReg <- names(keepRunning)[i]
			rawScores[[curReg]] <- c(rawScores[[curReg]],newScores[[ctr]])
			curMeans[i] <- mean(rawScores[[curReg]],na.rm=T)
			curSds[i] <- sd(rawScores[[curReg]],na.rm=T)
		}
		
		#browser()
		keepRunning <- (abs(curMeans-prevMeans) >= tolerance*abs(prevMeans)) | (abs(curSds-prevSds) >= tolerance*prevSds)
		keepRunning[is.nan(keepRunning)] <- TRUE
		keepRunning[is.na(keepRunning)] <- TRUE
		
		curRep <- curRep+1
	}
	if(verbose) {cat('\n')}
	
	return(list(means=curMeans,sds=curSds))
}


#' Return the means and SDs for the variances for each experimental condition
#'   Uses: getVarianceDist.bk
#' 
#' @param ratios  A refdata pointing to the ratios matrix
#' @param n  The number of genes
#' @param numSamples  The number of samples (DEFAULT: -1 i.e. autodetect)
#' @export
#' @usage varCutoff <- getVarianceMeanSD(ratios, n,numSamples = 10000)
getVarianceMeanSD.bk <- function(ratios, n, numSamples = -1) {
	varDist <- getVarianceDist(ratios, n,numSamples)
	means <- apply(varDist,2, mean )
	sds <- apply(varDist,2, sd )
	return(list(means=means,sds=sds))
}

#' Return the mean and SD for the variances of all experimental conditions
#' 
#' @param ratios  A a list of ratios matrices
#' @param n  The number of genes
#' @param numSamples  The number of samples (DEFAULT: -1 i.e. autodetect)
#' @export
#' @usage varCutoff <- getVarianceMeanSD.all(ratios, n,numSamples = 10000)
getVarianceMeanSD.all <- function(ratios, n, numSamples = -1) {
	varDist <- NULL
	for (i in 1:length(ratios)) {
		varDist<-c(varDist,as.numeric(getVarianceDist(refdata(e$ratios[[i]]), n,numSamples)))
	}
	return(data.frame(means=mean(varDist),sds=sd(varDist)))
}


#' OBSOLETE
#' Return a dataframe containing the means and standard deviations for all ratios matrixes in list ratios
#' 
#' @param ratios  A list of ratios matrixes
#' @param n  The number of genes
#' @param numSamples  The number of samples (DEFAULT: -1, i.e. autodetect)
#' @param all  Set to false to get a different means and sd for each experimental condition (DEFAULT: F)
#' @export
#' @usage varCutoff <- getVarianceMeanSD.DF.bk(ratios, n,numSamples = 10000, all = F)
getVarianceMeanSD.DF.bk <- function(ratios, n, numSamples = -1, all = F) {
    expNames <- as.character(unlist(lapply(ratios, colnames)))
    if (n == 1) {
        sds <- means <- rep(0,length(expNames))
    } else {
    	if (all == T) {
    		means.sds<-getVarianceMeanSD(ratios, n, numSamples)
    		means<-rep(means.sds$means, sum(sapply(ratios,ncol)))
        	sds<-rep(means.sds$sds, sum(sapply(ratios,ncol)))
        } else {
        	means.sds<-sapply(ratios,function(x) {getVarianceMeanSD(refdata(x), n, numSamples = numSamples)})
        	means<-unlist(means.sds['means',],use.names=F)
        	sds<-unlist(means.sds['sds',],use.names=F)
        }
    }
    names(means) <- names(sds) <- expNames
    return(data.frame(means=means,sds=sds))
}

#' Return a dataframe containing the means and standard deviations for all ratios matrixes in list ratios
#' 
#' @param ratios  A list of ratios matrixes
#' @param n  The number of genes
#' @export
#' @usage varCutoff <- getVarianceMeanSD(ratios, n)
getVarianceMeanSD.DF <- function(ratios, n) {
    expNames <- as.character(unlist(lapply(ratios, colnames)))
    if (n <= 1 ) {
        sds <- means <- rep(0,length(expNames))
    } else {
    	means.sds<-sapply(ratios,function(x) {getVarianceMeanSD(refdata(x), n)})
       	means<-unlist(means.sds['means',],use.names=F)
       	sds<-unlist(means.sds['sds',],use.names=F)
    }
    names(means) <- names(sds) <- expNames
    df <- data.frame(means=means,sds=sds)
    rownames(df) <- expNames
    return(df)
}

#' UNTESTED Return a pVal or vector of pVals for cluster K.
#' Intended to replace cluster.resid
#' Assumes "rows","rats","cols" `as part of "..."
#' 
#' @param k  The cluster number to calculate the pValue for
#' @param rats.inds  If not set to combined, will return one pVal for each member of "ratios" (DEFAULT: COMBINED)
#' @param means.sds  a list containing one data.frame for each number of genes.  Each data frame has means and sds. (DEFAULT: list())
#' @export
#' @usage resids <- cluster.ratPval( k, rats.inds="COMBINED", means.sds=means.sds ) 
cluster.ratPval <- function( k, rats.inds="COMBINED", means.sds=list(), clusterStack = get("clusterStack"), ... ) {
	clust <- clusterStack[[k]]
	numGenesInClusters<- clust$nrows
	
	#Add to means.sds as necessary
	numGeneList<-numGenesInClusters[! numGenesInClusters %in% names(means.sds)]
	lGeneList<-length(numGeneList)
	if (lGeneList > 0 ) {
		numGeneLists<-split(1:lGeneList,cut(1:lGeneList,floor(lGeneList/e$parallel.cores)))
		cat("Calculating variance means and sd's for",lGeneList,"gene counts\n")
		for (geneListNums in numGeneLists) { #Embed in for loop to have a status monitor
			curNumGenes<-numGeneList[geneListNums]
			cat(curNumGenes,'',sep=",")
			flush.console()
			new.means.sds<- foreach (n=curNumGenes, .inorder=TRUE) %do% {
				getVarianceMeanSD.DF(rats, n) 
			} #for (n in numGeneList)
			names(new.means.sds)<-as.character(curNumGenes)
			means.sds<-c(means.sds,new.means.sds)
		}
		cat('\n')
	}
	
	#Get the average pValue
	curVarList <- means.sds[[as.character(numGenesInClusters)]]
	
	#Get the pValues for each experiment given the number of genes
	rw<-get( "row.weights" )
	pVals <- rep(NA,length(rw))
	names(pVals)<-names(rw)
	for (i in 1:length(pVals) ) {
		rats <- get ("ratios")
		curCols <- clust$cols[clust$cols %in% colnames(rats[[i]])]
		relRats<- rats[[i]][rownames(rats[[i]]) %in% clust$rows,curCols]
		vars <- apply(relRats,2,var,na.rm=T)
		curPs<-NULL
		for (x.idx in 1:length(vars)) {
			x<-vars[x.idx]
			curPs<-c(curPs, pnorm(x,curVarList[names(x),'means'],curVarList[names(x),'sds']) )
		}
		pVals[i]<-mean(curPs)
	}
	
	if ( rats.inds[ 1 ] == "COMBINED" ) pVals <- weighted.mean( pVals, row.weights[ inds ], na.rm=T )

	return(pVals)
}

#' Return a list of pVals for a list of genes under all conditions
#' 
#' @param e  The cMonkey environment.  Used for the ratios matrix.  Will also accept a list containing a list of ratios matrices
#' @param geneList  one data.frame for each number of genes.
#' @param mean.sd  A single elements of means.sds.  A DF with means and sds, one for each experimental condition.  
#' @export
#' @usage pValList <- getRatPvals(e, geneList, mean.sd=means.sds[["6"]])
getRatPvals <- function(e, geneList, mean.sd) {
	col.pVals<-NULL
	for (ratIdx in 1:length(e$ratios)) {
		curExps <- colnames(e$ratios[[ratIdx]])
		relRats<- e$ratios[[ratIdx]] [rownames(e$ratios[[ratIdx]]) %in% geneList,,drop=F]
		if (is.array(relRats)) { 
			vars <- apply(relRats,2,var,na.rm=T)
		} else {
			vars <- rep(0,length(relRats))
			names(vars)<-names(relRats)
		}

		#browser()
		pVals<-vars
		for (x.idx in 1:length(vars)) {
			x<-vars[x.idx]
			if ( any( !(names(x) %in% rownames(mean.sd)) ) ) { 
				cat('Warning:  Some experiments not in mean.sd\n')
				missingExp <- sprintf('%s,',names(x)[!(names(x) %in% rownames(mean.sd))])
				cat('\t',missingExp,'\n')
			}
			pVals[x.idx]<-pnorm(x,mean.sd[names(x),'means'],mean.sd[names(x),'sds']) 
		}

		col.pVals <- c(pVals,col.pVals)
	}
	col.pVals
}


#' Given a cMonkey environment, build a new clusterStack with different cols & pValues based on variance
#' Returns "newClusterStack" and "means.sds".  "means.sds" may be used in subsequent calls to avoid recomputation
#' 
#' @param e  The cMonkey environment
#' @param means.sds  a list containing one data.frame for each number of genes.  Each data frame has means and sds. (DEFAULT: generate all)
#' @param numSamples  The number of times to sample the background distribution to determine the pValues (DEFAULT: -1, i.e. autodetect)
#' @param pValCut  The pValue at which to cut-off bicluster inclusion.  (DEFAULT: 0.05)
#' @param bonferroni  Set to TRUE to turn bonferroni correction on.  Prelimiary tests show this to be too strict. (DEFAULT: FALSE)
#' @param aveNumCond  The average number of conditions to include in each cluster.  If set, will ignor pValCut   (DEFAULT: NULL)
#' @export
#' @usage clusterStack <- resplitClusters(e, means.sds=list(), numSamples = 10000, bonferroni = F, aveNumCond=NULL))
resplitClusters <- function(e, means.sds=list(), numSamples = -1, pValCut = 0.05, bonferroni = F, aveNumCond=NULL, all = T) {

	#Calculate the background distribution
	numGenesInClusters<- unique(sort(sapply(e$clusterStack,function(x) {x$nrows} )))

	#Load means.sds as necessary
	numGeneList<-numGenesInClusters[! numGenesInClusters %in% names(means.sds)]
	lGeneList<-length(numGeneList)
	if (lGeneList > 0 ) {
		numGeneLists<-split(1:lGeneList,cut(1:lGeneList,max(floor(lGeneList/e$parallel.cores),2)))
		cat("Calculating variance means and sd's for",lGeneList,"gene counts\n")
		for (geneListNums in numGeneLists) { #Embed in for loop to have a status monitor
			curNumGenes<-numGeneList[geneListNums]
			cat(curNumGenes,'',sep=",")
			flush.console()
			new.means.sds<- foreach (n=curNumGenes, .inorder=TRUE) %dopar% {
				getVarianceMeanSD.DF(e$ratios, n)
			} #for (n in numGeneList)
			names(new.means.sds)<-as.character(curNumGenes)
			means.sds<-c(means.sds,new.means.sds)
		}
		cat('\n')
	}
	
	#Dynamically calculate the pValue cutoffs for each numGenesInCluster
	pValCuts <- rep (pValCut,length(numGenesInClusters))
	names(pValCuts) <- as.character(numGenesInClusters)
	
	if (! is.null(aveNumCond)) {
	
		for ( i in 1:length(numGenesInClusters) ) {
			numGenes <- numGenesInClusters[i]
			curIdxs <- sapply(e$clusterStack,function(x) {x$nrows} ) == numGenes
			if (any(curIdxs) & (numGenes > 1)){
				mean.sd <- means.sds[[as.character(numGenes)]]
				varLists<-lapply(which(curIdxs), function(x) {getRatPvals(e, e$clusterStack[[x]]$rows, mean.sd)})
				pValCuts[as.character(numGenes)]<-as.numeric(sort(c(varLists,recursive=T))[aveNumCond*sum(curIdxs)])
			} else {
				pValCuts[as.character(numGenes)] <- 1
			}
		}
	} #if (pValCut <= 0)
	
	#Build the new clusterStack
	newClustStack <- foreach (clust=e$clusterStack, .inorder=T) %do% {
		if (clust$nrows > 1) {
			pValCut <- pValCuts[as.character(clust$nrows)]
			
			#Get the pValues for each experiment given the number of genes
			col.pVals <- getRatPvals(e, clust$rows, means.sds[[as.character(clust$nrows)]])

			#Bonferroni correction
			if (bonferroni == T) {pValCut <- pValCut / sum(sapply(e$ratios,ncol))}

			newClust<-clust
			newClust$cols <- names(col.pVals)[col.pVals < pValCut]
			if (length(newClust$cols) <= 2) {  newClust$cols <- names(col.pVals)[order(col.pVals)[1:2]] }
			newClust$ncols <- length(newClust$cols)

			#Update the residuals to be the mean pValue of the included clusters
			for (idx in 1:length(newClust$resid)) {
				if (!is.null(names(clust$resid[idx]))) {
					relCols <- newClust$cols[newClust$cols %in% colnames(e$ratios[[names(clust$resid[idx])]])]
				} else {				
					relCols <- newClust$cols[newClust$cols %in% colnames(e$ratios[[1]])]
				}
				newClust$resid[idx] <- mean(col.pVals[relCols])
			}
		} else {
			newClust<-clust
		}
		newClust
	}
	
	return(list(newClusterStack=newClustStack, means.sds=means.sds))
}

#' Return the column scores for a given cluster k
#' Default is to normalize (diff)^2 by mean expression level, similar to "index of dispersion"
#'    http://en.wikipedia.org/wiki/Index_of_dispersion
#' 
#' @param k  The cluster number or a vector with the row names.
#' @param for.cols  The column names to include. (DEFAULT: "all")
#' @param ratios  The rations matrix.  (DEFAULT: ratios[[1]])
#' @param norm.method  The normalization method. (DEFAULT: c("mean","all.colVars","none")[1])  
#' @param scoring.method  The scoring method.  res.like is resildual like, var.p is variance pValue  (DEFAULT: c("var.p","res.like")[1])
#' @param ...    
#' @export
#' @usage pValList <- get.col.scores(k, for.cols="all", ratios=ratios[[ 1 ]],norm.method=c("mean","all.colVars","none")[1],scoring.method=c("var.p","res.like")[1])
get.col.scores <- function( k, for.cols="all", ratios=ratios[[ 1 ]], 
                           norm.method=c("mean","all.colVars","none")[1], 
                           scoring.method=c("var.p","res.like")[1],... ) {
	## Compute scores for ALL cols (over just the rows IN the cluster)
	if ( length( k ) <= 0 ) return( NULL )
	if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k ) 
	else rows <- k
	if ( for.cols[ 1 ] == "all" ) for.cols <- colnames( ratios ) ##attr( ratios, "cnames" )
	rows <- rows[ rows %in% rownames( ratios ) ]
	if ( length( rows ) <= 1 ) return( rep( NA, length( for.cols ) ) )

	rats <- ratios[ rows, for.cols, drop=F ]
	
	if (scoring.method != "var.p") { "res.like"
		row.weights <- if ( exists( "get.row.weights" ) ) get.row.weights( rows, cols, ratios ) else NA

		if ( is.na( row.weights[ 1 ] ) ) { ## Default
			rats.mn <- matrix( colMeans( rats, na.rm=T ), nrow=nrow( rats ), ncol=ncol( rats ), byrow=T )
		} else { ## Custom row weights
			rats.mn <- matrix( apply( rats, 2, weighted.mean, w=row.weights[ rows ], na.rm=T ), ncol=ncol( rats ), byrow=T )
		}

		rats[,] <- ( rats[,] - rats.mn )^2 ## abs( Multiplying is faster than squaring
		rats <- colMeans( rats, na.rm=T )

		var.norm <- 0.99
		if ( norm.method == "all.colVars" ) {
			all.colVars <- attr( ratios, "all.colVars" )
			if ( ! is.null( all.colVars ) ) var.norm <- all.colVars[ for.cols ]
		} else if ( norm.method == "mean" ) {
			var.norm <- abs( rats.mn[ 1, ] ) ##0.99 ## Use the mean expr. level (higher expressed expected to have higher noise)
		}

		##col.weights <- get.col.weights( rows, cols )
		##if ( is.na( col.weights ) )
		rats <- rats / ( var.norm + 0.01 ) ## default
		##!else rats <- colMeans( rats, na.rm=T ) / ( var.norm * col.weights[ cols ] + 0.01 ) ## customized col. weights

		##return( log( rats + 1e-99 ) )
	} else { #var.p,  Use the bicluster.by.variance code to get a pValue
		numGenes <- nrow(rats)
		#means.sds must be pulled from the environment.  It is a class variable
		#browser()
		if (is.null(means.sds[[as.character(numGenes)]])) {
			cat('\tRecalculating means and sds for',numGenes,'\n')
			means.sds[[as.character(numGenes)]] <- getVarianceMeanSD(refdata(get.cluster.matrix()), numGenes) 
		}
		mean.sd <- means.sds[[as.character(numGenes)]]
		
		ratList <- list()
		ratList[["ratios"]] <- list(rats)
		pVals <- getRatPvals(ratList, rows, mean.sd)  #env is usual first parameter.  Fake it with ratList
		
		rats <- log(pVals)
	}
	rats
}

#' Make sure that env$means.sds contains all necessary background distributions
#' 
#' @param env  The cMonkey environment.
#' @export
#' @usage env <- update.means.sds(env)
update.means.sds <- function( env ) {
	mc <- env$get.parallel()
	means.sds <- as.list(env$means.sds)
	if (is.null(env$scoring.method) || env$scoring.method == "var.p") {
		#Some conditions tested may have less than the number of genes in the cluster due to missing data
		maxNumGenes <- max ( as.numeric(sapply(env$clusterStack,function(x) {x$nrows} )), colSums(env$row.memb))
		maxNumGenes <- maxNumGenes + 1 #Include + one gene in case of growth
		numGenesInClusters <- 1:maxNumGenes
  	
  		#Load means.sds as necessary
		numGeneList <- numGenesInClusters[! numGenesInClusters %in% names(means.sds)]
		lGeneList <- length(numGeneList)
		if (lGeneList > 0) {
			if (mc$par == FALSE) {cores <- 1} else {cores <- mc$par}
			numGroups <- floor(lGeneList/cores)
			numGeneLists <- list(1:lGeneList)
			if (numGroups > 2) { numGeneLists <- split(numGeneLists[[1]],cut(numGeneLists[[1]],numGroups)) }
			cat("Calculating variance means and sd's for",lGeneList,"gene counts\n")
			for (geneListNums in numGeneLists) { #Embed in for loop to have a status monitor
				curNumGenes <- numGeneList[geneListNums]
				cat(curNumGenes,'',sep=",")
				flush.console()
				#browser()
				new.means.sds <- mc$apply( curNumGenes, function( n ) {
					getVarianceMeanSD.DF(env$ratios, n)
				} ) #for (n in numGeneList)
				names(new.means.sds) <- as.character(curNumGenes)
				means.sds <- c(means.sds,new.means.sds)
			}
		}
		cat('\n')
	} #if (scoring.method == "var.p") {
	if ( ! is.null( env ) ) assign( "means.sds", means.sds, envir=env )
	invisible(env)
}

#' Seed the biclusters
#'   Modified 03-15-11 so that get.col.scores will use scoring.method="res.like"
#' 
#' @param k.clust  The number of clusters to create.
#' @param seed.method  . (DEFAULT: "rnd")
#' @param col.method  . (DEFAULT: "rnd")
#' @export
#' @usage env <- seed.clusters(env)
seed.clusters <- function( k.clust, seed.method="rnd", col.method="rnd" ) {
  ## Allow it to be overridden by a custom function if desired (e.g. to seed from a gene list -- no that's a bad
  ## example -- there's the "list=" option below)
  if ( seed.method == "custom" && exists( "seed.clusters.custom" ) )
    return( seed.clusters.custom( k.clust, col.method ) )
  if ( substr( seed.method, 1, 3 ) == "net" && length( networks ) <= 0 ) {
    cat( "Seed method is", seed.method, ", but no networks -- using 'kmeans' instead.\n" )
    seed.method <- "kmeans"
  }
  if ( seed.method == "rnd" ) { ## RANDOM, ALL GENES ASSIGNED TO k.clust CLUSTERS
    row.membership <- t( sapply( 1:attr( ratios, "nrow" ),
                                function( i ) sample( 1:k.clust, n.clust.per.row[ 1 ],
                                                     replace=n.clust.per.row[ 1 ] > attr( ratios, "nrow" ) ) ) )
  } else if ( substr( seed.method, 1, 5 ) == "list=" ) { ## File with sets of genes - one set of genes per line
    row.membership <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=n.clust.per.row[ 1 ] ) ## separated by space,tab,or comma
    rownames( row.membership ) <- attr( ratios, "rnames" )
    fname <- strsplit( seed.method, "=" )[[ 1 ]][ 2 ]
    if ( exists( fname ) ) lists <- get( fname ) ## List of vectors already exists in memory
    else if ( file.exists( fname ) ) lists <- strsplit( readLines( fname ), split="[ \t,]", perl=T )
    for ( k in 1:min( c( k.clust, length( lists ) ) ) ) {
      probes <- unlist( lapply( get.synonyms( lists[[ k ]] ), function( i ) i %in% rownames( row.membership ) ) )
      row.membership[ probes[ row.membership[ probes, 1 ] == 0 ], 1 ] <- k
      row.membership[ probes[ row.membership[ probes, 1 ] != 0 ], 2 ] <- k
    }
    if ( length( lists ) < k.clust ) { ## Fill in additional clusters randomly
      for ( k in ( length( lists ) + 1 ):k.clust ) {
        rnames <- attr( ratios, "rnames" )[ ! attr( ratios, "rnames" ) %in% unlist( lists ) ]
        rows <- sample( rnames, 5 )
        row.membership[ rows[ row.membership[ rows, 1 ] == 0 ], 1 ] <- k
        row.membership[ rows[ row.membership[ rows, 1 ] != 0 & row.membership[ rows, 2 ] == 0 ], 2 ] <- k
      }
    }
  } else if ( substr( seed.method, 1, 4 ) == "rnd=" ) { ## RANDOM SAMPLED N per cluster
    n.samp <- as.integer( strsplit( seed.method, "=" )[[ 1 ]][ 2 ] )
    row.membership <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=n.clust.per.row[ 1 ] )
    rownames( row.membership ) <- attr( ratios, "rnames" )
    for ( i in 1:n.clust.per.row ) {
      sampled <- rep( FALSE, attr( ratios, "nrow" ) ); names( sampled ) <- attr( ratios, "rnames" )
      for ( k in 1:k.clust ) {
        g <- sample( attr( ratios, "rnames" )[ ! sampled ], n.samp )
        row.membership[ g, 1 ] <- k
        sampled[ g ] <- TRUE
      }
    }
  } else if ( seed.method == "kmeans" ) { ## Kmeans seeded
    if ( ! exists( "ratios" ) ) stop( "kmeans seed method but no ratios" )
    tmp.rat <- get.cluster.matrix() ##ratios;
    tmp.rat[ is.na( tmp.rat ) ] <- 0
    ##cat(dim(tmp.rat),k.clust,"\n")
    km <- kmeans
    row.membership <- km( tmp.rat, centers=k.clust, iter.max=20, nstart=2 )$cluster
    names( row.membership ) <- attr( ratios, "rnames" )
    if ( n.clust.per.row[ 1 ] > 1 ) row.membership <-
      cbind( row.membership, matrix( rep( 0, attr( ratios, "nrow" ) * ( n.clust.per.row[ 1 ] - 1 ) ),
                                         ncol=n.clust.per.row[ 1 ] - 1 ) )
  }
  
  if ( is.vector( row.membership ) ) row.membership <- t( row.membership ) ## probably n.clust.per.row == 1
  if ( nrow( row.membership ) == 1 ) row.membership <- t( row.membership ) ## probably n.clust.per.row == 1
 
  if ( col.method == "rnd" ) {
    col.membership <- t( sapply( 1:attr( ratios, "ncol" ), function( i )
                                sample( 1:k.clust, n.clust.per.col[ 1 ],
                                       replace=n.clust.per.col[ 1 ] > k.clust ) ) ) ##attr( ratios, "ncol" ) ) ) )
  } else if ( col.method == "best" ) {
    if ( ! exists( "ratios" ) ) stop( "best col seed method but no ratios" )
    all.rats <- get.cluster.matrix()
    attr( all.rats, "all.colVars" ) <- apply( all.rats, 2, var, use="pair", na.rm=T )
    col.scores <- -sapply( 1:k.clust, function( k )
                          if ( sum( row.membership == k, na.rm=T ) <= 0 ) rep( NA, attr( ratios, "ncol" ) ) else 
                          get.col.scores( k=get.rows( k, row.membership ), ratios=all.rats, method="orig" ,scoring.method="res.like") ) 
    col.membership <- t( apply( col.scores, 1, function( i ) order( i, decreasing=T )[ 1:n.clust.per.col[ 1 ] ] ) )
  }

  rownames( row.membership ) <- attr( ratios, "rnames" ) 
  rownames( col.membership ) <- attr( ratios, "cnames" )
  list( row.membership=row.membership, col.membership=col.membership )
}

#' Run one iteration of the cMonkey algorithm
#'   BIG TODO: make it so all big.memory are updated in-place
#'   Updated on 03-15-11 to make sure that all necessary background ditributions are calculated for var.p column scoring
#' 
#' @param env  The cMonkey environment.
#' @param dont.update  Set to T to not update the iteration counter & other iteration related functions. (DEFAULT: F)
#' @export
#' @usage env <- cmonkey.one.iter(env, dont.update=F)
cmonkey.one.iter <- function( env, dont.update=F, ... ) {  
  ##if ( ! exists( "clust.changed", envir=env ) ) env$clust.changed <- rep( FALSE, k.clust )
  if ( ! exists( "row.memb", envir=env ) ) {
    env$row.memb <- t( apply( row.membership[], 1, function( i ) 1:k.clust %in% i ) )
    env$col.memb <- t( apply( col.membership[], 1, function( i ) 1:k.clust %in% i ) )
    env$row.memb <- matrix.reference( env$row.memb, backingfile="row.memb", backingpath=env$cmonkey.filename )
    env$col.memb <- matrix.reference( env$col.memb, backingfile="col.memb", backingpath=env$cmonkey.filename )
  } else {
    env$row.memb[,] <- t( apply( row.membership[], 1, function( i ) 1:k.clust %in% i ) )
    env$col.memb[,] <- t( apply( col.membership[], 1, function( i ) 1:k.clust %in% i ) )
  }

  #browser()
  env <- update.means.sds(env) #SD 03-15-11
  tmp <- get.all.scores( ... )
  env$row.scores <- tmp$r##[,];
  env$mot.scores <- tmp$m; env$net.scores <- tmp$n; env$col.scores <- tmp$c
  env$meme.scores <- tmp$ms
  if ( ! is.null( tmp$cns ) ) env$cluster.net.scores <- tmp$cns

  tmp <- get.combined.scores( quant=T )
  env$r.scores <- tmp$r; env$c.scores <- tmp$c; env$n.scores <- tmp$n; env$m.scores <- tmp$m
  if ( ! is.null( env$row.scores ) ) attr( env$row.scores, "changed" ) <- FALSE
  if ( ! is.null( env$col.scores ) ) attr( env$col.scores, "changed" ) <- FALSE
  if ( ! is.null( env$mot.scores ) ) attr( env$mot.scores, "changed" ) <- FALSE
  if ( ! is.null( env$net.scores ) ) attr( env$net.scores, "changed" ) <- FALSE

  if ( length( tmp$scalings ) > 0 ) {
    env$row.scaling[ iter ] <- tmp$scalings[ "row" ]
    env$mot.scaling[ iter ] <- tmp$scalings[ "mot" ]
    env$net.scaling[ iter ] <- tmp$scalings[ "net" ]
  }

  ## Fuzzify scores a bit for stochasticity! (fuzz should be between 0.2 and 0 (decreasing with iter)
  if ( fuzzy.index[ iter ] > 1e-5 ) {
    env$r.scores[,] <- env$r.scores[,] +
      rnorm( length( env$r.scores[,] ), sd=sd( env$r.scores[,][ row.memb[,] == 1 ], na.rm=T ) * fuzzy.index[ iter ] )
    if ( ! is.null( env$c.scores ) ) env$c.scores[,] <- env$c.scores[,] +
      rnorm( length( env$c.scores[,] ), sd=sd( env$c.scores[,][ col.memb[,] == 1 ], na.rm=T ) * fuzzy.index[ iter ] )
  }

  tmp <- get.density.scores( ks=1:k.clust ) ## r.scores, col.scores, 
  env$rr.scores <- tmp$r; env$cc.scores <- tmp$c ##; rm( tmp )

  ## OUTPUT
  if ( iter %in% stats.iters ) { 
    env$clusterStack <- get.clusterStack( ks=1:k.clust ) 
    env$stats <- rbind( stats, get.stats() )
    cat( organism, as.matrix( stats[ nrow( stats ), ] ), '\n' )
  } else {
    cat( sprintf( "==> %04d %.3f %.3f %.3f\n", iter, ##%5d sum( row.memb != old.row.memb, na.rm=T ),
                 mean( row.scores[,][ row.memb[,] == 1 ], na.rm=T ), ##mean( col.scores[ col.memb ], na.rm=T, trim=0.05 ),
                 if ( ! is.null( mot.scores ) ) mean( mot.scores[,][ row.memb[,] == 1 & mot.scores[,] < 0 ], na.rm=T, trim=0.05 )
                 else NA,
                 if ( ! is.null( net.scores ) ) mean( net.scores[,][ row.memb[,] == 1 ##& net.scores < 0
                                                                 ], na.rm=T, trim=0.05 ) else NA ) ) ##, "\n" )
  }
    
  ## NEW - will it work? -- help shrink big clusters, grow small clusters, both in rows and cols
  size.compensation.func.rows <- function( n ) exp( -n / ( attr( ratios, "nrow" ) * n.clust.per.row / k.clust ) )
  size.compensation.func.cols <- function( n ) exp( -n / ( attr( ratios, "ncol" ) * n.clust.per.col / k.clust ) )
  for ( k in 1:k.clust ) {
    tmp <- sum( row.memb[ ,k ] )
    if ( tmp > 0 ) env$rr.scores[ ,k ] <- env$rr.scores[ ,k ] * size.compensation.func.rows( tmp ) 
    else env$rr.scores[ ,k ] <- env$rr.scores[ ,k ] * size.compensation.func.rows( cluster.rows.allowed[ 1 ] )
    if ( ! is.null( env$cc.scores ) ) {
      tmp <- sum( col.memb[ ,k ] )
      if ( tmp > 0 ) env$cc.scores[ ,k ] <- env$cc.scores[ ,k ] * size.compensation.func.cols( tmp ) 
      else env$cc.scores[ ,k ] <- env$cc.scores[ ,k ] * size.compensation.func.cols( attr( ratios, "ncol" ) / 10 )
    }
  }
  
  ## Fuzzify it along the lines of fuzzy c-means clustering
  ##   -- see http://en.wikipedia.org/wiki/Data_clustering#Fuzzy_c-means_clustering
  ## No -- it doesnt affect things - same ordering (and updated memberships are based on ordering)
  ##   but should use these scores to weight the centroids that are selected in the next iteration.
  ##   for this, fuzz should vary between e.g. 10 and 2
  ##rr.scores <- ( rr.scores / sum( rr.scores, na.rm=T ) )^( 2 / ( fuzz - 1 ) )
  ##cc.scores <- ( cc.scores / sum( cc.scores, na.rm=T ) )^( 2 / ( fuzz - 1 ) )

  if ( ! dont.update ) {
  
    ## Make a matrix of m[i,k] = whether row/col i is in cluster k
    if ( exists( "row.membership" ) ) { 
      env$old.row.membership <- row.membership 
      env$old.col.membership <- col.membership 
    }

    tmp <- get.updated.memberships() ## rr.scores, cc.scores )
    env$row.membership <- tmp$r; env$col.membership <- tmp$c
    if ( ! is.null( tmp ) ) { env$row.membership <- tmp$r; env$col.membership <- tmp$c }
    
    
    ## PLOTTING
    if ( ! is.na( plot.iters ) && iter %in% plot.iters ) {
      env$clusterStack <- get.clusterStack( ks=1:k.clust ) 
      try( plotStats( iter, plot.clust=env$favorite.cluster(), new.dev=T ), silent=T ) ## Can be set for your given organism
    }
  
    if ( exists( "cm.func.each.iter" ) ) try( cm.func.each.iter(), silent=T ) ## User-defined func. to run each iteration
    
    ## Allow temp source file to be sourced (e.g. to change a param in the middle of a run, or print or plot or save
    ## some intermediate results). If first line of file is '## QUIET' then this is done quietly.
    if ( any( cm.script.each.iter != "" ) ) {
      for ( f in cm.script.each.iter ) {
        if ( file.exists( f ) && file.info( f )$size > 1 ) {
          tmp <- readLines( f )
          if ( all( substr( tmp, 1, 1 ) == "#" ) ) next ## All commented-out code
          if ( tmp[ 1 ] != "## QUIET" ) cat( "Sourcing the script '", f, "' ...\n", sep="" )
          try( source( f, echo=tmp[ 1 ] != "## QUIET", local=T ), silent=T )
        }
      }
    }
  }

  
  ## Note: with the above code, when the env is saved via save.image(), all ff obj's are "closed" but their
  ##   filestores still exist, so you can "open.ff(x)" each of them after the env is re-loaded,
  ##   and that will reconnect them with their files.
  
  if ( get.parallel()$mc ) {
    if ( getDoParName() == "doMC" ) { ##require( multicore, quietly=T ) ) { ## Clean up any multicore spawned processes (as doc'ed in mclapply help)
      chld <- multicore::children()
      if ( length( chld ) > 0 ) { try( { multicore::kill( chld ); tmp <- multicore::collect( chld ) }, silent=T ) }
    } else if ( getDoParName() == "doSNOW" && "data" %in% ls( pos=foreach:::.foreachGlobals ) ) {
      cl <- get( "data", pos=foreach:::.foreachGlobals ) ## Tricky, eh?
      if ( ! is.null( data ) ) stopCluster( cl )
    }
  }
  
  if ( ! dont.update ) env$iter <- env$iter + 1
  invisible( env )
} #cmonkey.one.iter <- function


#' Change the cluster breaks so that inclusion is based on variance background distribution
#'  This should usually be called at the end of a cMonkey run.
#' 
#' @param e  The cMonkey environment.
#' @param pValCut  The pValue at which to cut-off bicluster inclusion.  (DEFAULT: 0.05)
#' @param colMap  Include to sort the conditions by time.  (DEFAULT: NULL)
#' @export
#' @usage env <- resplitClusters.by.var( e, pValCut = 0.05, colMap=NULL )
resplitClusters.by.var <- function( e, pValCut = 0.05, colMap=NULL ) {
	e$MPV <- getMPV(e) #Save the MPV so that it isn't lost by splitting the cluster

	means.sds<-list()
	if (! is.null(e$means.sds)) { means.sds <- e$means.sds }

	newClusterStack <- resplitClusters(e,means.sds=means.sds)
	e$means.sds <- newClusterStack$means.sds
	
	row.col.membership.from.clusterStack <- function( cs, row.membership, col.membership ) {
		row.memb <- row.membership * 0
		col.memb <- col.membership * 0
		for ( k in 1:length( cs ) ) {
	  		if ( k > ncol( row.memb ) ) row.memb <- cbind( row.memb, rep( 0, nrow( row.memb ) ) )
	  		rows <- cs[[ k ]]$rows; rows <- rows[ ! is.na( rows ) ]
	  		row.memb[ rows, k ] <- k
	  		if ( k > ncol( col.memb ) ) col.memb <- cbind( col.memb, rep( 0, nrow( col.memb ) ) )
	  		cols <- cs[[ k ]]$cols; cols <- cols[ ! is.na( cols ) ]
	  		col.memb[ cols, k ] <- k
	  	}
	  	row.memb <- t( apply( row.memb, 1, function( i ) c( i[ i != 0 ], i[ i == 0 ] ) ) )
	  	row.memb <- row.memb[ ,apply( row.memb, 2, sum ) != 0, drop=F ]
	  	colnames( row.memb ) <- NULL
	  	col.memb <- t( apply( col.memb, 1, function( i ) c( i[ i != 0 ], i[ i == 0 ] ) ) )
	  	col.memb <- col.memb[ ,apply( col.memb, 2, sum ) != 0, drop=F ]
	  	colnames( col.memb ) <- NULL
	  	list( r=row.memb, c=col.memb )
	}
	
	##environment( row.col.membership.from.clusterStack ) <- e
	tmp <- row.col.membership.from.clusterStack( newClusterStack$newClusterStack, e$row.membership, e$col.membership ) ##e$clusterStack )
	e$row.membership <- tmp$r
	e$col.membership <- tmp$c
	e$row.memb <- t( apply( e$row.membership, 1, function( i ) 1:e$k.clust %in% i ) )
	e$col.memb <- t( apply( e$col.membership, 1, function( i ) 1:e$k.clust %in% i ) )
	if (!is.null(colMap)) {
		condNames <- rownames(colMap)[order(colMap$time)]
		condNames <- condNames[condNames %in% rownames(e$col.memb)]
		e$col.memb <- e$col.memb[condNames,]
	}
	
	e$clusterStack <- e$get.clusterStack( force=T )
	e$stats <- rbind( e$stats, e$get.stats() )
	return(e)
}

#' Get the pValues for all conditions in a given cluster.  OLD VERSION
#' 
#' @param e  The cMonkey environment
#' @param i  The cluster index to get the score for
#' @export
#' @usage pVals <- getClusterPVals( e, i )
getClusterPVals <- function( e, i ) {
	geneList <- e$clusterStack[[i]]$rows
	mean.sd <- e$means.sds[[as.character(length(geneList))]]
	pValList <- getRatPvals(e, geneList, mean.sd)
	return(pValList)
}

#' Get the pValues for all conditions in a given cluster
#' 
#' @param e  The cMonkey environment
#' @param i  The cluster index to get the score for
#' @export
#' @usage pVals <- getClusterPVals( e, i )
getClusterPVals <- function( e, i ) {
	geneList <- e$clusterStack[[i]]$rows
	if (length(geneList) == 0) {
		if (! is.null(attributes(e$ratios))) {
			pValList <- rep(1,attributes(e$ratios)$ncol)
			names(pValList) <- attributes(e$ratios)$cnames
		} else {
			expNames <- unique(unlist(sapply(e$clusterStack, function(x) {x$cols})))
			pValList <- rep(1,length(expNames))
			names(pValList) <- expNames
		}		
	} else {
		mean.sd <- e$means.sds[[as.character(length(geneList))]]
		pValList <- getRatPvals(e, geneList, mean.sd)
	}
	return(pValList)
}

#' Get the Mean pValue for all included conditions
#' 
#' @param e  A cMonkey environment
#' @param VERBOSE  Set to T to display if NA values are removed (Default: T)
#' @export
#' @usage env <- getMPV( e, VERBOSE=T )
getMPV <- function( e, VERBOSE=T ) {
	e <- update.means.sds(e)  #Make sure all back distributions exist

	relPvals <- NULL
	for (i in 1:length(e$clusterStack)){
		curConds <- e$clusterStack[[i]]$cols
		pVals <- getClusterPVals( e, i )
		relPvals <- c(relPvals, pVals[curConds])
	}
	if (any(is.na(relPvals))) {
		numOnes <- sum(sapply(e$clusterStack, function(x) {x$ncols})==1)
		if (VERBOSE) { cat(sum(is.na(relPvals)),"NAs removed\t", numOnes,"Size one clusters\n") }
	}
	MPV <- mean(relPvals,na.rm=T)
	return(MPV)
}

#' Return a data frame with "ClusterSize" "Residual" "pValue"
#' 
#' @param e  A cMonkey environment
#' @param VERBOSE  Set to T to display if NA values are removed (Default: T)
#' @export
#' @usage sizeDF <- getScoresVsize( e, VERBOSE=T )
getScoresVsize <- function( e, VERBOSE=T ) {
	e <- update.means.sds(e)  #Make sure all back distributions exist

	sizeDF <- NULL
	for (i in 1:length(e$clusterStack)){
		curConds <- e$clusterStack[[i]]$cols
		pVals <- getClusterPVals( e, i )
		pValue <- mean(pVals[curConds],na.rm=T)
		res <- e$clusterStack[[i]]$resid
		n <- length(e$clusterStack[[i]]$rows)
		sizeDF <- rbind(sizeDF, data.frame(n=n, res=res, pValue=pValue))
	}
	return(sizeDF)
}



#' Determine the background distributions for experiments not in the original data set.
#' 
#' @param clusterStack  A cMonkey clusterStack
#' @param newRatios  The ratios matrix for the new experimental conditions
#' @export
#' @usage new.means.sds <- getNew.means.sds( e$clusterStack, newRatios )
getNew.means.sds <- function( clusterStack, newRatios) {
	
	#Note: need to iterate through all possible cluster sizes.
	clusterSizes <- sort(unique(unlist(lapply(clusterStack, function(x) {x$nrows}))))
	new.means.sds <- foreach (n = clusterSizes) %dopar% {
		cat(n,',',sep="")
		getVarianceMeanSD.DF(list(delRatios=newRatios), n)
	}
	cat('\n')
	names(new.means.sds) <- clusterSizes
	return(new.means.sds)
}

#' Return the probablities that new experiments are from the probability distribution 
#'   of the training conditions
#' 
#' @param e  A cMonkey environment
#' @param newRatios  The ratios matrix for the new experimental conditions
#' @param condNames  The ratios names to use for the comparison (DEFAULT: c("oleate.lte90m","oleate.gt90m"))
#' @param useMTC  Set to true to multiply pVals my multiple testing correction (DEFAULT: F)
#' @param newClusterStack  Include a new clusterStack to modify gene membership (DEFAULT: NULL)
#' @export
#' @usage diffPvals <- getNewPvals( e, newRatios, condNames=c("oleate.lte90m","oleate.gt90m"), useMTC=F )
getNewPvals <- function( e, newRatios, condNames=c("oleate.lte90m","oleate.gt90m"), useMTC=F, newClusterStack=NULL) {
	
	if ( is.null(newClusterStack) ) {
		newClusterStack <- e$clusterStack
	} else {
		#Removing genes from a cluster should not cause means.sds recalc.
		if ( max(unlist(lapply(newClusterStack, function(x) {length(x$rows)}))) > length(e$means.sds) ) {
			cat("clusterStack contains clusters with more genes than calculated in means.sds\n")
			cat("Tell the lazy programmer to recalc means.sds in this case\n")
			return(NULL)
		}
	}
	if( any(!(colnames(newRatios) %in% rownames(e$means.sds[[1]]))) ) {
		new.means.sds <- getNew.means.sds( newClusterStack, newRatios )
	}

	#Calculate the background distributions
	#new.means.sds <- getNew.means.sds( newClusterStack, newRatios )
	
	#Calculate the pVals for each condition on each cluster
	origE <- list(ratios=e$ratios[condNames], clusterStack=newClusterStack, means.sds=e$means.sds)
	newE <- list(ratios=list(newRatios), clusterStack=newClusterStack, means.sds=new.means.sds)
	ks <- sapply(e$clusterStack,function(x) {x$k})
	diffPvals <- foreach (k = ks) %dopar% {
		pVals.orig <- getClusterPVals( origE, k )
		pVals.new <- getClusterPVals( newE, k )
		
		#Calculate the prob of being part of the same dist
		dPvals <- sapply(pVals.new,function(x) {mean(x < pVals.orig)})
		if (useMTC==T) { 
			dPvals <- sapply(dPvals,function(x) {x*length(dPvals)})
			dPvals[dPvals>1] <-1
		}
		return(dPvals)
	}
	names(diffPvals) <- as.character(ks)
	return(diffPvals)
}
