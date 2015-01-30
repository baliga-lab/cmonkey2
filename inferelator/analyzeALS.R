source('cmonkey-python.R')
source('bicluster.by.variance.R')
if (require(doMC, quietly = T) && !getDoParRegistered() && require(multicore) ) { registerDoMC(cores = multicore:::detectCores()) }

clusterStack <- read.cmonkey.sqlite('../out/cmonkey_run.db')

#Before resplitting the clusters, find out if any of them naturally has a good split between HEMI and NCAR
ratios <- read.table(gzfile('../out/ratios.tsv.gz'), row.names=1)
nHEMI <- sum(grepl('HEMI', colnames(ratios)))
nNCAR <- sum(grepl('NCAR', colnames(ratios)))

ncarCount <- hemiCount <- rep(0,length(clusterStack))
names(ncarCount) <- names(hemiCount) <- sapply(clusterStack, function(x) {x$k})
for (i in 1:length(clusterStack)) {
	hemiCount[as.character(i)] <- sum(grepl('HEMI', clusterStack[[i]]$cols))
	ncarCount[as.character(i)] <- sum(grepl('NCAR', clusterStack[[i]]$cols))
}
difs <- hemiCount - ncarCount

nrows <- sapply(clusterStack, function(x) {x$nrows})
names(nrows) <- sapply(clusterStack, function(x) {x$k})
realClusters <- names(nrows)[nrows>0]
difs <- difs[realClusters]
difs <- difs[order(abs(difs))]


#646 100 149 354 360 381 413 634   6  86 174 322 475 882   4  61 407 499 837 572
# 12 -13  13  13 -13 -13 -13  13 -14 -14 -14 -14 -14 -14 -15 -15 -16 -17 -17 -18
#407 has KLF4 binfing site, which is a big stem-cell related protein.


#Resplit clusters
e <- list(clusterStack=clusterStack, ratios=list(ratios=ratios))
new.means.sds <- getNew.means.sds( clusterStack, ratios )
clusterStack.resplit <- resplitClusters(e, means.sds=new.means.sds) #Most of the clusters include all conditions.  Reasons unknown

clusterStack.bk <- clusterStack
clusterStack <- clusterStack.resplit$newClusterStack

#Find the clusters that are most different between normal and disease
ncarCount <- hemiCount <- rep(0,length(clusterStack))
names(ncarCount) <- names(hemiCount) <- sapply(clusterStack, function(x) {x$k})
for (i in 1:length(clusterStack)) {
	hemiCount[as.character(i)] <- sum(grepl('HEMI', clusterStack[[i]]$cols))
	ncarCount[as.character(i)] <- sum(grepl('NCAR', clusterStack[[i]]$cols))
}
difs <- hemiCount - ncarCount

nrows <- sapply(clusterStack, function(x) {x$nrows})
names(nrows) <- sapply(clusterStack, function(x) {x$k})
realClusters <- names(nrows)[nrows>0]
difs <- difs[realClusters]
difs <- difs[order(abs(difs))]

#
#

#Calculate the expected number of clusters different between normal and disease.

#Pick some example to say something biological