library("lars")
library("glmnet")
library("Matrix")
require("multicore")
require("foreach")
require("doMC")
require("ff")
require("igraph")
require('doMC')

DEBUG <- T
source('SD.DR.Inferelator.pkg.R')
source('cmonkey-python.R')
source('extract-biclusters.R')

inferelate <- function (ratios, clusterStack, tf.file) {
  tau <- 0
  tanayData <- ratios
  predictors <- readLines(tf.file)
  colMap <- NULL
  n.clust <- 1:length(clusterStack)

  coeffs <- runnit(ks=n.clust, data=tanayData, col.map=colMap, predictors=predictors,
                   clusterStack=clusterStack, tau = as.numeric(tau),
                   filter.pred.by.col=F, plot=F, tf.groups=999, alpha=0.8,
                   cv.choose="min+3se", aic.filter=NA, r.cutoff=2, funcs=NA, n.boot=2)
  coeffs
}

# put 3 steps into a single call
doit <- function(result.dir, tf.file, network.id, species.id) {
  ratio.filename <- paste(result.dir, 'ratios.tsv.gz', sep='/')
  result.filename <- paste(result.dir, 'cmonkey_run.db', sep='/')

  ratios <- read.table(gzfile(ratio.filename), header=T, as.is=T, row.names=1)
  clusterStack <- read.cmonkey.sqlite(result.filename)
  coeffs <- inferelate(ratios, clusterStack, tf.file)
  extract.influences(con=NULL, coeffs, network.id, species.id)
}
