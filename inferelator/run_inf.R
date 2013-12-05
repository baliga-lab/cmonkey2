#!/usr/bin/Rscript
library('getopt')

#####################################################################
#### run_inf.R - a command line script for Inferelator
####
#### This is a command line tool that takes takes required data from
#### the command line arguments, runs Inferelator and returns
#### a set of bicluster -> influence mappings in JSON format.
####
#### It can be used in 2 ways:
####
#### 1. provide a list of transcription factors and a cmonkey-python
####    output directory.
#### Example: ./run_inf.R --tfsfile <tfsfile> --resultdir <resultdir> --outfile <outfile>
####
#### 2. provide an R environment as an RData file that contains the
####    following parameters:
####    - ratios: the gene expression matrix
####    - predictors: the list of transcription factors
####    - e: an environment with the following data
####      - clusterStack: cMonkey standard cluster stack format
#### Example: ./run_inf.R --tfsfile <tfsfile> --rdata <rdatafile> --outfile <outfile>
####
#### Note: This script assumes that the necessary R packages and their
#### ----- dependencies are installed, namely
####       - cMonkeyNwInf
####       - RJSONIO
####       - RSQLite
####       - getopt
####
######################################################################

run.inferelator <- function(tfsfile, outfile, resultdir=NULL, rdata=NULL) {
  source('cmonkey-python.R')
  library('cMonkeyNwInf')
  library('RJSONIO')

  ge <- globalenv()
  ge$predictors <- readLines(tfsfile)

  if (!is.null(resultdir)) {
    message('Creating cluster stack from result directory...')
    ratio.filename <- paste(resultdir, 'ratios.tsv.gz', sep='/')
    result.filename <- paste(resultdir, 'cmonkey_run.db', sep='/')
    ge$ratios <- read.table(gzfile(ratio.filename), header=T, as.is=T, row.names=1)
    e <- environment()
    e$clusterStack <- read.cmonkey.sqlite(result.filename)
    e$k.clust <- length(e$clusterStack)
    e$envMap <- NULL
    e$colMap <- NULL
    ge$e <- e
    message('done.')
  } else {
    message('Loading cluster stack and ratios from environment...')
    load(rdata, envir=ge)
    # if these required variables are not defined, we simply set them
    # to meaningful defaults
    if (!exists('k.clust', envir=e)) e$k.clust <- length(e$clusterStack)
    if (!exists('envMap', envir=e))  e$envMap <- NULL
    if (!exists('colMap', envir=e))  e$colMap <- NULL
    e <- ge$e
    message('done.')
  }
  coeffs <- runnit.wrapper.halo(e, cv.choose="min+4se", tf.groups=999, alpha=0.8,
                                n.boot=1, tau=10,
                                r.cutoff=2, r.filter=0.8, weighted=T, aic.filter=25, plot=F)
  influences <- sapply(coeffs,
                       function(c) { if (class(c) != 'character') c$coeffs else emptyNamedList })
  names(influences) <- 1:length(influences)
  fileConn <- file(outfile)
  writeLines(toJSON(influences), fileConn)
  close(fileConn)
}

spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help',    'h', 0, "logical",
  'tfsfile', 'tfs', 1, "character",
  'resultdir', 'res', 2, "character",
  'rdata', 'rdat', 2, "character",
  'outfile', 'o', 1, "character"
  ), byrow=TRUE, ncol=4)

opt <- getopt(spec, usage=FALSE)
if (is.null(opt$verbose)) opt$verbose = 0
if (is.null(opt$help)) opt$help = FALSE

if (is.null(opt$tfsfile) || is.null(opt$outfile) ||
    is.null(opt$resultdir) && is.null(opt$rdata)) {
  print(getopt(spec, usage=TRUE))
} else {
  run.inferelator(opt$tfsfile, opt$outfile, opt$resultdir, opt$rdata)
}
