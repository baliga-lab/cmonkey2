# A script to convert cmonkey-python's JSON result into
# an RData environment which can be processed by inferelator
# Inferelator needs
#
# 1. input ratio matrix (2D-matrix, input to cmonkey)
# 2. clusterStack (output from cmonkey)
# 3. predictors (a list of genes)
#
source('cmonkey_python.R')

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  ratio.filename <- args[1]
  result.filename <- args[2]
  
  ratios <- read.table(gzfile(ratio.filename), header=T, as.is=T, row.names=1)
  clusterStack <- read.cmonkey.json(result.filename)
  print("Writing objects to R file...")
  save(ratios, clusterStack, file='output.RData')
  print("Done.")
} else {
  print("Usage: convert <ratio-file>")
  print(args)
}
