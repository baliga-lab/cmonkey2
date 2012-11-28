# Helper functionality to work with cMonkey/Python data in R
#
# Explanation of the clusterStack object
# --------------------------------------
# clusterStack is a list containing |num clusters| elements
# Each element of this list is a list of
# 
# 1. $nrows: number of rows in cluster
# 2. $ncols: number of columns in cluster
# 3. $rows: row names in the cluster
# 4. $cols: column names in the cluster
# 5. $k: cluster number
# 6. $p.clust p value of motif scoring
#      vector of number, s elements, each element named like sequence type
# 7. $e.val e value of motif scoring
#      2xs matrix dimensions column name sequence type
#      s number of sequence types, columns named
# 8. $resid cluster residual
#      vector of number, 1 element, named 'ratios'
library('rjson')

read.cmonkey.json <- function(result.filename) {
  cmresult.file <- gzfile(result.filename)
  json.str <- scan(cmresult.file, character(0), sep='\n')
  json <- fromJSON(json.str)
  close(cmresult.file)

  # TODO: replace with a map/apply call
  result <- list()
  num.clusters <- length(json$rows)

  for (cluster in 1:num.clusters) {
    cluster.data <- list()
    
    # the real work
    cluster.data$nrows <- length(json$rows[[cluster]])
    cluster.data$ncols <- length(json$columns[[cluster]])

    # can do an unlist() to convert list to vector
    cluster.data$rows <- json$rows[[cluster]]
    cluster.data$cols <- json$columns[[cluster]]
    cluster.data$k <- cluster
    cluster.data$p.clust <- NULL
    cluster.data$e.val <- NULL
    cluster.data$resid <- json$residuals[[cluster]]

    # append the cluster data
    result[[cluster]] <- cluster.data
  }
  result
}
