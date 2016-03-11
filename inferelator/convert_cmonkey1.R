#!/usr/bin/env Rscript
suppressMessages(library('RSQLite'))
suppressMessages(library('getopt'))

# convert_cmonkey1.R
# This tools converts a cmonkey 1 run to a cmonkey2 sqlite result file.


convert.rowscols <- function (con, e) {
  # copy rows
  db.rows <- data.frame(attr(e$ratios, 'rnames'))
  names(db.rows) <- c('name')
  db.rows$order_num <- seq(1, dim(db.rows)[1])
  dbWriteTable(con, 'row_names', db.rows, row.names=FALSE)

  # copy columns
  db.cols <- data.frame(attr(e$ratios, 'cnames'))
  names(db.cols) <- c('name')
  db.cols$order_num <- seq(1, dim(db.cols)[1])
  dbWriteTable(con, 'column_names', db.cols, row.names=FALSE)
  list(rows=db.rows, cols=db.cols)
}

convert.run_info <- function (con, e) {
  run_info <- data.frame(t(c(e$DATE, e$date.run, e$n.iter, e$iter, e$organism,
                             e$rsat.species, e$taxon.id,
                             dim(e$ratios[[1]])[1], dim(e$ratios[[1]])[2],
                             e$k.clust, NA)))
  names(run_info) <- c('start_time', 'finish_time', 'num_iterations',
                       'last_iteration', 'organism', 'species', 'ncbi_code',
                       'num_rows', 'num_columns', 'num_clusters', 'git_sha')
  dbWriteTable(con, 'run_infos', run_info, row.names=FALSE)
}

convert.clusters <- function (con, e, db.rows, db.cols) {
  cluster.stack <- e$get.clusterStack()
  for (cluster in cluster.stack) {
    #tmp <- c(cluster$k, cluster$resid, cluster$rows, cluster$cols, cluster$resid)
    # cluster$p.clust
    # cluster$e.val matrix (2 motifs * sequence types)

    # row members first
    order.nums <- db.rows[db.rows$name %in% cluster$rows,]$order_num
    for (on in order.nums) {
      rs <- dbSendPreparedQuery(con,
                                'insert into row_members (iteration,cluster,order_num) values (:iteration,:cluster,:order_num)',
                                data.frame(iteration=e$iter, cluster=cluster$k, order_num=on))
      dbClearResult(rs)
    }

    # column members
    order.nums <- db.cols[db.cols$name %in% cluster$cols,]$order_num
    for (on in order.nums) {
      rs <- dbSendPreparedQuery(con,
                                'insert into column_members (iteration,cluster,order_num) values (:iteration,:cluster,:order_num)',
                                data.frame(iteration=e$iter, cluster=cluster$k, order_num=on))
      dbClearResult(rs)
    }
  }
}

create.tables <- function (con) {
  rs <- dbSendQuery(con, 'create table row_members (iteration int,cluster int, order_num int)')
  dbClearResult(rs)
  rs <- dbSendQuery(con, 'create table column_members (iteration int,cluster int, order_num int)')
  dbClearResult(rs)
}

make.rundb <- function (db.filename, e) {
  sqlite <- dbDriver("SQLite")
  con <- dbConnect(sqlite, dbname=db.filename)
  create.tables(con)
  convert.run_info(con, e)
  rowcols <- convert.rowscols(con, e)
  convert.clusters(con, e, rowcols$rows, rowcols$cols)

  dbDisconnect(con)
}

spec = matrix(c(
  'outdb', 'o',  2, "character",
  'rdata', 'i',  2, "character"
  ), byrow=TRUE, ncol=4)

opt <- getopt(spec)
if (is.null(opt$rdata) || is.null(opt$outdb)) {
  message(getopt(spec, usage=TRUE))
  q(status=1)
} else {
  suppressMessages(load(opt$rdata))
  make.rundb(opt$outdb, e)
}
