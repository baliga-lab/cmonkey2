#!/usr/bin/env Rscript
suppressMessages(library('RSQLite'))
suppressMessages(library('getopt'))

# convert_cmonkey1.R
# This tools converts a cmonkey 1 run to a cmonkey2 sqlite result file.

create.tables <- function (con) {
  message("creating tables...")
  rs <- dbSendQuery(con, 'create table run_infos (start_time timestamp, finish_time timestamp, num_iterations int, last_iteration int, organism text, species text, ncbi_code int, num_rows int, num_columns int, num_clusters int, git_sha text)')
  dbClearResult(rs)
  rs <- dbSendQuery(con, 'create table row_members (iteration int,cluster int, order_num int)')
  dbClearResult(rs)
  rs <- dbSendQuery(con, 'create table column_members (iteration int,cluster int,order_num int)')
  dbClearResult(rs)
  rs <- dbSendQuery(con, 'create table cluster_stats (iteration int,cluster int,num_rows int,num_cols int,residual decimal)')
  dbClearResult(rs)
  rs <- dbSendQuery(con, 'create table row_names (order_num int, name text);')
  dbClearResult(rs)
  rs <- dbSendQuery(con, 'create table column_names (order_num int, name text);')
  dbClearResult(rs)
}

convert.rowscols <- function (con, e) {
  message("writing row names...")
  # copy rows
  db.rows <- data.frame(attr(e$ratios, 'rnames'))
  names(db.rows) <- c('name')
  db.rows$order_num <- seq(1, dim(db.rows)[1])
  apply(db.rows, 1, function (row) {
    rs <- dbSendPreparedQuery(con,
                              'insert into row_names (order_num,name) values (:order_num,:name)',
                              data.frame(order_num=as.numeric(row['order_num']), name=row['name']))
    dbClearResult(rs)
  })

  # copy columns
  message("writing column names...")
  db.cols <- data.frame(attr(e$ratios, 'cnames'))
  names(db.cols) <- c('name')
  db.cols$order_num <- seq(1, dim(db.cols)[1])
  for (col in db.cols) {
    rs <- dbSendPreparedQuery(con,
                              'insert into column_names (order_num,name) values (:order_num,:name)',
                              data.frame(order_num=as.numeric(col['order_num']), name=col['name']))
    dbClearResult(rs)
  }
  list(rows=db.rows, cols=db.cols)
}

convert.run_info <- function (con, e) {
  message("writing run info...")

  rs <- dbSendPreparedQuery(con,
                            'insert into run_infos (start_time,finish_time,num_iterations,last_iteration,organism,species,ncbi_code,num_rows,num_columns,num_clusters) values (:start_time,:finish_time,:num_iterations,:last_iteration,:organism,:species,:ncbi_code,:num_rows,:num_cols,:num_clusters)',                            
                            data.frame(start_time=e$DATE, finish_time=e$date.run, num_iterations=e$n.iter, last_iteration=e$iter,
                                       organism=e$organism, species=e$rsat.species, ncbi_code=e$taxon.id, num_rows=dim(e$ratios[[1]])[1],
                                       num_cols=dim(e$ratios[[1]])[2],
                                       num_clusters=e$k.clust))
  dbClearResult(rs)
}

convert.clusters <- function (con, e, db.rows, db.cols) {
  message("converting clusters...")
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
    
    # cluster stats
    rs <- dbSendPreparedQuery(con,
                              'insert into cluster_stats (iteration,cluster,num_rows,num_cols,residual) values (:iteration,:cluster,:num_rows,:num_cols,:residual)',
                              data.frame(iteration=e$iter, cluster=cluster$k, num_rows=cluster$nrows, num_cols=cluster$ncols, residual=cluster$resid))
    dbClearResult(rs)
  }
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
  'ratios', 'r',  2, "character",
  'rdata', 'i',  2, "character"
  ), byrow=TRUE, ncol=4)

opt <- getopt(spec)
if (is.null(opt$rdata) || is.null(opt$outdb) || is.null(opt$ratios)) {
  message(getopt(spec, usage=TRUE))
  q(status=1)
} else {
  message("loading input data...")
  suppressMessages(load(opt$rdata))
  make.rundb(opt$outdb, e)
  write.table(e$ratios, file=opt$ratios, sep='\t')
  message("Done.")
}
