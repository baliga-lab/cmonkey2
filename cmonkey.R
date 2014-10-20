# An integration module to run cmonkey.py from within R
library('RSQLite')
cmonkey.currentiter <- function(db.filename) {
  sqlite <- dbDriver("SQLite")
  dbConnect(sqlite, dbname=db.filename)
  run.info <- dbGetQuery(con, "select last_iteration, num_clusters from run_infos")[1,]
  dbDisconnect(con)
  run.info[,1]
}

run.cmonkey <- function(organism, ratios) {
    command <- paste(c("./cmonkey.py", "--organism", organism, "--ratios", ratios), collapse=" ")
    message(command)
    system(command, wait=FALSE, ignore.stdout=TRUE, ignore.stderr=TRUE)
}

