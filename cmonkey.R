# An integration module to run cmonkey.py from within R
run.cmonkey <- function(organism, ratios) {
    command <- paste(c("./cmonkey.py", "--organism", organism, "--ratios", ratios), collapse=" ")
    message(command)
    system(command, wait=FALSE, ignore.stdout=TRUE, ignore.stderr=TRUE)
}
