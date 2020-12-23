dummy <- function()
{
  require(plyr)
  require(doMC)
  doMC::registerDoMC(cores = 14)
  require(gdata)
  require(data.table)
  result <- rbindlist(alply(sentinels[1,], 1, function(obs) {
          tryCatch({
                  x <- as.numeric(obs)
                  return(x)
          }, error = function(e) {
                  dummyresult <- dummy
                  return(dummyresult)
          })}, .progress = "none", .parallel = TRUE))

  result <- as.data.frame(result)
  withRestarts(invokeRestart("foo", 1, 2), foo = function(x, y) {x + y})
}

