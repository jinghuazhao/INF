test <- function()
{
  library(gwascat)
  cur = makeCurrentGwascat()  # result varies by day
  data(cur)
  cur
  library(rtracklayer)
  path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
  ch = import.chain(path)
  ch
  seqlevelsStyle(cur) = "UCSC"  # necessary
  cur19 = liftOver(cur, ch)
  class(cur19)
}

plry_doMC <- function()
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

oneKG <- function()
{
  for chr in {1..22}
  do
    echo ${chr}
    export chr1KG=~/rds/public_databases/1000G/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
    qctool -g $chr1KG -og 1KG-${chr} -ofiletype binary_ped
 done
}
