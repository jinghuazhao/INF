# 27-12-2022 JHZ

options(width=200)

prune <- function()
{
  p <- Sys.getenv("p")
  r <- Sys.getenv("r")
  pr <- paste0(p,"-",r)

  tbl <- read.delim(paste0(pr,".z"),sep=" ")
  suppressMessages(library(gap))
  z <- suppressMessages(cs(tbl))
  write(z[["MarkerName"]],file=paste0("~/INF/cs/prune/",pr,".cs"),nrow(z))
  write(format(z[["ppa"]],digits=3,scientific=TRUE),file=paste0("~/INF/cs/prune/",pr,".ppa"),nrow(z))
}

unprune <- function()
{
  p <- Sys.getenv("p")
  r <- Sys.getenv("r")
  pr <- paste0(p,"-",r)

  tbl <- read.delim(paste0(pr,".tbl.gz"))
  suppressMessages(library(gap))
  z <- suppressMessages(cs(tbl))
  write(z[["MarkerName"]],file=paste0("~/INF/cs/unprune/",pr,".cs"),nrow(z))
  write(format(z[["ppa"]],digits=3,scientific=TRUE),file=paste0("~/INF/cs/unprune/",pr,".ppa"),nrow(z))
}

# prune()
# tbl <- within(tbl,{logp <- logp(Effect/StdErr)})
# l <- cs(tbl,log_p="logp")

# library(corrcoverage)
# pp <- with(tbl,z0_pp(z=Effect/StdErr,f=MAF,type="quant",s=0.1,N=N))
# cs_pp <- credset(pp,thr=0.99)
# idx <- with(cs_pp,credset)
# write(tbl[["MarkerName"]][idx],file=paste0(pr,".chk"),length(idx))
