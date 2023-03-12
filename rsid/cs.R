# 15-1-2023 JHZ

options(width=200)

cs <- function(tbl, b="Effect", se="StdErr", log_p=NULL, cutoff=0.95)
# credible set based on METAL sumstats, copied here from R/gap
{
  requireNamespace("matrixStats")
  tbl <- within(tbl, {
           if (is.null(log_p)) z <- tbl[[b]]/tbl[[se]]
           else z <- qnorm(tbl[[log_p]]-log(2), lower.tail=FALSE, log.p=TRUE)
           z2 <- z * z / 2
           d <- matrixStats::logSumExp(z2)
           log_ppa <- z2 - d
           ppa <- exp(log_ppa)
        })
  ord <- with(tbl, order(ppa,decreasing = TRUE))
  tbl_ord <- within(tbl[ord,], {cppa <- cumsum(ppa)})
  last <- which(with(tbl_ord,cppa) >= cutoff)[1]
  tbl[ord[1:last],]
}

prune <- function()
{
  p <- Sys.getenv("p")
  r <- Sys.getenv("r")
  pr <- paste0(p,"-",r)

  tbl <- read.delim(paste0(pr,".z"),sep=" ")
  z <- suppressMessages(cs(tbl))
  write(z[["MarkerName"]],file=paste0("~/INF/cs/prune/",pr,".cs"),nrow(z))
  write(format(z[["ppa"]],digits=3,scientific=TRUE),file=paste0("~/INF/cs/prune/",pr,".ppa"),nrow(z))
}

unprune <- function()
{
  p <- Sys.getenv("p")
  r <- Sys.getenv("r")
  pr <- paste0(p,"-",r)

  suppressMessages(library(dplyr))
  tbl <- read.delim(paste0(pr,".tbl.gz")) %>%
         mutate(MAF=if_else(Freq1 <=0.05, Freq1, 1-Freq1),z=Effect/StdErr) %>%
         rename(f=MAF)
  z <- suppressMessages(cs(tbl))
  write(z[["MarkerName"]],file=paste0("~/INF/cs/unprune/",pr,".cs"),nrow(z))
  write(format(z[["ppa"]],digits=3,scientific=TRUE),file=paste0("~/INF/cs/unprune/",pr,".ppa"),nrow(z))
}

# prune()

test <- function()
{
  tbl <- within(tbl,{logp <- logp(Effect/StdErr)})
  l <- cs(tbl,log_p="logp")
  library(dplyr)
  library(corrcoverage)
  attach(tbl)
  pp <- z0_pp(z,f,"quant",N,0.2)
  detach(tbl)
  cs_pp <- credset(pp,thr=0.99)
  idx <- with(cs_pp,credset)
  write(tbl[["MarkerName"]][idx],file=paste0(pr,".chk"),length(idx))
}
