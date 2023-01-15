# Meyer data
function meyer_test()
{
set.seed(1234567)
meyer <- within(meyer,
        {
          yNa <- y
           g1 <- ifelse(generation == 1, 1, 0)
           g2 <- ifelse(generation == 2, 1, 0)
           id <- animal
       animal <- ifelse(!is.na(animal), animal, 0)
          dam <- ifelse(!is.na(dam), dam, 0)
         sire <- ifelse(!is.na(sire), sire, 0)
        })
G <- kin.morgan(meyer)$kin.matrix * 2
library("regress")
r <- regress(y ~ -1 + g1 + g2, ~ G, data = meyer)
r
attach(meyer)
X <- as.matrix(meyer[c("g1", "g2")])
m <- BLR(yNa, XF = X, GF = list(ID = 1:nrow(G), A = G),
         prior = list(varE = list(df = 1, S = 0.25),
         varU = list(df = 1, S = 0.63)),
         nIter = 5000, burnIn = 500, thin = 1, saveAt = "meyer.BLR")
with(r, h2G(sigma, sigma.cov))
}

library("BLR")
library(gap)

# INTERVAL data
ReadGRMBin <- function(prefix, AllN=FALSE, size=4)
{
  BinFileName <- paste(prefix,".grm.bin",sep="")
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  id <- read.table(IDFileName)
  n <- dim(id)[1]
  BinFile <- file(BinFileName, "rb");
  grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  close(BinFile)
  NFile <- file(NFileName, "rb");
  if(AllN) N <- readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  else N <- readBin(NFile, n=1, what=numeric(0), size=size)
  close(NFile)
  i <- sapply(1:n, function(i) i*(i+1)/2)
  GRM <- matrix(NA,n,n)
  GRM[upper.tri(GRM,diag=TRUE)] <- grm
  GRM[lower.tri(GRM)] <- t(GRM)[lower.tri(GRM)]
  invisible(list(grm=grm, id=id, N=N, GRM=GRM))
}
INF <- Sys.getenv("INF")
system("sed '2d' /home/jhz22/INF/INTERVAL/o5000-inf1-outlier_in-r2.sample > s.sample")
f <- file.path(INF,"INTERVAL","o5000-inf1-outlier_in-r2.sample")
s <- read.table("s.sample",header=TRUE)
covar <- c("sexPulse","season","plate")
qcovar <- c("age","bleed_to_process_time",paste0("PC",1:20))
big3 <- c("CCL25","IL.18R1","TNFB")
header3 <- c("CCL25___O15444","IL.18R1___Q13478","TNFB___P01374")
prefix <- file.path(INF,"h2","INTERVAL")
G <- ReadGRMBin(prefix)
GRM <- with(G,GRM)
m <- as.formula(paste("CCL25___O15444","~",paste(covar,collapse="+"),"+",paste(qcovar,collapse="+"),"+","GRM"))
i <- regress(m, data = s)
i
attach(s)
X <- as.matrix(s[c(covar, qcovar)])
keep <- !is.na(apply(X,1,sum))
y <- CCL25___O15444[keep]
XF <- X[keep,]
A <- with(G,GRM)[keep,keep]+0.02*diag(length(keep[keep]))
m <- BLR(y, XF = XF, GF = list(ID = 1:nrow(A), A = A),
         prior = list(varE = list(df = 1, S = 0.25),
         varU = list(df = 1, S = 0.63)),
         nIter = 5000, burnIn = 500, thin = 1, saveAt = "CCL25___O15444.BLR")
with(m, h2G(sigma, sigma.cov))
detach(s)
