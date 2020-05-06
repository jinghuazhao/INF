# 6-5-2019 JHZ

options(scipen=20, width=2000)
data_type <- Sys.getenv("data_type")
f <- Sys.getenv("pr")
s <- Sys.getenv('study')
n <- as.numeric(Sys.getenv("N"))
# summary statistics
z <- read.table(paste0(f,".z"), as.is=TRUE, header=TRUE)
# reference data
if (data_type == "bgen")
{
  require(rbgen)
  samples <- matrix(scan(paste0(s,".id"), what=c("","")), ncol=2, byrow=TRUE)
  b <- bgen.load(paste0(f,"-jam.bgen"), rsids=scan(paste0(f,"-jam.incl"), what=""))
  dimnames(b$data)[[2]] <- samples[,1]
  vid <- as.character(with(with(b, variants), rsid))
  X.ref <- t(apply(b$data,1:2,"%*%",2:0))
} else if (data_type == "binary_ped" ) {
  require(plink2R)
  p <- read_plink(paste0(f,"-jam"))
  vid <- with(with(p, bim), V2)
  X.ref <- with(p, as.matrix(2-bed))
}
for (i in 1:ncol(X.ref)) X.ref[is.na(X.ref[,i]), i] <- median(X.ref[,i], na.rm = TRUE)
sumstats <- subset(z, rsid %in% vid)
# JAM
require(R2BGLiMS)
snp <- make.names(with(sumstats,rsid))
priors <- list("a"=1, "b"=nrow(sumstats), "Variables"=snp)
X <- with(sumstats,beta)
names(X) <- colnames(X.ref) <- snp
jam <- JAM(marginal.betas=X, n=n, X.ref=X.ref, n.mil=5, tau=n, full.mcmc.sampling = FALSE, model.space.priors=priors)
save(X,X.ref,n,priors,jam,file=paste0(f,"-jam.rda"))
pst <- slot(jam, "posterior.summary.table")
tm <- TopModels(jam)
cs <- CredibleSet(jam, credible.percentile.threshold=0.75)
msbf <- ModelSizeBayesFactors(jam)[[1]]
# xlsx
require(openxlsx)
xlsx <- paste0(f,"-jam.xlsx")
wb <- createWorkbook(xlsx)
addWorksheet(wb, "TopModels", zoom=150)
writeDataTable(wb, "TopModels", as.data.frame(tm), rowNames=TRUE)
addWorksheet(wb, "CredibleSet", zoom=150)
writeDataTable(wb, "CredibleSet", as.data.frame(subset(pst,rownames(pst)%in%cs)), rowNames=TRUE)
addWorksheet(wb, "ModelSizeBayesFactors", zoom=150)
writeDataTable(wb, "ModelSizeBayesFactors", as.data.frame(msbf), rowNames=TRUE)
addWorksheet(wb, "posterior.summary.table", zoom=150)
writeDataTable(wb, "posterior.summary.table", as.data.frame(pst), rowNames=TRUE)
saveWorkbook(wb, file=xlsx, overwrite=TRUE)
