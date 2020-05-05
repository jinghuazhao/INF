# 5-5-2019 JHZ

options(scipen=20, width=2000)
f <- Sys.getenv("pr")
s <- Sys.getenv('study')
n <- as.numeric(Sys.getenv("N"))
data_type <- Sys.getenv("data_type")
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
  X.ref <- t(apply(d$data,1:2,"%*%",2:0))
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
save(jam,file=paste0(f,"-jam.rda"))
pst <- slot(jam, "posterior.summary.table")
tm <- TopModels(jam)
sr <- data.frame(snp=snp, rsid=with(sumstats, rsid))
cs <- CredibleSet(jam, credible.percentile.threshold=0.75)
msbf <- ModelSizeBayesFactors(jam)[[1]]
n.col <- ncol(tm)
n.snps <- n.col-1
post.prob <- tm[,n.col]
n.sel <- apply(tm[,1:n.snps],1,sum)
tm1 <- tm[1,-n.col]
selected <- names(tm1[tm1==1])
if(n.sel[1]>0&n.sel[1]!=n.snps)
{
   PostProb_model <- rep(post.prob[1],n.sel[1])
   t <- cbind(subset(sr,snp%in%selected), PostProb_model, subset(pst,rownames(pst)%in%selected))
   write.table(t,paste0(f,".sel"),row.names=FALSE,quote=FALSE)
}
# xlsx
require(openxlsx)
xlsx <- paste0(f,"-jam.xlsx")
wb <- createWorkbook(xlsx)
addWorksheet(wb, "ID")
writeDataTable(wb, "ID", sr)
addWorksheet(wb, "TopModels")
writeDataTable(wb, "TopModels", as.data.frame(tm))
addWorksheet(wb, "Model.1")
PostProb_model <- rep(post.prob[1],n.sel[1])
writeDataTable(wb, "Model.1", cbind(subset(sr,snp%in%selected),PostProb_model,subset(pst,rownames(pst)%in%selected)))
addWorksheet(wb, "CredibleSet")
writeDataTable(wb, "CredibleSet", cbind(subset(sr,snp%in%cs),subset(pst,rownames(pst)%in%cs)))
addWorksheet(wb, "ModelSizeBayesFactors")
writeDataTable(wb, "ModelSizeBayesFactors", as.data.frame(msbf))
addWorksheet(wb, "posterior.summary.table")
writeDataTable(wb, "posterior.summary.table", cbind(ID=rownames(pst), as.data.frame(pst)))
saveWorkbook(wb, file=xlsx, overwrite=TRUE)
