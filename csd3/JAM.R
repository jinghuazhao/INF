# 8-8-2019 JHZ

options(scipen=20, width=2000)
f <- Sys.getenv("pr")
s <- Sys.getenv('study')
data_type <- Sys.getenv("data_type")
# summary statistics
sumstats <- read.table(paste0(f,".z"), as.is=TRUE, header=TRUE)
# snpid-rsid mapping
load(paste0(f,".rda"))
sumstats <- merge(sumstats, snpid_rsid, by="rsid")
# reference data
if (data_type == "bgen")
{
  require(rbgen)
  samples <- matrix(scan(paste0(s,".id"), what=c("","")), ncol=2, byrow=TRUE)
  d <- bgen.load(paste0(f,".bgen"), rsids=scan(paste0(f,".incl"), what=""))
  dimnames(d$data)[[2]] <- samples[,1]
  vid <- as.character(with(with(d, variants), rsid))
  R <- t(apply(d$data,1:2,"%*%",2:0))
} else if (data_type == "binary_ped" ) {
  library(plink2R)
  p <- read_plink(f)
  vid <- with(with(p, bim), V2)
  R <- with(p, as.data.frame(2-bed))
}
for (i in 1:ncol(R)) R[is.na(R[,i]), i] <- mean(R[,i], na.rm = TRUE)
ref <- list(rsid=vid,R=R)
X.ref <- with(ref, R)
ss <- subset(sumstats, rsid %in% with(ref, rsid))
beta <- with(ss, beta)
rsid <- with(ss, name)
snpid <- with(ss, rsid)
require(R2BGLiMS)
ssnpid <- paste0("snp", 1:length(beta))
names(beta) <- colnames(X.ref) <- ssnpid
priors <- list("a"=1, "b"=length(beta), "Variables"=ssnpid) # ssnpid
n <- as.numeric(Sys.getenv("N"))
j <- JAM(marginal.betas=beta, n=n, X.ref=X.ref, n.mil=5, tau=n, full.mcmc.sampling = FALSE, model.space.priors=priors)
save(j,file=paste0(f,".j"))
pst <- slot(j, "posterior.summary.table")
tm <- TopModels(j)
ssr <- data.frame(ssnpid=ssnpid, snpid=snpid, rsid=rsid)
cs <- CredibleSet(j, credible.percentile.threshold=0.75)
msbf <- ModelSizeBayesFactors(j)[[1]]
n.col <- ncol(tm)
n.snps <- n.col-1
post.prob <- tm[,n.col]
n.sel <- apply(tm[,1:n.snps],1,sum)
tm1 <- tm[1,-n.col]
selected <- names(tm1[tm1==1])
if(n.sel[1]>0&n.sel[1]!=n.snps)
{
   PostProb_model <- rep(post.prob[1],n.sel[1])
   t <- cbind(subset(ssr,ssnpid%in%selected), PostProb_model, subset(pst,rownames(pst)%in%selected))
   write.table(t,paste0(f,".sel"),row.names=FALSE,quote=FALSE)
}
require(openxlsx)
xlsx <- paste0(f,"-JAM.xlsx")
wb <- createWorkbook(xlsx)
addWorksheet(wb, "ID")
writeDataTable(wb, "ID", ssr)
addWorksheet(wb, "TopModels")
writeDataTable(wb, "TopModels", as.data.frame(tm))
addWorksheet(wb, "Model.1")
PostProb_model <- rep(post.prob[1],n.sel[1])
writeDataTable(wb, "Model.1", cbind(subset(ssr,ssnpid%in%selected),PostProb_model,subset(pst,rownames(pst)%in%selected)))
addWorksheet(wb, "CredibleSet")
writeDataTable(wb, "CredibleSet", cbind(subset(ssr,ssnpid%in%cs),subset(pst,rownames(pst)%in%cs)))
addWorksheet(wb, "ModelSizeBayesFactors")
writeDataTable(wb, "ModelSizeBayesFactors", as.data.frame(msbf))
addWorksheet(wb, "posterior.summary.table")
writeDataTable(wb, "posterior.summary.table", cbind(ID=rownames(pst), as.data.frame(pst)))
saveWorkbook(wb, file=xlsx, overwrite=TRUE)
