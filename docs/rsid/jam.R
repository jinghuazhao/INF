# 18-1-2021 JHZ

options(scipen=20, width=2000)
f <- Sys.getenv("pr")
s <- Sys.getenv('study')
n <- as.numeric(Sys.getenv("N"))
k <- as.integer(Sys.getenv("k"))
# summary statistics
z <- read.table(paste0(f,".z"), as.is=TRUE, header=TRUE)
# reference data
require(rbgen)
samples <- matrix(scan(paste0(s,".id"), what=c("","")), ncol=2, byrow=TRUE)
b <- bgen.load(paste0(f,"-jam.bgen"), rsids=scan(paste0(f,"-jam.incl"), what=""))
dimnames(b$data)[[2]] <- samples[,1]
variants <- with(b,variants)
map <- data.frame(variants,ord=1:nrow(variants))
m <- merge(z,map,by="rsid")
X.ref <- t(apply(b$data,1:2,"%*%",0:2))
for (i in 1:ncol(X.ref)) X.ref[is.na(X.ref[,i]), i] <- median(X.ref[,i], na.rm = TRUE)
sumstats <- m[with(m,ord),]
# JAM
require(R2BGLiMS)
snp <- make.names(with(sumstats,rsid))
priors <- list("a"=1, "b"=nrow(sumstats), "Variables"=snp)
X <- with(sumstats,beta)
names(X) <- colnames(X.ref) <- snp
jam <- JAM(marginal.betas=X, n=n, X.ref=X.ref, n.mil=5, tau=n, model.space.priors=priors, trait.variance=1)
ref <- within(b, {variants <- subset(variants,rsid%in%snp);data <- subset(data,rownames(data)%in%snp)})
save(X,X.ref,ref,n,priors,jam,file=paste0(f,"-jam.rda"))
pst <- slot(jam, "posterior.summary.table")
tm <- TopModels(jam)
cs <- CredibleSet(jam, credible.percentile.threshold=0.95)
msbf <- ModelSizeBayesFactors(jam)[[1]]
png(paste0(f,"-jam.png"),height=12,width=8,units="in",res=300)
ManhattanPlot(jam)
dev.off()
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
addWorksheet(wb, "Manhattan plot", zoom=150)
insertImage(wb, "Manhattan plot", paste0(f,"-jam.png"),height=12,width=8)
saveWorkbook(wb, file=xlsx, overwrite=TRUE)
