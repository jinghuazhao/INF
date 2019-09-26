# 26-9-2019 JHZ

outfile <- "INTERVAL/o5000-inf1-outlier_out-r2.sample"
colnames <- read.table(outfile, as.is=TRUE, nrows=1)
out <- read.table(outfile,skip=2,as.is=TRUE,col.names=colnames)
dim(out)
pdf("INTERVAL-pca.pdf")
ppc <- with(out, princomp(na.omit(out), cor=TRUE))
screeplot(ppc, npcs=20, type="lines", main="INTERVAL proteins PCA(invnorm(protein)) screeplot")
with(ppc, plot(scores[,1:2], main="INTERVAL PCA(invnorm(protein)) PC1 -- PC2"))
scores <- out[paste0("PC",1:20)]
names(scores) <- paste0("Comp.",1:20)
row.names(scores) <- with(out, ID_1)
dim(scores)
gpc <- list()
gpc <- within(gpc, {scores <- scores})
with(gpc, plot(scores[,1:2], main="INTERVAL genotypes PC1 -- PC2"))
sscores <- subset(scores, Comp.1 <= 0.005 & Comp.2 >= -0.02)
plot(sscores[,1:2], main="INTERVAL genotypes PC1 (<=0.005) -- PC2 (>=-0.02)")
dev.off()
dim(sscores)
id <- setdiff(with(out,ID_1),rownames(sscores))
write.table(cbind(id,id), file="INTERVAL/INTERVAL-excl.id", row.names=FALSE, col.names=FALSE)
