# 3-10-2018 JHZ

KORA <- "/data/jinhua/KORA/Affy AxiomPhase3_n3775_CodeAX1KG3_V2_LU9220"

inf1 <- scan("doc/olink.prot.list.txt","character",skip=184)
UniProt <- unlist(lapply(strsplit(inf1,"___"),"[",2))
iid <- scan(paste0(KORA,"/","individuals.txt"),"character")
n <- length(iid)
p <- matrix(rnorm(n*92),n,92)
colnames(p) <- UniProt
PCA <- matrix(runif(n*4),n,4)
colnames(PCA) <- paste0("PC",1:4)
pheno <- cbind(ID_1=1:n,ID_2=iid,missing=runif(n)/n,sex=1+(runif(n)>0.5),PCA,p)
l2 <- c(rep("0",3),rep("C",5),rep("P",92))
write.table(rbind(l2,pheno),file="KORA.txt",quote=FALSE,row.names=FALSE)
