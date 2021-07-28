INF <- Sys.getenv("INF")
f <- file.path(INF,"work","INF1.merge")
INF1_merge <- read.delim(f)[c("Chrom","Start","End","prot","MarkerName")]

library(dplyr)
cvt <- read.csv(file.path(INF,"work","/INF1.merge.cis.vs.trans"),as.is=TRUE) %>%
                rename(MarkerName=SNP) %>% 
                mutate(chr=paste0("hs",Chr),start=bp-1,end=bp,value=-log10p,
                       color=ifelse(cis,"color=red","color=blue"),
                       chrbp=paste(Chr,bp,sep=":"))

geneinfo <- gene <- vector("character")
for(i in 1:180)
{
  geneinfo[i] <- with(ieugwasr::variants_chrpos(cvt$chrbp[i],5000),geneinfo)
  gene[i]=gsub(":([0-9])*","",geneinfo[i])
}
annotate <- within(cvt,{gene=gene})
is.cis <- with(annotate,cis)
annotate[is.cis,"gene"] <- annotate[is.cis,"p.gene"]
INF1_merge_cvt <- merge(INF1_merge,annotate,by=c("prot","MarkerName")) %>%
                  arrange(Chr,bp)

pQTLs <- select(INF1_merge_cvt,chr,start,end,value)
write.table(pQTLs,file="pQTLs.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
pQTL_labels <- select(INF1_merge_cvt,chr,start,end,gene,color)
write.table(pQTL_labels,file="pQTL_labels.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
