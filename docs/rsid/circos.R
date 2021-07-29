INF <- Sys.getenv("INF")
f <- file.path(INF,"work","INF1.merge")
INF1_merge <- read.delim(f)[c("Chrom","Start","End","prot","MarkerName")]

library(dplyr)
cvt <- read.csv(file.path(INF,"work","/INF1.merge.cis.vs.trans"),as.is=TRUE) %>%
                rename(MarkerName=SNP) %>% 
                mutate(chr=Chr,chrom=paste0("hs",Chr),start=bp,end=bp,p.chrom=paste0("hs",p.chr),value=-log10p,
                       fcolor=ifelse(cis,"color=vdred","color=vdblue"),
                       lcolor=ifelse(cis,"color=lred","color=lblue"),
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
                  left_join(gap::inf1[c("prot","target.short")]) %>%
                  arrange(Chr,bp)

pQTLs <- select(INF1_merge_cvt,chrom,start,end,value,fcolor)
write.table(pQTLs,file="pQTLs.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
pQTL_labels <- filter(INF1_merge_cvt,gene!=".") %>%
               mutate(gene=paste0(gene,"/",prot)) %>%
               select(chrom,start,end,gene,fcolor)
write.table(pQTL_labels,file="pQTL_labels.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
pQTL_links <- filter(INF1_merge_cvt,!cis) %>%
              select(p.chrom,p.start,p.end,chrom,start,end,lcolor)
write.table(pQTL_links,file="pQTL_links.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)

## not used:

pQTL_source_genes <- filter(INF1_merge_cvt,cis) %>%
                     select(p.chr,p.chrom,,p.start,p.end,p.gene) %>%
                     rename(chr=p.chr,chrom=p.chrom,start=p.start,end=p.end,gene=p.gene) %>%
                     distinct()
pQTL_target_genes <- filter(INF1_merge_cvt,!cis) %>%
                     select(chr,chrom,start,end,gene)
pQTL_genes <- rbind(pQTL_source_genes,pQTL_target_genes) %>%
              filter(gene!=".") %>%
              arrange(chr,start) %>%
              select(-chr) %>%
              distinct()
write.table(pQTL_genes,file="pQTL_genes.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
