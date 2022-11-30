INF <- Sys.getenv("INF")
suppressMessages(library(dplyr))
library(stringr)
mr <- read.delim(file.path(INF,"mr","gsmr","mr-efo-mr.tsv")) %>%
      mutate(disease=gsub("\\s\\(oligoarticular or rheumatoid factor-negative polyarticular\\)","",disease)) %>%
      mutate(disease=gsub("\\s\\(non-Lofgren's syndrome\\)","",disease)) %>%
      mutate(disease=gsub("_PSORIASIS","",disease)) %>%
      mutate(outcome=disease,
             exposure=gene,
             or=exp(b),
             group=cut(or,breaks=c(0,0.49,0.99,1.49,2))) %>%
      select(gene,outcome,or,b,se,pval,nsnp,fdr,group) %>%
      mutate(or=ifelse(!is.na(or) & or<=2,or,NA))
options(width=200)
subset(mr,fdr<=0.05)
gene <- unique(with(mr,gene))
outcome <- unique(with(mr,outcome))
n <- length(gene)
m <- length(outcome)
mr_mat <- matrix(NA,m,n)
colnames(mr_mat) <- gene
rownames(mr_mat) <- outcome
mr_mat_fdr <- mr_mat
for(k in 1:nrow(mr))
{
   t <- mr[k,c("gene","outcome","or","group","fdr")]
   i <- t[["outcome"]]
   j <- t[["gene"]]
   mr_mat[i,j] <- t[["or"]]
   mr_mat_fdr[i,j] <- t[["fdr"]]
}
rownames(mr_mat) <- gsub("\\b(^[a-z])","\\U\\1",rownames(mr_mat),perl=TRUE)
rm(gene,outcome)
library(grid)
library(pheatmap)
png(file.path(INF,"mr","gsmr","mr-efo.png"),res=300,width=30,height=18,units="in")
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))),
        action="prepend")
pheatmap(mr_mat,cluster_rows=FALSE,cluster_cols=FALSE,angle_col="315",fontsize_row=30,fontsize_col=30,
         display_numbers = matrix(ifelse(!is.na(mr_mat) & abs(mr_mat_fdr) <= 0.05, "*", ""), nrow(mr_mat)), fontsize_number=20)
setHook("grid.newpage", NULL, "replace")
grid.text("Proteins", y=-0.07, gp=gpar(fontsize=48))
grid.text("Immune-mediated outcomes", x=-0.07, rot=90, gp=gpar(fontsize=48))
dev.off()
write.table(colnames(mr_mat),quote=FALSE,row.names=FALSE)
