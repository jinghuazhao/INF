gsmr_efo <- function(root="gsmr-efo")
{
  INF <- Sys.getenv("INF")
  suppressMessages(library(dplyr))
  library(stringr)
  gsmr <- read.delim(file.path(INF,"mr","gsmr",paste0(root,".txt"))) %>%
          mutate(outcome=paste0(Disease),
                 exposure=protein,
                 z=bxy/se,
                 or=exp(bxy)
          ) %>%
          left_join(gap.datasets::inf1[c("target.short","gene")],by=c("exposure"="target.short")) %>%
          select(gene,outcome,z,or,bxy,se,p,nsnp,fdr)
  gene <- sort(unique(with(gsmr,gene)))
  outcome <- sort(unique(with(gsmr,outcome)))
  n <- length(gene)
  m <- length(outcome)
  gsmr_mat <- matrix(NA,m,n)
  colnames(gsmr_mat) <- gene
  rownames(gsmr_mat) <- outcome
  gsmr_mat_id <- gsmr_mat_fdr <- gsmr_mat
  for(k in 1:nrow(gsmr))
  {
     t <- gsmr[k,c("gene","outcome","or","fdr")]
     i <- t[["outcome"]]
     j <- t[["gene"]]
     gsmr_mat[i,j] <- t[["or"]]
     gsmr_mat_fdr[i,j] <- t[["fdr"]]
     gsmr_mat_id[i,j] <- k
  }
  rownames(gsmr_mat) <- gsub("\\b(^[a-z])","\\U\\1",rownames(gsmr_mat),perl=TRUE)
  rm(gene,outcome)
  options(width=200)
  gsmr_fdr01 <- subset(gsmr,fdr<=0.01)
  deprioritise <- c(FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,
                    FALSE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE)
  alist <- (1:nrow(gsmr_fdr01))[!deprioritise]
  blist <- (1:nrow(gsmr_fdr01))[deprioritise]
  library(grid)
  library(pheatmap)
  png(file.path(INF,"mr","gsmr",paste0(root,".png")),res=300,width=30,height=18,units="in")
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))),
          action="prepend")
  dn <- matrix(ifelse(!is.na(gsmr_mat) & abs(gsmr_mat_fdr) <= 0.01,
               ifelse(gsmr_mat_id%in%alist,unicode2[1],unicode2[2]), ""), nrow(gsmr_mat))
  pheatmap(gsmr_mat,cluster_rows=FALSE,cluster_cols=FALSE,angle_col="315",fontsize_row=35,fontsize_col=28,
           display_numbers = dn, fontsize_number=15)
  setHook("grid.newpage", NULL, "replace")
  grid.text("Proteins", y=-0.07, gp=gpar(fontsize=48))
  grid.text("Immune-mediated outcomes", x=-0.07, rot=90, gp=gpar(fontsize=48))
  dev.off()
  write.table(colnames(gsmr_mat),quote=FALSE,row.names=FALSE)
}

# https://www.rapidtables.com/code/text/unicode-characters.html (Page 38)
unicode1 <- c("\u25FC","\u25FB")
unicode2 <- c("\u26AB","\u25EF")
gsmr_efo()
gsmr_efo("gsmr-efo-reduce")
