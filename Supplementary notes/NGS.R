#
suppressMessages(library(dplyr))
suppressMessages(library(openxlsx))

INF <- Sys.getenv("INF")
metal <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
         mutate(uniprot_rsid=paste0(uniprot,"-",rsid))
load(file.path(INF,"work","novelpqtls.rda"))
f <- 'https://www.biorxiv.org/content/biorxiv/early/2022/02/20/2022.02.18.481034/DC3/embed/media-3.xlsx'
S1 <- read.xlsx(f)
f <- 'https://www.biorxiv.org/content/biorxiv/early/2022/02/20/2022.02.18.481034/DC12/embed/media-12.xlsx'
S3 <- read.xlsx(f)
save(S1,S3,file=file.path(INF,"ukb","tables.rda"))
ss <- merge(S1,S3) %>%
      mutate(uniprot_rsid=paste0(uniprot,"-",rsName)) %>%
      filter(uniprot_rsid%in%metal$uniprot_rsid) %>%
      left_join(metal) %>%
      left_join(select(gap.datasets::inf1,prot,target,target.short)) %>%
      select(target.short,rsid,MarkerName,uniprot,cis.trans,Log10.pval.gc.cor.unadj) %>%
      transmute(Protein=target.short,Sentinels=rsid,UniProt=uniprot,SNPid=MarkerName,"cis/trans"=cis.trans,
                Proxies="as sentinel",r2=1,p=10^-Log10.pval.gc.cor.unadj,Source="Eldjarn, et al. (2022)",PMID="",Comment="")
options(width=200)
write.table(ss,file=file.path(INF,"ukb","NGS.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
