#!/usr/bin/bash

export proteomics_results=~/rds/results/public/proteomics/Fenland
export all=${proteomics_results}/all.grch37.tabix.gz
export v4=SomaLogicv4.tsv

function replication()
(
  gunzip -c ${all} | \
  head -1
  join -13 -22 <(cut -f1,4,7 --output-delimiter=' ' ${INF}/deCODE/${v4} | sort -k3,3) \
               <(cut -f2,4,5,20 --output-delimiter=' ' ${INF}/work/INF1.METAL | awk '{print $2":"$3,$4,$1}' | sort -k2,2) | \
  parallel -C' ' -j20 --env proteomics_results '
    tabix ${all} {4} | grep -w {2} | grep -w {5}
  '
) > ${INF}/Fenland/replication.tsv

function fenland()
(
echo -e "Protein\tSentinels\tUniProt\tSNPid\tcis/trans\tProxies\tr2\tp\tTarget Full Name\tSource\tPMID\tComment"
join -12 -21 <(awk 'NR>1 && $14<=5e-8' ${INF}/Fenland/replication.tsv | cut -f4,11,14 | sort -k2,2) \
             <(cut -f1,3,5,6 --complement --output-delimiter=' ' ${INF}/Fenland/${v4} | sed 's/-/_/' | sort -k1,1) | \
awk '{print $4"-"$2,$0}' | \
sort -k1,1 | \
join - <(awk 'NR>1 {print $20"-"$1,$3,$2,$21}' ${INF}/work/INF1.METAL | sort -k1,1) | \
awk 'a[$1]++==0' | \
cut -d' ' -f1,2 --complement | \
sort -k4,4 | \
tr ' ' '\t'| \
join -24 <(Rscript -e 'write.table(pQTLtools::inf1[c("prot","target.short")],col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")' | \
           sort -t$'\t' -k1,1) - | \
awk -vOFS='\t' '{print $2,$6,$5,$3,$7,"as sentinels",1,$4}' | \
sort -t$'\t' -k1,1 | \
join -t$'\t' - \
             <(Rscript -e 'write.table(data.frame(pQTLtools::inf1[c("target.short","target")],Source="Pietzner et al. (2021)",PMID="",Comment=""),
                                       col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")' | \
               sort -t$'\t' -k1,1)
)

fenland > ${INF}/Fenland/Fenland.tsv

# Selection by region

(
  gunzip -c ${all} | \
  head -1 | \
  awk -v OFS='\t' '{print $0,"Prot","MarkerName","sentinel"}'
  join -13 -21 <(cut -f1,4,7 --output-delimiter=' ' ${INF}/deCODE/${v4} | sort -k3,3) \
               <(Rscript -e '
                   INF <- Sys.getenv("INF")
                   suppressMessages(library(dplyr))
                   load(file.path(INF,"work","novel_data.rda"))
                   novel_data <- novel_data %>%
                                 mutate(region=paste0(Chromosome,":",Position-1e6,"-",Position+1e6)) %>%
                                        select(uniprot,region,prot,MarkerName,rsid)
                   write.table(novel_data,row.names=FALSE,col.names=FALSE,quote=FALSE)
                 ' | \
                 sort -k1,1) | \
  parallel -C' ' -j20 --env all 'tabix ${all} {4} | grep -w {2} | awk -v prot={5} -vsnpid={6} -vrsid={7} -vOFS="\t" "{print \$0,prot,snpid,rsid}"'
) | \
awk 'NR==1||$14<5e-8' | cut -f3,4,11-14,21-23> ${INF}/Fenland/region.tsv
Rscript -e '
  INF <- Sys.getenv("INF")
  suppressMessages(library(dplyr))
  region <- read.delim(file.path(INF,"Fenland","region.tsv")) %>%
            mutate(rsid=gsub("_[A-Z]*_[A-Z]*","",rsid))
  key <- Sys.getenv("LDLINK_TOKEN")
  sentinels <- unique(region$sentinel)
  blocks <- r <- list()
  for(i in 1:length(sentinels))
  {
    blocks[[i]] <- subset(region,sentinel==sentinels[i])
    sentinel_and_snps <- c(sentinels[i],blocks[[i]]$rsid[grepl("^rs",blocks[[i]]$rsid)])
    r[[i]] <- LDmatrix(snps=sentinel_and_snps,pop="EUR",r2d="r2",token=key)
    r2 <- subset(r[[i]],RS_number==sentinels[i])
    sel <- !is.na(r2) & r2>=0.8
    cat(sentinels[i],"\n")
    print(names(r2)[sel])
    print(r2[sel])
  }
'
# > sentinels
# [1] "rs28735437" "rs7213460"  "rs3184504"
# rs28735437
# [1] "RS_number"  "rs28735437"
# [1] "rs28735437" "1"
#
# rs7213460
# [1] "RS_number" "rs7213460"
# [1] "rs7213460" "1"
#
# rs3184504
# [1] "RS_number" "rs3184504" "rs597808"
# [1] "rs3184504" "1"         "0.96"
# > subset(region,rsid=="rs597808")
#         rsid          MarkerName  Somamer Effect StdErr    Pvalue  Prot        MarkerName.1  sentinel
# 626 rs597808 chr12:111973358_A_G 9188_119 0.0747 0.0136 3.797e-08 CXCL9 chr12:111884608_C_T rs3184504

# legacy code

# ieugwasr::ld_matrix() requires explicit names
# r <- sapply(1:nrow(region),function(x) {z=LDlinkR::LDpair(region[["rsid"]][x],region[["sentinel"]][x],token=key);r2=z$r2})

function panel_bqc()
{
if [ ! -d ${INF}/Fenland ]; then mkdir ${INF}/Fenland; fi
R --no-save -q <<END
  options(width=200)
  f <- file.path(Sys.getenv("proteomics_results"),"bqc19_jgh_prt_soma_list.xlsx")
  SomaLogicv4 <- openxlsx::read.xlsx(f,sheet=1,startRow=3)
  out <- file.path(Sys.getenv("INF"),"Fenland",Sys.getenv("v4"))
  write.table(SomaLogicv4,file=out,quote=FALSE,row.names=FALSE,sep="\t")
END
}

function replication_bqc()
(
  gunzip -c ${all} | \
  head -1
  join -j2 <(cut -f1,3,6 --complement --output-delimiter=' ' ${INF}/Fenland/${v4} | sed 's/-/_/' | sort -k2,2) \
           <(cut -f2,4,5,20 --output-delimiter=' ' ${INF}/work/INF1.METAL | awk '{print $2":"$3,$4,$1}' | sort -k2,2) | \
  parallel -C' ' -j20 --env proteomics_results '
    tabix ${all} {4} | grep -w {2} | grep -w {5}
  '
) > ${INF}/Fenland/replication_bqc.tsv

