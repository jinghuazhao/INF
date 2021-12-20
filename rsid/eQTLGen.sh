#!/usr/bin/bash

export eQTLGen=~/rds/public_databases/eQTLGen

function eQTLGen_tabix()
{
  export TMPDIR=${HPC_WORK}/work
  export eQTLGen=~/rds/public_databases/eQTLGen
  export eQTLGen_tabix=${eQTLGen}/tabix
# MAF
  gunzip -c ${eQTLGen}/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz | \
  bgzip -f > ${eQTLGen_tabix}/AF_AC_MAF_pos.txt.gz
  tabix -S1 -s2 -b3 -e3 -f ${eQTLGen_tabix}/AF_AC_MAF_pos.txt.gz
  ln -sf ${eQTLGen_tabix}/AF_AC_MAF_pos.txt.gz ${eQTLGen}/AF_AC_MAF_pos.txt.gz
  ln -sf ${eQTLGen_tabix}/AF_AC_MAF_pos.txt.gz.tbi ${eQTLGen}/AF_AC_MAF_pos.txt.gz.tbi
# cis/trans/cis.full
  export cis=${eQTLGen}/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
  export trans=${eQTLGen}/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz
  export cis_full=${eQTLGen}/cis-eQTLs_full_20180905.txt.gz
  export txt_gz=($cis $trans $cis_full)
  export type=(cis trans cis_full)
  for i in {0..2}
  do
    cat <(gunzip -c ${txt_gz[$i]} | head -1) <(gunzip -c ${txt_gz[$i]} | sed '1d' | sort -k3,3n -k4,4n) | \
    bgzip -f > ${eQTLGen_tabix}/${type[$i]}.txt.gz
    tabix -S1 -s3 -b4 -e4 -f ${eQTLGen_tabix}/${type[$i]}.txt.gz
    ln -sf ${eQTLGen_tabix}/${type[$i]} ${eQTLGen}/${type[$i]}.txt.gz
    ln -sf ${eQTLGen_tabix}/${type[$i]}.txt.gz.tbi ${eQTLGen}/${type[$i]}.txt.gz.tbi
  done
}

function cistrans_python()
# select named columns from .csv
{
python3 <<END
#!/usr/bin/python3
import csv 

with open('work/INF1.merge.cis.vs.trans', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        print(row['uniprot'], row['SNP'], row['p.gene'])
END
}

function lookup_merge()
{
  for cistrans in cis trans
  do
    export cistrans_rsid=${INF}/work/INF1.${cistrans}-rsid
    grep -f ${INF}/work/INF1.${cistrans} ${INF}/work/INTERVAL.rsid > ${cistrans_rsid}
    cistrans_python | \
    sort -k2,2 | \
    join -12 - <(sort -k1,1 ${cistrans_rsid}) | \
    sort -k4,4 | \
    join -14 -22 - <(cut -d' ' -f2 ${cistrans_rsid} | zgrep -f - -w ${eQTLGen}/${cistrans}.txt.gz | cut -f1,2,5,6,9 | sort -k2,2) | \
    awk '$4==$8' > ${INF}/work/eQTLGen.${cistrans}
    # all from eQTLGen but with LD
    cistrans_python | \
    sort -k3,3 | \
    join -13 -25 - <(gunzip -c ${eQTLGen}/${cistrans}.txt.gz | cut -f1,2,5,6,9 | sort -k5,5) > ${INF}/work/eQTLGen.${cistrans}-all
  done

  Rscript -e '
    ys1 <- c(paste0("Yes-",1:22),paste0("No1-",1:62238))
    ys2 <- c(paste0("Yes-",1:22),paste0("No2-",1:37))
    ys <- list(ys1,ys2)
    names(ys) <- c("eQTL","pQTL")
    INF <- Sys.getenv("INF")
    VennDiagram::venn.diagram(x = ys, filename=file.path(INF,"work","eQTLGen-cis.png"), disable.logging=TRUE, imagetype="png", output=TRUE,
                              height=12, width=12, units="cm", resolution=500,
                              fill=c("yellow","purple"), cat.pos=c(-30,30), rotation.degree = 0)
    ys1 <- c(paste0("Yes-",1:3),paste0("No1-",1:1129))
    ys2 <- c(paste0("Yes-",1:3),paste0("No2-",1:118))
    ys <- list(ys1,ys2)
    names(ys) <- c("eQTL","pQTL")
    VennDiagram::venn.diagram(x = ys, filename=file.path(INF,"work","eQTLGen-trans.png"), disable.logging=TRUE, imagetype="png", output=TRUE,
                              height=12, width=12, units="cm", resolution=500,
                              fill=c("yellow","purple"), cat.pos=c(-30,30), rotation.degree = 0)
  '
}

function lookup_jma()
{
  cd ${INF}
  module load ceuadmin/stata 
  stata <<\ \ END
  local INF : env INF
  insheet using "`INF'/sentinels/INF1.jma-rsid.cis.vs.trans.", case clear delim(" ")
  l uniprot SNP pgene
  outsheet uniprot SNP pgene using "`INF'/eQTLGen/uniprot-rsid-gene.", delim(" ") noquote noname replace
  outsheet SNP using "`INF'/eQTLGen/INF1.cis-rsid" if cistrans=="cis", noname noquote replace
  outsheet SNP using "`INF'/eQTLGen/INF1.trans-rsid" if cistrans=="trans", noname noquote replace
  END
  cd -

  for cistrans in cis trans
  do
    export cistrans_rsid=${INF}/eQTLGen/INF1.${cistrans}-rsid
  # cistrans_stata
    sort -k2,2 ${INF}/eQTLGen/uniprot-rsid-gene | \
    join -12 - <(sort -k1,1 ${cistrans_rsid}) | \
    sort -k1,1 | \
    join -11 -22 - <(cut -d' ' -f2 ${cistrans_rsid} | zgrep -f - -w ${eQTLGen}/${cistrans}.txt.gz | cut -f1,2,5,6,9 | sort -k2,2) | \
    awk '$3==$7' > ${INF}/eQTLGen/eQTLGen.${cistrans}
  # all from eQTLGen but with LD
    sort -k3,3 ${INF}/eQTLGen/uniprot-rsid-gene | \
    join -13 -25 - <(gunzip -c ${eQTLGen}/${cistrans}.txt.gz | cut -f1,2,5,6,9 | sort -k5,5) > ${INF}/eQTLGen/eQTLGen.${cistrans}-all
  done

  Rscripe -e '
    ys1 <- c(paste0("Yes-",1:37),paste0("No1-",1:81931))
    ys2 <- c(paste0("Yes-",1:37),paste0("No2-",1:62))
    ys <- list(ys1,ys2)
    names(ys) <- c("eQTL","pQTL")
    INF <- Sys.getenv("INF")
    VennDiagram::venn.diagram(x = ys, filename=file.path(INF,"eQTLGen","eQTLGen-cis.png"), imagetype="png", output=TRUE,
                              disable.logging=TRUE, height=12, width=12, units="cm", resolution=500,
                              fill=c("yellow","purple"), cat.pos=c(-30,30), rotation.degree = 0)
    ys1 <- c(paste0("Yes-",1:14),paste0("No1-",1:1462))
    ys2 <- c(paste0("Yes-",1:14),paste0("No2-",1:114))
    ys <- list(ys1,ys2)
    names(ys) <- c("eQTL","pQTL")
    VennDiagram::venn.diagram(x = ys, filename=file.path(INF,"eQTLGen","eQTLGen-trans.png"), imagetype="png", output=TRUE,
                              disable.logging=TRUE,  height=12, width=12, units="cm", resolution=500,
                              fill=c("yellow","purple"), cat.pos=c(-30,30), rotation.degree = 0)
  '
}

function run_coloc()
{
  Rscript -e '
   options(width=200)
   INF <- Sys.getenv("INF")
   M <- 1e6
   r <- as.integer(Sys.getenv("r"))
   prot <- Sys.getenv("prot")
   cvt <- Sys.getenv("cvt")
   pkgs <- c("dplyr", "ggplot2", "readr", "coloc", "GenomicRanges","seqminer")
   invisible(suppressMessages(lapply(pkgs, require, character.only=TRUE)))
   HPC_WORK <- Sys.getenv("HPC_WORK")
   gwasvcf::set_bcftools(file.path(HPC_WORK,"bin","bcftools"))
   sentinels <- subset(read.csv(cvt),cis)
   sentinel <- sentinels[r,]
   ensRegion <- with(subset(pQTLtools::inf1,prot==sentinel[["prot"]]),
                     {
                       start <- start-M
                       if (start<0) start <- 0
                       end <- end+M
                       paste0(chr,":",start,"-",end)
                     })
   gene <- subset(pQTLtools::inf1,prot==sentinel[["prot"]])[["gene"]]
   # "GWAS sumstats
   vcf <- file.path(INF,"METAL","gwas2vcf",paste0(sentinel[["prot"]],".vcf.gz"))
   gs <- gwasvcf::query_gwas(vcf, chrompos = ensRegion)
   gwas_stats <- gwasvcf::vcf_to_tibble(gs) %>%
                 filter(!is.na(ES)) %>%
                 mutate(SNP=ID,chr=as.integer(seqnames),pos=end,MAF=if_else(AF<0.5,AF,1-AF),sdY=1)
   dup1 <- duplicated(with(gwas_stats,SNP))
   # eQTLGen
   cis_eQTL <- Sys.getenv("cis_eQTL")
   eqtl_info <- seqminer::tabix.read.table(tabixFile=file.path(cis_eQTL,"AF_AC_MAF_pos.txt.gz"),tabixRange=ensRegion,stringsAsFactors=FALSE)
   names(eqtl_info) <- c("SNP","chr","pos","AlleleA","AlleleB","allA_total","allAB_total","allB_total","f")
   es <- seqminer::tabix.read.table(tabixFile=file.path(cis_eQTL,"cis_full.txt.gz"),tabixRange=ensRegion,stringsAsFactors=FALSE)
   names(es) <- c("Pvalue","SNP","chr","pos","z","A1","A2","Gene","GeneSymbol","GeneChr","GenePos","NrCohorts","N","FDR")
   es <- left_join(select(es,SNP,chr,pos,A1,A2,z,GeneSymbol,N),select(eqtl_info,SNP,chr,pos,f)) %>%
         mutate(maf=if_else(f<0.5,f,1-f),d=sqrt(2*f*(1-f)*(z^2+N)),beta=z/d,se=1/d) %>%
         select(-d,-f) %>%
         left_join(select(gwas_stats,SNP,chr,pos,REF,ALT)) %>%
         distinct() %>%
         mutate(sign=if_else(A1==ALT,1,-1),beta=sign*beta,sdY=1)
   eqtl_stats <- subset(es, GeneSymbol==gene & !is.na(beta) & !is.na(maf))
   dup2 <- duplicated(with(eqtl_stats,SNP))
   run_coloc <- function(eqtl_sumstats=eqtl_stats[!dup2,], gwas_sumstats=gwas_stats[!dup1,])
   {
     eQTL_dataset <- with(eqtl_sumstats, list(beta=beta,varbeta=se^2,N=N,MAF=maf,type="quant",snp=SNP))
     gwas_dataset <- with(gwas_sumstats, list(beta=ES,varbeta=SE^2,type="quant",snp=ID,MAF=MAF,N=SS))
     coloc_res <- coloc::coloc.abf(dataset1=eQTL_dataset, dataset2=gwas_dataset, p1=1e-4, p2=1e-4, p12=1e-5)
     res_formatted <- dplyr::as_tibble(t(as.data.frame(coloc_res$summary)))
   }
   res <- run_coloc()
   write.table(res,file=file.path(INF,"eQTLGen",paste0(r,"-",prot,"-",gene,".out")))
 ' 
}

function coloc()
{
  cd ${INF}/eQTLGen
  for r in {1..59}
  do
     export r=${r}
     export cvt=${INF}/work/INF1.merge.cis.vs.trans
     export cis_eQTL=~/rds/public_databases/eQTLGen
     read prot MarkerName < <(awk -vFS="," '$14=="cis"' ${cvt} | awk -vFS="," -vr=${r} 'NR==r{print $2,$5}')
     echo ${r} - ${prot} - ${MarkerName}
     export prot=${prot}
     export MarkerName=${MarkerName}
     run_coloc 2>&1 | \
     tee ${INF}/eQTLGen/log/${r}-${prot}-${MarkerName}.log
  done
  cat <(cat *.out | head -1 | awk '{gsub(/\"/,"",$0);print "ID","prot","gene",$0}' ) \
      <(grep -v PP.H *.out | awk '{gsub(/.out:"1"/,"",$1);split($1,a,"-");$1=a[1]" "a[2]" "a[3];print}' | sort -k1,1n) > coloc.txt
  cd -
}

if [ ! -d ${INF}/eQTLGen/log ]; then mkdir -p ${INF}/eQTLGen/log; fi

# lookup_merge
# lookup_jma

coloc
