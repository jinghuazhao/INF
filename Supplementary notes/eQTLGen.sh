#!/usr/bin/bash

if [ ! -d ${INF}/eQTLGen ]; then mkdir ${INF}/eQTLGen; fi
export eQTLGen=~/rds/public_databases/eQTLGen

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
   cvt <- Sys.getenv("cvt")
   pkgs <- c("dplyr", "ggplot2", "readr", "coloc", "GenomicRanges","seqminer")
   invisible(suppressMessages(lapply(pkgs, require, character.only=TRUE)))
   sentinels <- subset(read.csv(cvt),cis)
   sentinel <- sentinels[r,]
   ensRegion <- with(subset(pQTLtools::inf1,prot==sentinel[["prot"]]),
                     {
                       start <- start-M
                       if (start<0) start <- 0
                       end <- end+M
                       paste0(chr,":",start,"-",end)
                     })
   ensGene <- subset(pQTLtools::inf1,prot==sentinel[["prot"]])[["gene"]]
   HPC_WORK <- Sys.getenv("HPC_WORK")
   gwasvcf::set_bcftools(file.path(HPC_WORK,"bin","bcftools"))
   cat("GWAS sumstats\n")
   vcf <- file.path(INF,"METAL","gwas2vcf",paste0(sentinel[["prot"]],".vcf.gz"))
   gwas_stats <- gwasvcf::query_gwas(vcf, chrompos = ensRegion)
   gwas_stats <- gwasvcf::vcf_to_granges(gwas_stats)
   # eQTLGen
   safe_import <- purrr::safely(import_eQTLCatalogue)
   summary_list <- purrr::map(ftp_path_list, ~safe_import(., region38, selected_gene_id = ensGene, column_names))
   result_list <- purrr::map(summary_list, ~.$result)
   result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
   result_filtered <- purrr::map(result_list[lapply(result_list,nrow)!=0], ~dplyr::filter(., !is.na(se)))
   purrr::map_df(result_filtered, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")
   coloc.abf()
   # etc
 ' 
}

function eQTLGen_tabix()
{
  export TMPDIR=${HPC_WORK}/work
  export eQTLGen=~/rds/public_databases/eQTLGen
  export eQTLGen_tabix=${eQTLGen}/tabix
# MAF
  gunzip -c ${eQTLGen}/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz | \
  bgzip -f > ${eQTLGen_tabix}/AF_AC_MAF_pos.txt.gz
  tabix -S1 -s2 -b3 -e3 -f ${eQTLGen_tabix}/AF_AC_MAF_pos.txt.gz
# cis/trans/cis.full
  export cis=${eQTLGen}/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
  export trans=${eQTLGen}/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz
  export cis.full=${eQTLGen}/cis-eQTLs_full_20180905.txt.gz
  export txt_gz=($cis $trans $cis.full)
  export type=(cis trans cis.full)
  for i in {0..2}
  do
    cat <(gunzip -c ${txt_gz[$i]} | head -1) <(gunzip -c ${txt_gz[$i]} | sed '1d' | sort -k3,3n -k4,4n) | \
    bgzip -f > ${eQTLGen_tabix}/${type[$i]}.txt.gz
    tabix -S1 -s3 -b4 -e4 -f ${eQTLGen_tabix}/${type[$i]}.txt.gz
    ln -sf ${eQTLGen_tabix}/${type[$i]} ${eQTLGen}/${type[$i]}.txt.gz
    ln -sf ${eQTLGen_tabix}/${type[$i]}.txt.gz.tbi ${eQTLGen}/${type[$i]}.txt.gz.tbi
  done
}

function coloc()
{
for r in {1..59}
do
   export r=${r}
   export cvt=${INF}/work/INF1.merge.cis.vs.trans
   read prot MarkerName < \
                        <(awk -vFS="," '$14=="cis"' ${cvt} | \
                          awk -vFS="," -vr=${r} 'NR==r{print $2,$5}')
   echo ${r} - ${prot} - ${MarkerName}
   export prot=${prot}
   export MarkerName=${MarkerName}
   if [ ! -f ${INF}/coloc/${prot}-${MarkerName}.pdf ] || \
      [ ! -f ${INF}/coloc/${prot}-${MarkerName}.RDS ]; then
     cd ${INF}/eQTLGen
#    R --no-save < ${INF}/rsid/coloc.R 2>&1 | \
#    tee ${prot}-${MarkerName}.log
#    ls *tbi | xargs -I {} bash -c "rm {}"
     cd -
   fi
done
}
