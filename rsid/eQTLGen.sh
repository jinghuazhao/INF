#!/usr/bin/bash

export eQTLGen=~/rds/public_databases/eQTLGen
export TMPDIR=${HPC_WORK}/work

function eQTLGen_tabix()
{
  export eQTLGen_tabix=tabix
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
  cd ${eQTLGen}
  for i in {0..2}
  do
    cat <(gunzip -c ${txt_gz[$i]} | head -1) <(gunzip -c ${txt_gz[$i]} | sed '1d' | sort -k3,3n -k4,4n) | \
    bgzip -f > ${eQTLGen_tabix}/${type[$i]}.txt.gz
    tabix -S1 -s3 -b4 -e4 -f ${eQTLGen_tabix}/${type[$i]}.txt.gz
    ln -sf ${eQTLGen_tabix}/${type[$i]}.txt.gz
    ln -sf ${eQTLGen_tabix}/${type[$i]}.txt.gz.tbi
  done
  cd -
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
  export jma=${INF}/sentinels/INF1.jma-rsid.cis.vs.trans
  cd ${INF}
  module load ceuadmin/stata 
  stata <<\ \ END
  local jma : env jma
  insheet using "`jma'.", case clear delim(" ")
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
  # full details
    export f=~/rds/public_databases/eQTLGen/${cistrans}.txt.gz
    cat ${INF}/eQTLGen/eQTLGen.${cistrans} | \
    parallel -C' ' --env INF --env f '
      grep -w {1} <(Rscript -e "subset(read.table(Sys.getenv(\"jma\"),header=TRUE),select=c(uniprot,SNP,p.gene,Chr,bp))") | \
      grep -w {2} | \
      awk "{cmd=sprintf(\"tabix %s %d:%d-%d\", ENVIRON[\"f\"], \$5,\$6,\$6);system(cmd)}" | \
      awk -v gene={3} "gene==\$9"
    ' | \
    sort -k3,3n -k4,4n > ${INF}/eQTLGen/eQTLGen.${cistrans}-detailed
  done

  Rscript -e '
    suppressMessages(library(dplyr))
    library(seqminer)
    INF <- Sys.getenv("INF")
    jma <- Sys.getenv("jma")
    eQTLGen <- Sys.getenv("eQTLGen")
    cis_full_cols <- c("P_eqtl","SNP","chr","pos","z","A1","A2","Gene","GeneSymbol","GeneChr","GenePos","NrCohorts","N","FDR")
    cols <- c("P_eqtl","SNP","chr","pos","Allele1","Allele2","z","Gene","GeneSymbol","GeneChr","GenePos","NrCohorts","N","FDR","BonferroniP")
    reverse_strand <- function(allele,forward=c("A","C","G","T"),reverse=c("T","G","C","A")) reverse[(1:4)[toupper(allele)==forward]]
    for(cistrans in c("cis","trans"))
    {
       qtls <- read.table(file.path(INF,"eQTLGen",paste0("eQTLGen.",cistrans)),
                          col.names=c("rsid","uniprot","gene","p1","A1","A2","gene2")) %>%
               left_join(read.table(jma,header=TRUE),by=c("uniprot"="uniprot","gene"="p.gene")) %>%
               rename(chr=Chr,pos=bp) %>%
               left_join(read.table(file.path(INF,"eQTLGen",paste0("eQTLGen.",cistrans,"-detailed")),col.names=cols),
                         by=c('chr'='chr','pos'='pos','gene'='GeneSymbol','SNP'='SNP')) %>%
               filter(!is.na(z)) %>%
               select(chr,pos,prot,gene,SNP,b,se,p,bJ,bJ_se,pJ,z,P_eqtl,N,A1,A2,Allele1,Allele2) %>%
               distinct()
       for(i in 1:nrow(qtls))
       {
         snpRegion <- with(qtls[i,],paste0(chr,":",pos,"-",pos))
         eqtl_info <- tabix.read.table(tabixFile=file.path(eQTLGen,"AF_AC_MAF_pos.txt.gz"),tabixRange=snpRegion,stringsAsFactors=FALSE)
         names(eqtl_info) <- c("SNP","chr","pos","AlleleA","AlleleB","allA_total","allAB_total","allB_total","f")
         d <- left_join(qtls[i,],eqtl_info)
         qtls[i,c("b_eqtl","se_eqtl")] <- with(d,signif(gap::get_b_se(f,N,z),2))
         qtls[i,"cis.trans"] <- cistrans
       }
       qtls <- filter(qtls,A1==Allele1) %>%
               select(-A1,-A2,-Allele1,-Allele2)
       assign(paste0(cistrans,".pqtl"),arrange(qtls,chr,pos,prot))
    }
    pqtls <- bind_rows(cis.pqtl,trans.pqtl)
    write.table(pqtls,file=file.path(INF,"eQTLGen","SCALLOP-INF-eQTLGen.txt"),quote=FALSE,row.names=FALSE)
  '
#              mutate(A1A2=if_else(A1<A2,paste0(A1,A2),paste0(A2,A1)),
#                     rA1=reverse_strand(A1),rA2=reverse_strand(A2),rA1A2=if_else(rA1<rA2,paste0(rA1,rA2),paste0(rA2,rA1)),
#                     Allele1Allele2=if_else(Allele1<Allele2,paste0(Allele1,Allele2),paste0(Allele2,Allele1))) %>%
#              filter((A1A2==Allele1Allele2 | rA1A2==Allele1Allele2))
#   qtls <- mutate(qtls,be=if_else(A1A2==Allele1Allele2 & A1==Allele1,be,-be),be=if_else(rA1A2=Allele1Allele2 & rA1==Allele1,-be,be))

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
   ensRegion <- with(subset(pQTLdata::inf1,prot==sentinel[["prot"]]),
                     {
                       start <- start-M
                       if (start<0) start <- 0
                       end <- end+M
                       paste0(chr,":",start,"-",end)
                     })
   gene <- subset(pQTLdata::inf1,prot==sentinel[["prot"]])[["gene"]]
   # "GWAS sumstats
   vcf <- file.path(INF,"METAL","gwas2vcf",paste0(sentinel[["prot"]],".vcf.gz"))
   gs <- gwasvcf::query_gwas(vcf, chrompos = ensRegion)
   gwas_stats <- gwasvcf::vcf_to_tibble(gs) %>%
                 filter(!is.na(ES)) %>%
                 mutate(SNP=ID,chr=as.integer(seqnames),pos=end,MAF=if_else(AF<0.5,AF,1-AF),sdY=1,
                        snpid=gap::chr_pos_a1_a2(chr,pos,ALT,REF))
   # eQTLGen
   cis_eQTL <- Sys.getenv("cis_eQTL")
   eqtl_info <- seqminer::tabix.read.table(tabixFile=file.path(cis_eQTL,"AF_AC_MAF_pos.txt.gz"),
                                           tabixRange=ensRegion,stringsAsFactors=FALSE) %>%
                setNames(c("SNP","chr","pos","AlleleA","AlleleB","allA_total","allAB_total","allB_total","f")) %>%
                filter(!is.na(AlleleA) & !is.na(AlleleB)) %>%
                mutate(snpid=gap::chr_pos_a1_a2(chr,pos,AlleleA,AlleleB))
   es <- seqminer::tabix.read.table(tabixFile=file.path(cis_eQTL,"cis_full.txt.gz"),tabixRange=ensRegion,stringsAsFactors=FALSE) %>%
         setNames(c("Pvalue","SNP","chr","pos","z","A1","A2","Gene","GeneSymbol","GeneChr","GenePos","NrCohorts","N","FDR")) %>%
         mutate(snpid=gap::chr_pos_a1_a2(chr,pos,A1,A2)) %>%
         select(SNP,chr,pos,A1,A2,z,GeneSymbol,N,snpid) %>%
         filter(GeneSymbol==gene) %>%
         left_join(select(eqtl_info,SNP,chr,pos,f,snpid)) %>%
         mutate(maf=if_else(f<0.5,f,1-f),d=sqrt(2*f*(1-f)*(z^2+N)),beta=z/d,se=1/d) %>%
         select(-d,-f) %>%
         left_join(select(gwas_stats,SNP,chr,pos,REF,ALT,snpid)) %>%
         distinct() %>%
         mutate(sign=if_else(A1==ALT,1,-1),beta=sign*beta,sdY=1)
   eqtl_stats <- subset(es, GeneSymbol==gene & !is.na(beta) & !is.na(maf))
   dup1 <- duplicated(with(gwas_stats,snpid))
   dup2 <- duplicated(with(eqtl_stats,snpid))
   run_coloc <- function(gwas_sumstats=gwas_stats[!dup1,],eqtl_sumstats=eqtl_stats[!dup2,])
   {
     eQTL_dataset <- with(eqtl_sumstats, list(beta=beta,varbeta=se^2,N=N,MAF=maf,type="quant",snp=snpid))
     gwas_dataset <- with(gwas_sumstats, list(beta=ES,varbeta=SE^2,type="quant",snp=snpid,MAF=MAF,N=SS))
     coloc_res <- coloc::coloc.abf(dataset1=gwas_dataset, dataset2=eQTL_dataset, p1=1e-4, p2=1e-4, p12=1e-5)
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

function cis_lst()
{
  export M=250000
  Rscript -e '
    suppressMessages(library(dplyr))
    INF <- Sys.getenv("INF")
    M <- Sys.getenv("M") %>% as.integer
    INF1_METAL <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
                  left_join(select(pQTLdata::inf1,prot,gene,chr,start,end,target.short)) %>%
                  filter(cis.trans=="cis") %>%
                  mutate(start=if_else(start-M<0,0,start-M),end=end+M,region=paste0(chr,":",start,"-",end)) %>%
                  select(uniprot,prot,gene,region,chr,rsid,target.short)
    write.table(INF1_METAL,file=file.path(INF,"eQTLGen","cis.lst"),row.names=FALSE,col.names=FALSE,quote=FALSE)
  '
}

function lz()
{
  export cis_full=${eQTLGen}/cis-eQTLs_full_20180905.txt.gz
  export eQTLGen_tabix=${eQTLGen}/tabix
  module load python/2.7
  export dir=${INF}/eQTLGen
# 500kb for CCL4
# echo P13236 CCL4 CCL4 17 33930983 34933014 rs8064426 CCL4 | \
  cut -d ' ' -f1-4,6,7 ${dir}/cis.lst | awk '{split($4,a,":|-");print $1,$2,$3,a[1],a[2],a[3],$5,$6}' | \
  parallel -j8 -C ' ' --env dir '
  echo "{1}-{2}-{3}"
# eQTLGen
  cat <(echo -e "snpid rsid chr pos a1 a2 mlog10p") \
      <(tabix ${eQTLGen_tabix}/cis_full.txt.gz {4}:{5}-{6} | \
        awk -vgene={3} "\$9==gene' \| \
        cut -f2-7 | \
        awk "
        {
          if (\$5<\$6) snpid=\"chr\"\$2\":\"\$3\"_\"\$5\"_\"\$6;
          else snpid=\"chr\"\$2\":\"\$3\"_\"\$6\"_\"\$5
          print snpid, \$1, \$2, \$3, \$5, \$6, \$4
        }") | \
  Rscript -e "
     suppressMessages(library(dplyr))
     z <- within(read.table(\"stdin\",header=TRUE),{mlog10p <- -gap::log10p(mlog10p);mlog10p=if_else(mlog10p>20,20,mlog10p)});
     write.table(z,row.names=FALSE,quote=FALSE,sep=\"\t\")
  " | \
  gzip -f > ${dir}/eQTLGen-{1}-{2}-{3}.tsv.gz
# SCALLOP/INF
  cat <(echo -e "snpid rsid chr pos a1 a2 mlog10p") \
      <(tabix ${INF}/METAL/{2}-1.tbl.gz {4}:{5}-{6} | \
        awk "
        {
          split(\$3,a,\"_\")
          print a[1],\$1,\$2,\$10/\$11,\$3,toupper(\$4),toupper(\$5)
        }" | \
        sort -k1,1 | \
        join -12 -21 <(grep chr{4} ${INF}/work/snp_pos) - | \
        awk -vOFS="\t" "{print \$6, \$2, \$3, \$4, \$7, \$8, \$5}") | \
  Rscript -e "
     z <- within(read.table(\"stdin\",header=TRUE),{mlog10p <- -gap::log10p(mlog10p)});
     write.table(z,row.names=FALSE,quote=FALSE,sep=\"\t\")
  " | \
  gzip -f > ${dir}/INF-{1}-{2}-{3}.tsv.gz
  export nlines=$(gunzip -c ${dir}/eQTLGen-{1}-{2}-{3}.tsv.gz | wc -l | cut -d" " -f1)
  if [ ${nlines} -eq 1 ]; then
    echo -e "{1}\t{2}\t{3}\tremoved"
    rm ${dir}/eQTLGen-{1}-{2}-{3}.tsv.gz ${dir}/INF-{1}-{2}-{3}.tsv.gz
  else
    (
      echo -e "chr\tpos\trsid\tmlog10P"
      gunzip -c ${dir}//eQTLGen-{1}-{2}-{3}.tsv.gz | \
      awk -v OFS="\t" "NR>1 {print \$3,\$4,\$2,\$7}" | \
      sort -k1,1n -k2,2n
    ) > ${dir}/eQTLGen-{1}-{2}-{3}.lz
    locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${dir}/eQTLGen-{1}-{2}-{3}.lz \
              --delim tab title="eQTLGen: {8} ({3})-{7}" \
              --markercol rsid --pvalcol mlog10P --no-transform --chr {4} --start {5} --end {6} --cache None \
              --no-date --plotonly --prefix=eQTLGen-{3} --rundir ${dir} --refsnp {7}
    (
      echo -e "chr\tpos\trsid\tmlog10P"
      gunzip -c ${INF}/eQTLGen/INF-{1}-{2}-{3}.tsv.gz | \
      awk -v OFS="\t" "NR>1 {print \$3,\$4,\$2,\$7}" | \
      sort -k1,1n -k2,2n
    ) > ${dir}/INF-{1}-{2}-{3}.lz
    locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${dir}/INF-{1}-{2}-{3}.lz \
              --delim tab title="SCALLOP: {8} ({3})-{7}" \
              --markercol rsid --pvalcol mlog10P --no-transform --chr {4} --start {5} --end {6} --cache None \
              --no-date --plotonly --prefix=INF-{3} --rundir ${dir} --refsnp {7}
    pdftopng -f 1 -l 1 -r 300 ${dir}/eQTLGen-{3}_{7}.pdf ${dir}/eQTLGen-{3}_{7}
    pdftopng -f 1 -l 1 -r 300 ${dir}/INF-{3}_{7}.pdf ${dir}/INF-{3}_{7}
    mv ${dir}/eQTLGen-{3}_{7}-000001.png ${dir}/eQTLGen-{3}_{7}.png
    mv ${dir}/INF-{3}_{7}-000001.png ${dir}/INF-{3}_{7}.png
    convert -append ${dir}/eQTLGen-{3}_{7}.png ${dir}/INF-{3}_{7}.png -resize x500 -density 300 ${dir}/{1}-{2}-{3}.png
    convert ${dir}/{1}-{2}-{3}.png -quality 0 ${dir}/{1}-{2}-{3}.jp2
    convert ${dir}/{1}-{2}-{3}.jp2 ${dir}/{2}-{1}-{3}-lz.pdf
    rm ${dir}/eQTLGen-{1}-{2}-{3}.lz ${dir}/INF-{1}-{2}-{3}.lz
    rm ${dir}/eQTLGen-{3}_{7}.pdf ${dir}/INF-{3}_{7}.pdf
    rm ${dir}/eQTLGen-{3}_{7}.png ${dir}/INF-{3}_{7}.png
    rm ${dir}/{1}-{2}-{3}.jp2 ${dir}/{1}-{2}-{3}.png
  fi
  '
  qpdf --empty --pages $(ls ${dir}/*-lz.pdf) -- ${dir}/eQTLGen-protein-lz-20.pdf
# rm ${dir}/*-lz.pdf
}
# ls -l eQTLGen/eQTLGen* -S | awk -vOFS="\t" '$(NF-4)==51 {split($NF,a,"-");split(a[3],b,".");print a[2],a[3],b[1]}' | xsel -i
# convert -density 300 ${dir}/eQTLGen-{2}_{7}.pdf[0] ${dir}/eQTLGen-{2}_{7}.png
# convert -density 300 ${dir}/INF-{2}_{7}.pdf[0] ${dir}/INF-{2}_{7}.png
# gs -sDEVICE=jpeg -r300 -dNOPAUSE -dBATCH -sOutputFile=eQTLGen-{2}_{7}.jpg -dFirstPage=1 -dLastPage=1 eQTLGen-{2}_{7}.pdf
# gs -sDEVICE=jpeg -r300 -dNOPAUSE -dBATCH -sOutputFile=INF-{2}_{7}.jpg -dFirstPage=1 -dLastPage=1 INF-{2}_{7}.pdf

if [ ! -d ${INF}/eQTLGen/log ]; then mkdir -p ${INF}/eQTLGen/log; fi

# lookup_merge
# lookup_jma
# coloc
lz

function check_lz()
{
  cd ${INF}/eQTLGen
  join -11 -22 -v2 \
       <(ls *pdf | grep -v protein | sed 's/-lz.pdf//;s/-/\t/g' | sort -k1,1) \
       <(ls eQTLGen*tsv.gz | sed 's/eQTLGen-//;s/-/\t/g;s/.tsv.gz//'| sort -k2,2)
}

function addition()
{
# 20 22 25 34
cat << 'EOL' > ll
CCL23 CCL23
CCL4 CCL4
IL.10 IL10
MIP.1.alpha CCL3
EOL
(
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
  done
) | grep -f <(cut -d' ' -f1 ll) - -w
}
