#!/usr/bin/bash

export prot=TNFB
export chr=12
export start=6300000
export end=6700000
export region=${chr}:${start}-${end}
export dir=~/rds/results/public/gwas/blood_cell_traits/chen_2020
export TMPDIR=/rds/user/jhz22/hpc-work/work
export v3=~/rds/results/public/gwas/multiple_sclerosis/discovery_metav3.0.meta.gz

function lz()
{
# gunzip -c discovery_metav3.0.meta.gz | grep rs1800693
# CHR BP SNP A1 A2 N P OR
# 12 6440009 rs1800693 T C 14 1.017e-13 0.8808
# gunzip -c Final-metaanalysis-echip.txt.gz | grep 1800693
#  CHR         BP            SNP  A1  A2   N           P        P(R)      OR   OR(R)       Q       I
#  12     6440009  exm-rs1800693   T   C  13   1.003e-26    2.11e-09  1.0332  1.0300  0.0351   46.00
  module load python/2.7
  ls ${dir}/tsv/*gz | xargs -I{} basename {} .gz | \
  grep -e wbc -e mono -e neut -e lymph -e eo -e baso | \
  parallel -C' ' --env INF --env dir --env region --env chr --env start --env end  '
  (
    echo snpid rsid chromosome base_pair_location effect_allele other_allele effect_allele_frequency beta standard_error p_value
    join -23 ${INF}/work/INTERVAL.rsid <(tabix ${dir}/tsv/{}.gz ${region}) | \
    sort -k3,3n -k4,4n
  ) | \
  awk "!index(\$2,\":\")" | \
  tr " " "\t" > ${INF}/MS/{}.lz
  rm -rf ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${INF}/MS/{}.lz --delim tab title="{}" \
            --markercol rsid --pvalcol p_value --chr ${chr} --start ${start} --end ${end} \
            --no-date --plotonly --prefix="{}" --rundir .
  qpdf {}_chr${chr}_${start}-${end}.pdf --pages . 1 -- {}-lz.pdf
  '
  qpdf --empty --pages *lz.pdf -- blood-cell-traits.pdf
  (
    echo -e "SNPid\tSNP\tchr\tpos\ta1\ta2\tb\tse\tp"
    Rscript -e 'require(dplyr)
                write.table(read.table(Sys.getenv("v3"),header=TRUE) %>%
                            mutate(A1=toupper(A1),A2=toupper(A2),
                                   SNP=paste0("chr",CHR,":",BP,"_",if_else(A1>A2,paste0(A2,"_",A1),paste0(A1,"_",A2))),
                                   b=log(OR),se=TwoSampleMR::get_se(b,P)) %>%
                            filter(CHR==as.integer(Sys.getenv("chr")) & BP>=as.integer(Sys.getenv("start")) & BP < as.integer(Sys.getenv("end"))) %>%
                            filter(!is.na(b) & !is.na(se) & b!=0 & se!=0) %>%
                            arrange(BP) %>%
                            select(CHR,BP,SNP,A1,A2,b,se,P),
                            col.names=FALSE,row.names=FALSE,quote=FALSE)' | \
    sort -k1,1 | \
    join -23 ${INF}/work/INTERVAL.rsid - | \
    tr ' ' '\t'
  ) > ${INF}/MS/v3.lz
  rm -rf ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${INF}/MS/v3.lz --delim tab title="MS" \
            --markercol SNP --pvalcol p --chr ${chr} --start ${start} --end ${end} \
            --no-date --plotonly --prefix="v3" --rundir .
  qpdf v3_chr${chr}_${start}-${end}.pdf --pages . 1 -- v3-lz.pdf
}

function ma()
{
  echo wbc mono neut lymph eo baso | \
  tr ' ' '\n' | \
  parallel -C' ' --env INF '
  (
    echo SNP A1 A2 ref b se p N
    awk -vN=562132 "NR>1 && \$10 <= 5e-8 {print \$1,\$5,\$6,\$7,\$8,\$9,\$10,N}" ${INF}/MS/EUR-{}.lz
  ) > ${INF}/MS/EUR-{}.ma
  sed "1d" ${INF}/MS/EUR-{}.lz | \
  cut -f1 > ${INF}/MS/EUR-{}.snpid
  sed "1d" ${INF}/MS/EUR-{}.lz | \
  sort -k10,10g | \
  awk "NR==1{print \$1}" > ${INF}/MS/EUR-{}.top
  '
  (
    echo SNP A1 A2 freq b se p N
    awk 'NR>1{print $1,$5,$6,$7,$8,$9}' ${INF}/MS/v3.lz | \
    join - <(awk 'NR>1{print $2,$3,$4,$5}' ${INF}/work/LTBR.frq) | \
    awk -v N=115803 '
    {
      if($1==$7) freq=$9; else freq=1-$9
      print $1,$2,$3,freq,$4,$5,$6,N
    }'
  ) > ${INF}/MS/EUR-v3.ma
  sed "1d" ${INF}/MS/EUR-v3.ma | \
  cut -d" " -f1 > ${INF}/MS/EUR-v3.snpid
  sed "1d" ${INF}/MS/EUR-v3.ma | \
  sort -k7,7g | \
  awk "NR==1{print \$1}" > ${INF}/MS/EUR-v3.top}

function cojo()
{
  module load plink/2.00-alpha
  for cell in wbc mono neut lymph eo baso v3
  do
      plink2 --bfile ${INF}/INTERVAL/cardio/INTERVAL --extract ${INF}/MS/EUR-${cell}.snpid \
             --geno 0.1 --mind 0.1 --maf 0.005 --indep-pairwise 1000kb 1 0.1 --out ${INF}/MS/EUR-${cell}
      if [ $(grep -w -f ${INF}/MS/EUR-${cell}.top ${INF}/MS/EUR-${cell}.prune.in | wc -l) -eq 0 ]; then
         export top=$(cat ${INF}/MS/EUR-${cell}.top)
         export i=$(grep -w -f ${INF}/MS/EUR-${cell}.prune.in ${INF}/INTERVAL/cardio/INTERVAL.bim | \
                    awk 'function abs(x) {if (x<0) return -x; else return x;} {
                         top=ENVIRON["top"];split(top,a,":");split(a[2],b,"_")
                         d=abs($4-b[1]);
                         print $1, $2, $4, d}' | \
                    sort -r -k4,4n | \
                    awk 'NR==1 {print $2}' \
                   )
         sed -i 's/'"$i"'/'"$top"'/g' ${INF}/MS/EUR-${cell}.prune.in
      fi
      sort ${INF}/MS/EUR-${cell}.prune.in > ${INF}/MS/EUR-${cell}.prune
      export P_threshold=5e-8
    # drop the --extract option with v3
      if [ "${cell}" == "v3" ]; then
      gcta-1.9 --bfile ${INF}/INTERVAL/cardio/INTERVAL \
               --cojo-file ${INF}/MS/EUR-${cell}.ma \
               --cojo-slct \
               --cojo-p ${P_threshold} \
               --maf 0.005 \
               --cojo-collinear 0.9 \
               --out ${INF}/MS/EUR-${cell}
      echo chr12:6440009_C_T ${INF}/MS/rs1800693/EUR-v3.top
      gcta-1.9 --bfile ${INF}/INTERVAL/cardio/INTERVAL \
               --cojo-file ${INF}/MS/EUR-${cell}.ma \
               --cojo-cond ${INF}/MS/rs1800693/EUR-${cell}.top \
               --cojo-p ${P_threshold} \
               --maf 0.005 \
               --cojo-collinear 0.9 \
               --out ${INF}/MS/rs1800693/EUR-${cell}
      echo chr12:6514963_A_C ${INF}/MS/rs2354485/EUR-v3.top
      gcta-1.9 --bfile ${INF}/INTERVAL/cardio/INTERVAL \
               --cojo-file ${INF}/MS/EUR-${cell}.ma \
               --cojo-cond ${INF}/MS/rs2364485/EUR-${cell}.top \
               --cojo-p ${P_threshold} \
               --maf 0.005 \
               --cojo-collinear 0.9 \
               --out ${INF}/MS/rs2364485/EUR-${cell}
      else
      gcta-1.9 --bfile ${INF}/INTERVAL/cardio/INTERVAL \
               --cojo-file ${INF}/MS/EUR-${cell}.ma \
               --extract ${INF}/MS/EUR-${cell}.prune \
               --cojo-slct \
               --cojo-p ${P_threshold} \
               --maf 0.005 \
               --cojo-collinear 0.9 \
               --out ${INF}/MS/EUR-${cell}
      gcta-1.9 --bfile ${INF}/INTERVAL/cardio/INTERVAL \
               --cojo-file ${INF}/MS/EUR-${cell}.ma \
               --extract ${INF}/MS/EUR-${cell}.prune \
               --cojo-cond ${INF}/MS/EUR-${cell}.top \
               --cojo-p ${P_threshold} \
               --maf 0.005 \
               --cojo-collinear 0.9 \
               --out ${INF}/MS/EUR-${cell}
      fi
  done
}

cojo

function blood_cell_traits()
{
  export ext=_EUR_buildGRCh37.tsv.gz
  ls ${dir}/*EUR* | xargs -I{} basename {} ${ext} | \
  parallel -C' ' --env dir --env ext '
  (
    echo chromosome base_pair_location snpid effect_allele other_allele effect_allele_frequency beta standard_error p_value | \
    tr " " "\t"
    gunzip -c ${dir}/{}${ext} | \
    sed "1d" | \
    sort -k1,1n -k2,2n | \
    awk -vOFS="\t" "{
      if (\$3<\$4) snpid=\"chr\"\$1\":\"\$2\"_\"\$3\"_\"\$4;
      else snpid=\"chr\"\$1\":\"\$2\"_\"\$4\"_\"\$3;
      print \$1, \$2, snpid, \$3, \$4, \$5, \$6, \$7, \$8
    }"
  ) | \
  bgzip -f > ${dir}/tsv/EUR-{}.gz
  tabix -f -S1 -s1 -b2 -e2 ${dir}/tsv/EUR-{}.gz
  '
}

blood_cell_traits

function coloc()
{
  join <(sed '1d' ${INF}/MS/EUR-v3.ma | awk '{print $1,$2,$3,$4,$5,$6,$8}' | sort -k1,1) \
       <(Rscript -e '
              suppressMessages(library(dplyr))
              cis_pQTL <- merge(read.delim("eQTLGen.lz") %>% filter(GeneSymbol=="LTBR"),read.delim("eQTLGen.AF"),by="SNP") %>%
              mutate(data.frame(gap::get_b_se(AlleleB_all,NrSamples,Zscore)))
              write.table(cis_pQTL[c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","AlleleB_all","b","se","NrSamples")], sep = "\t",
                          row.names = FALSE, col.names = TRUE, quote=FALSE)
             ' | \
         awk '{
                 chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                 if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                 print snpid,A1,A2,$6,$7,$8,$9
              }' | \
         sort -k1,1 \
         ) | \
  join - <(sed '1d' ${INF}/work/TNFB.lz | \
           awk '{
                   chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                   if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                   print snpid,A1,A2,$6,$10,$11,$18,$1
                }' | \
           sort -k1,1 \
          ) | \
  awk -vOFS="\t" '
  {
    if($2!=$8) $11=-$11
    if($2!=$14) $17=-$17
    print
  }' | \
  awk 'a[$1]++==0' | \
  awk 'a[$20]++==0' > ${INF}/work/LTBR.coloc
  Rscript -e '
    options(width=200)
    INF <- Sys.getenv("INF")
    MS_LTBR_TNFB <- read.table(file.path(INF,"work","LTBR.coloc"))
    MS <- MS_LTBR_TNFB[c(1:7,20)]
    LTBR <- MS_LTBR_TNFB[c(1,8:13,20)]
    TNFB <- MS_LTBR_TNFB[c(1,14:20)]
    names(MS) <- names(LTBR) <- names(TNFB) <- c("snpid","A1","A2","MAF","beta","se","N","rsid")
    MS <- within(MS,{varbeta=se^2})
    LTBR <- within(LTBR,{varbeta=se^2})
    TNFB <- within(TNFB,{varbeta=se^2})
    MS <- c(as.list(MS),type="quant")
    LTBR <- c(as.list(LTBR),type="quant")
    TNFB <- c(as.list(TNFB),type="quant")
    require(coloc)
    abf12 <- coloc.abf(MS,LTBR)
    abf13 <- coloc.abf(MS,TNFB)
    abf23 <- coloc.abf(LTBR,TNFB)
  '
}

function gsmr()
{
# MS
  cat ${INF}/MS/EUR-v3.ma | gzip -f > ${INF}/MS/gsmr_MS.ma.gz
# prot
  (
    echo SNP A1 A2 freq b se p N
    zcat ${INF}/METAL/${prot}-1.tbl.gz | \
    awk -vchr=${chr} -vstart=${start} -vend=${end} 'NR>1 && $1==chr && $2 >= start && $2 < end {print $3,toupper($4),toupper($5),$6,$10,$11,10^$12,$18}'
  ) | gzip -f > ${INF}/MS/gsmr_${prot}.ma.gz
# control files
  if [ ! -f ${INF}/MS/gsmr_ref_data ]; then echo ${INF}/work/INTERVAL > ${INF}/MS/gsmr_ref_data; fi
  if [ ! -f ${INF}/MS/gsmr_MS ]; then echo MS ${INF}/MS/gsmr_MS.ma.gz > ${INF}/MS/gsmr_MS; fi
  if [ ! -f ${INF}/MS/gsmr_${prot} ]; then echo ${prot} ${INF}/MS/gsmr_${prot}.ma.gz > ${INF}/MS/gsmr_${prot}; fi

  gcta-1.9 --mbfile ${INF}/MS/gsmr_ref_data --gsmr-file ${INF}/MS/gsmr_MS ${INF}/MS/gsmr_${prot} \
           --gsmr-direction 0 \
           --clump-r2 0.05 --gwas-thresh 1e-5 --diff-freq 0.4 --heidi-thresh 0.05 --gsmr-snp-min 10 --effect-plot \
           --out ${INF}/MS/gsmr_${prot}-MS

  Rscript -e '
    prot <- Sys.getenv("prot")
    source("http://cnsgenomics.com/software/gcta/res/gsmr_plot.r")
    gsmr_data <- read_gsmr_data(paste0("MS/gsmr_",prot,"-MS.eff_plot.gz"))
    gsmr_summary(gsmr_data)
    pdf(paste0("MS/gsmr_",prot,"-MS.eff_plot.pdf"))
    plot_gsmr_effect(gsmr_data, prot, "MS", colors()[75])
    dev.off()
  '
}

function wgs()
{
# all significant SNPs
  (
    zcat ${INF}/METAL/${prot}-1.tbl.gz | head -1
    zcat ${INF}/METAL/${prot}-1.tbl.gz | \
    awk 'NR>1 && length($4)==1 && length($5)==1 && $12<=-5' | \
    sort -k1,1n -k2,2n | \
    gzip -f > ${INF}/MS/${prot}.p.gz
# HLA
    zcat ${INF}/MS/${prot}.p.gz | \
    awk 'NR>1 && !($1 == 6 && $2 >= 25392021 && $2 < 33392022)'
    zcat ${INF}/MS/${prot}.p.gz | \
    awk '$1 == 6 && $2 >= 25392021 && $2 < 33392022' | \
    sort -k12,12g | \
    awk 'NR==1'
  ) > ${INF}/MS/${prot}.p
  export lines=$(expr $(wc -l ${INF}/MS/${prot}.p | cut -d' ' -f1) - 1)
  if [ $lines -eq 0 ]; then rm ${INF}/MS/${prot}.p; fi
  (
    awk -vprot=${prot} -vOFS="\t" 'BEGIN{print prot, "rsid", "chr", "pos", "beta", "se", "snpid", "A1", "A2", "EAF", "P", "N"}'
    awk -vFS="\t" '{split($3,a,"_"); print a[1],$1,$2,$10,$11,$3,toupper($4),toupper($5),$6,-$12,$18}' ${INF}/MS/${prot}.p | \
    sort -k1,1 | \
    join -12 -21 ${INF}/work/snp_pos - | \
    awk 'a[$7]++==0' | \
    awk -vprot=${prot} -vOFS="\t" '{print prot, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' | \
    sort -k3,3n -k4,4n
  ) > ${INF}/MS/${prot}-pQTL.dat
  R --no-save -q <<\ \ END
    options(echo=FALSE,width=200)
    INF <- Sys.getenv("INF")
    prot <- Sys.getenv("prot")
    rsid <- Sys.getenv("rsid")
    pQTL <- file.path(INF,"MS",paste0(prot,"-pQTL.dat"))
    library(TwoSampleMR)
    x <- read_exposure_data(pQTL,
         clump = FALSE,
         sep = "\t",
         phenotype_col = "TNFB",
         snp_col = "rsid",
         beta_col = "beta",
         se_col = "se",
         eaf_col = "EAF",
         effect_allele_col = "A1",
         other_allele_col = "A2",
         pval_col = "P",
         samplesize_col = "N",
         id_col = "rsid",
         log_pval = TRUE)
    pdf(file.path(INF,"MS","MS-forestplot.pdf"))
    for(id in c("ieu-b-18","ukb-b-17670","finn-a-G6_MS"))
    {
      cat("--",pQTL,"-",id,"--\n")
      y <- extract_outcome_data(with(x,SNP), id, proxies = TRUE, rsq = 0.8)
      xy <- mr(harmonise_data(x, y))
      forest_plot(xy)
      print(xy)
    }
    dev.off()
  END
}

wgs > ${INF}/MS/${prot}-MS-MR.log

function pqtl_qtl_rsid()
{
  export start=$(expr ${pos} - ${M})
  export end=$(expr ${pos} + ${M})
  if [ ${get_data} == "yes" ]; then
    (
      awk -vprot=${prot} -vOFS="\t" 'BEGIN{print prot, "rsid", "chr", "pos", "beta", "se", "snpid", "A1", "A2", "EAF", "P", "N"}'
      gunzip -c ${INF}/METAL/${prot}-1.tbl.gz | \
      awk -vFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} '(NR>1 && $1 == chr && $2 >= start && $2 <= end) {
          split($3,a,"_"); print a[1],$1,$2,$10,$11,$3,toupper($4),toupper($5),$6,-$12,$18
      }' | \
      sort -k1,1 | \
      join -12 -21 ${INF}/work/snp_pos - | \
      awk 'a[$7]++==0' | \
      awk -vprot=${prot} -vOFS="\t" '{print prot, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' | \
      sort -k3,3n -k4,4n
    ) > ${INF}/MS/${prot}-pQTL-${rsid}.dat
  fi
  R --no-save -q <<\ \ END
    options(echo=FALSE,width=200)
    INF <- Sys.getenv("INF")
    prot <- Sys.getenv("prot")
    rsid <- Sys.getenv("rsid")
    pQTL <- file.path(INF,"MS",paste0(prot,"-pQTL-",rsid,".dat"))
    library(TwoSampleMR)
    x <- subset(read_exposure_data(pQTL,
                clump = TRUE,
                sep = "\t",
                phenotype_col = "TNFB",
                snp_col = "rsid",
                beta_col = "beta",
                se_col = "se",
                eaf_col = "EAF",
                effect_allele_col = "A1",
                other_allele_col = "A2",
                pval_col = "P",
                samplesize_col = "N",
                id_col = "rsid",
                log_pval = TRUE), pval<1e-5)
    for(id in c("ieu-b-18","ukb-b-17670","finn-a-G6_MS"))
    {
      cat("##",pQTL,"-",id,"##\n")
      y <- extract_outcome_data(with(x,SNP), id, proxies = TRUE, rsq = 0.8)
      xy <- mr(harmonise_data(x, y))
      print(xy)
    }
  END
}

function pqtl()
# pQTLs
{
  (
     awk -vOFS="\t" 'BEGIN{print "Phenotype", "SNP", "chr", "pos", "beta", "se", "snpid", "effect_allele", "other_allele", "eaf", "pval", "N"}'
     zgrep -e chr12:6514963_A_C -e chr12:6440009_C_T -e chr6:31540757_A_C -e chr6:31073047_A_G ${INF}/METAL/${prot}-1.tbl.gz | \
     awk -vFS="\t" '{split($3,a,"_"); print a[1],$1,$2,$10,$11,$3,toupper($4),toupper($5),$6,10^$12,$18}' | \
     sort -k1,1 | \
     join -12 -21 ${INF}/work/snp_pos - | \
     awk 'a[$7]++==0' | \
     awk -vprot=${prot} -vOFS="\t" '{print prot, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' | \
     sort -k3,3n -k4,4n
  ) > ${INF}/MS/${prot}-pQTL.dat
  R --no-save -q <<\ \ END
    library(pQTLtools)
    INF <- Sys.getenv("INF")
    prot <- Sys.getenv("prot")
    pqtl <- file.path(INF,"MS",paste0(prot,"-pQTL.dat"))
    ivs <- read.table(pqtl,as.is=TRUE,header=TRUE)
    setwd(file.path(INF,"MS"))
    for(id in c("ieu-b-18","ukb-b-17670","finn-a-G6_MS"))
    {
      pqtlMR(subset(ivs,SNP%in%"rs2364485"),id,prefix=paste0(prot,"-",id,"-rs2364485"))
      pqtlMR(subset(ivs,SNP%in%"rs1800693"),id,prefix=paste0(prot,"-",id,"-rs1800693"))
      pqtlMR(subset(ivs,SNP%in%"rs9263621"),id,prefix=paste0(prot,"-",id,"-rs9263621"))
      pqtlMR(subset(ivs,SNP%in%"rs2229092"),id,prefix=paste0(prot,"-",id,"-rs2229092"))
    }
  END
}

function pqtl_flanking()
# +/- 0.5Mb
{
  (
  # cis pQTLs
  # chr6:31540757_A_C rs2229092
  # chr6:31073047_A_G rs9263621
    export chr=6
    export rsid=rs2229092
    export pos=31540757
    pqtl_qtl_rsid
    export rsid=rs9263621
    export pos=31073047
    pqtl_qtl_rsid
  # trans pQTL
  # chr12:6514963_A_C rs2364485
  # chr12:111865049_C_G rs7310615; MS:chr12:6440009_C_T rs1800693
  # r2(rs1800693,rs2364485)=0.0029, r2(rs1800693,rs7310615)=0.0023, r2(rs2364485,rs7310615)=0.0013
    export chr=12
    export rsid=rs2364485
    export pos=6514963
    pqtl_qtl_rsid
    export rsid=rs1800693
    export pos=6440009
    pqtl_qtl_rsid
    export rsid=rs7310615
    export pos=111865049
    pqtl_qtl_rsid
  ) >> ${INF}/MS/${prot}-MS-MR.log
}

export M=60000
export get_data=no

pqtl_flanking
pqtl

R --no-save -q <<END
  info <- function()
  {
    # https://gwas.mrcieu.ac.uk/datasets/?trait__icontains=multiple%20sclerosis
      ms_ids <- c("ukb-b-17670","ieu-b-18","ieu-a-1025","ebi-a-GCST005531","ieu-a-1024","ebi-a-GCST001198","ieu-a-820","ieu-a-821","finn-a-G6_MS")
      ms_gwasinfo <- ieugwasr::gwasinfo(ms_ids)
      cols <- c("id","year","ncase","ncontrol","nsnp","build","author","pmid","population")
      ms_gwasinfo[cols]
#                 id year ncase ncontrol     nsnp       build         author     pmid population
# 1      ukb-b-17670 2018  1679   461254  9851867 HG19/GRCh37   Ben Elsworth       NA   European
# 2         ieu-b-18 2019 47429    68374  6304359 HG19/GRCh37 Patsopoulos NA 31604244   European
# 3       ieu-a-1025 2013 14498    24091   156632 HG19/GRCh37        Beecham 24076602   European
# 4 ebi-a-GCST005531 2013 14498    24091   132089 HG19/GRCh37     Beecham AH 24076602   European
# 5       ieu-a-1024 2011  9722    17376   465435 HG19/GRCh37       Sawcer S 21833088   European
# 6 ebi-a-GCST001198 2011  9772    16849   463040 HG19/GRCh37       Sawcer S 21833088   European
# 7        ieu-a-820 2007   931     1862   327095 HG19/GRCh37      Hafler DA 17660530   European
# 8        ieu-a-821 2009   978      883   514572 HG19/GRCh37   Baranzini SE 19010793   European
# 9     finn-a-G6_MS 2020   378    44677 16152119 HG19/GRCh37             NA       NA   European
# 2,9 with no VCF
# https://gwas.mrcieu.ac.uk/files/ieu-a-1025/ieu-a-1025.vcf.gz
      mr_info <- as.data.frame(epigraphdb::mr(outcome_trait = "Multiple sclerosis", pval_threshold = 1e-8))
      select <- c("ukb-a-100","ukb-a-104","ieu-a-294","ieu-a-295","ieu-a-971")
      subset(mr_info,exposure.id%in%select,select=-c(outcome.trait,mr.selection,mr.method,mr.moescore))
# exposure.id                                             exposure.trait outcome.id       mr.b      mr.se      mr.pval
#   ieu-a-971                                         Ulcerative colitis ieu-a-1025 -0.4839367 0.06111146 2.395840e-15
#   ieu-a-295                                 Inflammatory bowel disease ieu-a-1025 -0.2593652 0.03799687 8.733802e-12
#   ieu-a-294                                 Inflammatory bowel disease ieu-a-1025  0.1232999 0.01697710 1.375433e-10
#   ukb-a-104 Non-cancer illness code  self-reported: ulcerative colitis ieu-a-1025 20.4956226 3.23034883 2.228468e-10
#   ukb-a-100          Non-cancer illness code  self-reported: psoriasis ieu-a-1025 11.0119963 1.86965907 3.865659e-09
      mr_summary <- epigraphdb::mr(outcome_trait="Multiple sclerosis")
      names(mr_summary)[c(2,3,5,6,7)] <- c("exposure","outcome","b","se","pval")
      tryx::volcano_plot(mr_summary)
  }
END

# https://www.gtexportal.org/home/gene/RP1-102E24.8
# no sample size (N)!!!
export ichip=~/rds/results/public/gwas/multiple_sclerosis/ImmunoChip_Results/Immunochip_FinalResults_LimitedDiscovery.txt
awk -vchr=${chr} -vstart=${start} -vend=${end} '$1==chr&&$3>=start&&$3<=end {print $4,$5,$6,$7,log($10),$11,$9}' ${ichip}
# CHR BPHG18 BPHG19 ImmunochipID Risk_Allele Ref_Allele Risk_Allele_Freq N P OR SE Q I Region
# 1 1108138 1118275 vh_1_1108138 G A 0.9586 11 7.62E-01 1.012 0.0404 0.085 39.51 none

# --- legacy ---

function ieu_id_ma()
  (
    echo SNP A1 A2 freq b se p N
    bcftools query -f "%CHROM %POS %ALT %REF %AF [%ES] [%SE] [%LP] [%SS]\n" -r ${region} \
                   ~/rds/results/public/gwas/multiple_sclerosis/${ieu_id}.vcf.gz | \
    awk '{if ($3<$4) snpid="chr"$1":"$2"_"$3"_"$4;else snpid="chr"$1":"$2"_"$4"_"$3;print snpid, $3, $4, $5, $6, $7, $8, $9}' | \
    awk 'a[$1]++==0 && $8<5 {$7=10^-$7};1'
  ) > ${INF}/MS/EUR-${ieu_id}.ma

function MS2013()
{
export chr=12
export start=6400000
export end=6520000
export M=0
export gene=LTBR
export prot=TNFB

echo Multiple sclerosis
# rs1800693 chr12:6440009
# rs2364485 chr12:6514963
# 2018, the rsid is incomplete
zgrep -w -e 6440009 -e 6514963 data/discovery_metav3.0.meta.gz
# 12 6440009 rs1800693 T C 14 1.017e-13 0.8808
# 12 6514963 rs2364485 C A 15 5.778e-06 0.9041
gunzip -c data/discovery_metav3.0.meta.gz | \
sed '1d' | \
awk '!/^$/' | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if ($1==chr && $2>=start-M && $2<=end+M)
  {
    if ($4<$5) snpid="chr" $1 ":" $2 "_" $4 "_" $5;
    else snpid="chr" $1 ":" $2 "_" $5 "_" $4
    print "chr"$1":"$2,$1,$2,$7,snpid,$4,$5
  }
}' | \
sort -k1,1 | \
join -12 -21 work/snp_pos - | \
awk -vOFS="\t" '{print $2,$3,$4,$5,$6,$7,$8}' > work/${prot}-QTL.lz

echo 2013
zgrep -w -e rs1800693 -e rs2364485 data/imsgc_2013_24076602_ms_efo0003885_1_ichip.sumstats.tsv.gz
zgrep -v -w -e NaN -e Infinity -e 0.0 data/imsgc_2013_24076602_ms_efo0003885_1_ichip.sumstats.tsv.gz | \
sed '1d' | \
awk '!/^$/' | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if ($1==chr && $2>=start-M && $2<=end+M)
  {
    if ($4<$5) snpid="chr" $1 ":" $2 "_" $4 "_" $5;
    else snpid="chr" $1 ":" $2 "_" $5 "_" $4
    print $3,$1,$2,$7/$8,snpid,$5,$4
  }
}' > work/${prot}-2013.lz

echo SCALLOP/INF
# chr12:6514963_A_C
gunzip -c METAL/${prot}-1.tbl.gz | \
sed '1d' | \
awk -vFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} '
{
  if ($1 == chr && $2 >= start-M && $2 <= end+M)
  {
    split($3,a,"_")
    print a[1],$1,$2,$10/$11,$3,toupper($4),toupper($5)
  }
}' | \
sort -k1,1 | \
join -12 -21 work/snp_pos - | \
awk 'a[$6]++==0' | \
awk -vOFS="\t" '{print $2, $3, $4, $5, $6, $7, $8}' > work/${prot}-pQTL.lz

cut -f5 work/${prot}-pQTL.lz > work/${prot}-pQTL.snpid
plink --bfile INTERVAL/cardio/INTERVAL --extract work/${prot}-pQTL.snpid \
      --r2 inter-chr yes-really --ld-snps chr12:6514963_A_C --ld-window-r2 0 --out work/${prot}-pQTL
(
  awk -vOFS="\t" 'BEGIN{print "snpid","rsid","chr","pos","z","A1","A2","r2"}'
  join -15 -t$'\t' work/${prot}-pQTL.lz <(awk -vOFS="\t" 'NR>1 {print $6,$7}' work/${prot}-pQTL.ld | sort -k1,1)
) > work/${prot}-pQTL.txt

echo eQTLGen
zgrep -w ${gene} data/eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if ($5<$6) snpid="chr" $3 ":" $4 "_" $5 "_" $6;
  else snpid="chr" $3 ":" $4 "_" $6 "_" $5
  if($3==chr && $4>=start-M && $4 <=end+M) print $2,$3,$4,$7,snpid,$5,$6
}' > work/${prot}-eQTL.lz

join -a2 -e "NA" -o2.5,2.1,2.2,2.3,1.4,1.6,1.7,2.1,2.2,2.3,2.4,2.6,2.7 \
     -j5 <(sort -k5,5 work/${prot}-QTL.lz) <(sort -k5,5 work/${prot}-pQTL.lz) | \
join -a1 -e "NA" -25 -o1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.1,2.2,2.3,2.4,2.6,2.7 \
     - <(sort -k5,5 work/${prot}-eQTL.lz) | \
awk -vOFS="\t" '
{
# 2013 with beta/se
# if($6!="NA" && $12!="NA" && $11!="NA" && $6!=12) $11=-$11
# if($6!="NA" && $18!="NA" && $17!="NA" && $6!=18) $17=-$17
  print $1,$2,$3,$4,$5,$11,$17
}' | \
awk 'a[$1]++==0' | \
awk 'a[$2]++==0' | \
sort -k3,3n -k4,4n > work/${prot}.z

cut -f1 work/${prot}.z > work/${prot}.snpid
plink --bfile INTERVAL/cardio/INTERVAL --extract work/${prot}.snpid --make-bed --out work/${prot}
cut -f2 work/${prot}.bim > work/${prot}.snpid
plink --bfile work/${prot} --extract work/${prot}.snpid --r square --out work/${prot}
grep -f work/${prot}.snpid work/${prot}.z > work/${prot}.gassoc

R --no-save -q <<END
  prot <- Sys.getenv("prot")
  d <- read.table(paste0(file.path("work",prot),".gassoc"),
                         col.names=c("snpid","marker","chr","pos","Multiple_sclerosis","pQTL","eQTL"))
  d <- within(d,{Multiple_sclerosis <- qnorm(Multiple_sclerosis/2)})
  markers <- d[c("marker","chr","pos")]
  z <- d[c("Multiple_sclerosis","pQTL","eQTL")]
  rownames(z) <- with(d,marker)
  ld <- read.table(paste0(file.path("work",prot),".ld"),col.names=with(d,marker),row.names=with(d,marker))
  library(gassocplot)
  for (rsid in c("rs1800693","rs2364485"))
  {
    sap <- stack_assoc_plot(markers, z, ld, top.marker=rsid, traits = c("Multiple_sclerosis","pQTL","eQTL"), ylab = "-log10(P)", legend=TRUE)
    stack_assoc_plot_save(sap, paste0(file.path("work",prot),"-",rsid,".png"), 3, width=8, height=13, dpi=300)
  }
  z <- within(z,{Multiple_sclerosis <- NA})
# rs1800693: P 2e-47 Z -14.46555
  pos1 <- with(d,marker=="rs1800693")
  z[pos1,"Multiple_sclerosis"] <- -14.46555
# rs2364485: P 2e-20 Z  -9.26234
  pos2 <- with(d,marker=="rs2364485")
  z[pos2,"Multiple_sclerosis"] <- -9.26234
  rsid <- "rs2364485"
  sap <- stack_assoc_plot(markers, z, ld, top.marker=rsid, traits = c("Multiple_sclerosis","pQTL","eQTL"), ylab = "-log10(P)", legend=TRUE)
  stack_assoc_plot_save(sap, paste0(file.path("work",prot),"-",rsid,"-fixed.png"), 3, width=8, height=13, dpi=300)
END

# Multiple_sclerosis

# 2018
# CHR BP SNP A1 A2 N P OR
# 1 11154 chr1:11154 C A 4 0.7911 0.9818

# 2013
# chrom	pos	rsid	other_allele	effect_allele	p	beta	se	OR	OR_lower	OR_upper
# 12      6440009 rs1800693       A       G       6.92e-16        0.13453089295760606     0.01666652509712173     1.144   1.1072334356041962      1.1819874273267836
}
