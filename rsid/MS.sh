#!/usr/bin/bash

export prot=TNFB
export chr=12
export start=6300000
export end=6700000
export region=${chr}:${start}-${end}
export dir=~/rds/results/public/gwas/blood_cell_traits/chen_2020
export TMPDIR=/rds/user/jhz22/hpc-work/work

function slct()
{
  module load plink/2.00-alpha
  for cell in wbc mono neut lymph eo baso ieu-a-1025
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
      if [ ${cell} != "ieu-a-1025" ]; then export P_threshold=5e-8; else export P_threshold=1e-5; fi
      gcta-1.9 --bfile ${INF}/INTERVAL/cardio/INTERVAL \
               --cojo-file ${INF}/MS/EUR-${cell}.ma \
               --extract ${INF}/MS/EUR-${cell}.prune \
               --cojo-slct \
               --cojo-p ${P_threshold} \
               --maf 0.005 \
               --cojo-collinear 0.9 \
               --out ${INF}/MS/EUR-${cell}
  done
}

slct

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
  export ieu_id=ieu-a-1025
  (
    echo SNP A1 A2 freq b se p N
    bcftools query -f "%CHROM %POS %ALT %REF %AF [%ES] [%SE] [%LP] [%SS]\n" -r ${region} \
                   ~/rds/results/public/gwas/multiple_sclerosis/${ieu_id}.vcf.gz | \
    awk '{if ($3<$4) snpid="chr"$1":"$2"_"$3"_"$4;else snpid="chr"$1":"$2"_"$4"_"$3;print snpid, $3, $4, $5, $6, $7, $8, $9}' | \
    awk 'a[$1]++==0 && $8<5 {$7=10^-$7};1'
  ) > ${INF}/MS/EUR-${ieu_id}.ma
  sed "1d" ${INF}/MS/EUR-${ieu_id}.ma | \
  cut -d" " -f1 > ${INF}/MS/EUR-${ieu_id}.snpid
  sed "1d" ${INF}/MS/EUR-${ieu_id}.ma | \
  sort -k7,7g | \
  awk "NR==1{print \$1}" > ${INF}/MS/EUR-${ieu_id}.top
}

ma

function lz()
{
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
  rm -f ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${INF}/MS/{}.lz --delim tab title="{}" \
            --markercol rsid --pvalcol p_value --chr ${chr} --start ${start} --end ${end} \
            --no-date --plotonly --prefix="{}" --rundir .
  qpdf {}_chr${chr}_${start}-${end}.pdf --pages . 1 -- {}-lz.pdf
  '
  qpdf --empty --pages *lz.pdf -- blood-cell-traits.pdf
}

lz

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

function gsmr()
{
# MS 
  export ieu_id=ieu-a-1025
  (
    echo SNP A1 A2 freq b se p N
    bcftools query -f "%CHROM %POS %ALT %REF %AF [%ES] [%SE] [%LP] [%SS]\n" ~/rds/results/public/gwas/multiple_sclerosis/${ieu_id}.vcf.gz | \
    awk '{if ($3<$4) snpid="chr"$1":"$2"_"$3"_"$4;else snpid="chr"$1":"$2"_"$4"_"$3;print snpid, $3, $4, $5, $6, $7, $8, $9}' | \
    awk 'a[$1]++==0'
  ) | gzip -f > ${INF}/MS/gsmr_MS-${ieu_id}.gz
  R --no-save -q <<\ \ END
    INF <- Sys.getenv("INF")
    ieu_id <- Sys.getenv("ieu_id")
    ma <- file.path(INF,"MS",paste0("gsmr_MS-",ieu_id,".gz"))
    ms <- within(read.table(ma,header=TRUE),{p <- 10^-p})
    f <- gzfile(file.path(INF,"MS","gsmr_MS.ma.gz"))
    write.table(ms,file=f,quote=FALSE,row.names=FALSE)
    unlink(ma)
  END
# prot
  (
    echo SNP A1 A2 freq b se p N
    zcat ${INF}/METAL/${prot}-1.tbl.gz | \
    awk 'NR>1 {print $3,toupper($4),toupper($5),$6,$10,$11,10^$12,$18}'
  ) | gzip -f > ${INF}/MS/gsmr_${prot}.ma.gz
# control files
  if [ ! -f ${INF}/MS/gsmr_ref_data ]; then echo ${INF}/work/INTERVAL > ${INF}/MS/gsmr_ref_data; fi
  if [ ! -f ${INF}/MS/gsmr_MS ]; then echo MS ${INF}/MS/gsmr_MS.ma.gz > ${INF}/MS/gsmr_MS; fi
  if [ ! -f ${INF}/MS/gsmr_${prot} ]; then echo ${prot} ${INF}/MS/gsmr_${prot}.ma.gz > ${INF}/MS/gsmr_${prot}; fi

  gcta-1.9 --mbfile ${INF}/MS/gsmr_ref_data --gsmr-file ${INF}/MS/gsmr_MS ${INF}/MS/gsmr_${prot} \
           --gsmr-direction 0 \
           --clump-r2 0.05 --gwas-thresh 1e-5 --diff-freq 0.4 --heidi-thresh 0.05 --gsmr-snp-min 5 --effect-plot \
           --out ${INF}/MS/gsmr_${prot}-MS

  R --no-save -q <<\ \ END
    prot <- Sys.getenv("prot")
    source("http://cnsgenomics.com/software/gcta/res/gsmr_plot.r")
    gsmr_data <- read_gsmr_data(paste0("MS/gsmr_",prot,"-MS.eff_plot.gz"))
    gsmr_summary(gsmr_data)
    pdf(paste0("MS/gsmr_",prot,"-MS.eff_plot.pdf"))
    plot_gsmr_effect(gsmr_data, prot, "MS", colors()[75])
    dev.off()
  END
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
    for(id in c("ieu-b-18","ieu-a-1025","ukb-b-17670","finn-a-G6_MS"))
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
    for(id in c("ieu-b-18","ieu-a-1025","ukb-b-17670","finn-a-G6_MS"))
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
    for(id in c("ieu-b-18","ieu-a-1025","ukb-b-17670","finn-a-G6_MS"))
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
  # r2=0.0029
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
