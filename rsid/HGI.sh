#!/usr/bin/bash

if [ ! -d ${INF}/HGI/mr ]; then mkdir ${INF}/HGI/mr; fi

# gs://covid19-hg-public/20210415/results/20210607/
# gs://covid19-hg-public/20210415/results/ancestry/20210622

function _outcome_r6()
{
  export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20210415/results/20210607
  for trait in A2 B1 B2 C2
  do
    export trait=${trait}
    export src=${HGI}/COVID19_HGI_${trait}_ALL_leave_23andme_20210607.b37.txt.gz
    cut -f5,6,8,9 --output-delimiter=' ' ${INF}/work/INF1.merge-rsid | \
    sed '1d' | \
    awk -vM=1e6 '{start=$4-M;if(start<0) start=1;end=$4+M;print $1,$2,$3":"start"-"end}' | \
    parallel -j15 --env src -C' ' 'tabix -h ${src} {3} | gzip -f > ${INF}/HGI/mr/r6-${trait}-{1}-{2}.gz'
  done
}
_outcome_r6

function _ins()
(
  echo SNP Phenotype Allele1 Allele2 EAF Effect StdErr P N
  sed '1d' ${INF}/work/INF1.METAL | \
  cut -f2,3,6,7,8-11,17,21 | \
  awk '{print $1,$2,toupper($3),toupper($4),$5,$6,$7,10^$8,$9}'
) > ${INF}/HGI/mr/INF.ins
_ins

function _exposure_tsv()
{
  for trait in A2 B1 B2 C2
  do
    export trait=${trait}
    cut -f5,6,8,9 --output-delimiter=' ' ${INF}/work/INF1.merge-rsid | \
    sed '1d' | \
    awk -vM=1e6 '{start=$4-M;if(start<0) start=1;end=$4+M;print $1,$2,$3":"start"-"end}' | \
    parallel -j15 --env INF --env trait -C' ' '
    (
      echo -e "prot\trsid\tChromosome\tPosition\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tP\tN"
      tabix ${INF}/METAL/gwas2vcf/{1}.tsv.gz {3} | \
      awk "{\$4=toupper(\$4);\$5=toupper(\$5);print}" | \
      sort -k3,3 | \
      join -23 ${INF}/work/INTERVAL.rsid - | \
      awk -v prot={1} "a[\$2]++==0{\$1=prot;print}" | \
      tr " " "\t"
    ) | gzip -f > ${INF}/HGI/mr/tsv/r6-${trait}-{1}-{2}.tsv.gz
    '
  done
}
_exposure_tsv

module load gcc/6

for trait in A2 B1 B2 C2
do
  export trait=${trait}
  sed '1d' ${INF}/HGI/mr/INF.ins | \
  parallel -j15 -C' ' '
    export prot={2}; export rsid={1};
    export Allele1={3}; export Allele2={4}; export EAF={5}
    export Effect={6}; export StdErr={7}; export P={8}; export N={9}
    R --no-save < ${INF}/rsid/HGI.R
  '
done

(
  awk -vOFS="\t" 'BEGIN{print "Method","Batch","Trait","Protein","pQTL","b","se","p"}'
  awk -vOFS="\t" '/Inverse/ {print FILENAME,$9,$10,$11}' ${INF}/HGI/mr/MR*r6*result*
  awk -vOFS="\t" '!/pval/{print FILENAME,$8,$9,$10}' ${INF}/HGI/mr/pqtl*r6*result*
) | \
sed 's|/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/HGI/mr/||;s|-result.txt||' | \
awk -vOFS="\t" '{gsub(/-/,"\t",$1)};1' > ${INF}/HGI/mr.tsv

R --no-save -q <<END
  options(width=200)
  library(dplyr)
  INF <- Sys.getenv("INF")
  METAL <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
           select(prot,rsid,cis.trans,Chromosome,Position) %>%
           rename(pQTL=rsid)
  tsv <- read.delim(file.path(INF,"HGI","mr.tsv")) %>%
         filter(Method=="pqtlMR") %>%
         mutate(fdr=p.adjust(p,method="fdr")) %>%
         left_join(METAL)
  txt <- left_join(tsv,gap.datasets::inf1[c("prot","target.short")],by=c("Protein"="prot")) %>%
         mutate(Protein=target.short) %>%
         select(-c(Method,Batch,prot,target.short)) %>%
         arrange(desc(fdr))
  subset(txt,fdr<=0.05)
  write.table(txt,file=file.path(INF,"HGI","pqtlMR.txt"),quote=FALSE,row.names=FALSE,sep="\t")
END

R --no-save -q <<END
  regions <- data.frame(chr="chr19",start=49206145,end=49206674)
  job <- submitGreatJob(regions, species="hg19", version="3.0.0")
  et <- getEnrichmentTables(job,download_by = 'tsv')
  tb <- do.call('rbind',et)
  write.table(tb,file=paste0("tb.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
END

# --- optional replacement ---

function etc()
{
# grep -e chr3 -e chr9 work/INF1.merge-rsid | cut -f6 | zgrep -f - -w $C2
  ls *-rs* | sed 's/-/ /g' | \
  parallel -C' ' '
      echo {1} - {2};
      grep -w {2} {1}-{2};
      grep -w -e {1} work/INF1.merge.cis.vs.trans-rsid | grep -w {2};
      grep -w -e {1} INF.ins | grep -e {2}
  '
  ls ${INF}/HGI/mr/no-clumping/MR*result* | \
  sed 's/-/ /g;' | \
  cut -d' ' -f5-7 | \
  parallel -C' ' 'echo {1} {2} {3}'
}

function _outcome_r4()
{
  export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20200915/results/20201020
  grep -e chr3 -e chr9 work/INF1.merge-rsid | \
  cut -f5,6,8,9 --output-delimiter=' ' | \
  parallel -C' ' '
     gunzip -c ${HGI}/COVID19_HGI_C2_ALL_leave_23andme_20201020.b37.txt.gz | \
     awk -v M=1e6 -v chr={3} -v pos={4} "NR==1||(\$1==chr && \$2>=pos-M && \$2<=pos+M)" > work/HGI/{1}-{2}
     cat ${HGI}/COVID19_HGI_C2_ALL_20201020.b37_1.0E-5.txt | \
     awk -v M=1e6 -v chr={3} -v pos={4} "NR==1||(\$1==chr && \$2>=pos-M && \$2<=pos+M)" > work/HGI/{1}-{2}_1e-5
  '
}

function _outcome_r5()
{
  export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20201215/results/20210107
  for trait in A2 B2 C2
  do
    export trait=${trait}
    cut -f5,6,8,9 --output-delimiter=' ' ${INF}/work/INF1.merge-rsid | \
    sed '1d' | \
    parallel --env trait -C' ' '
       gunzip -c ${HGI}/COVID19_HGI_${trait}_ALL_eur_leave_23andme_20210107.b37.txt.gz | \
       awk -v M=1e6 -v chr={3} -v pos={4} "NR==1||(\$1==chr && \$2>=pos-M && \$2<=pos+M)" > ${INF}/HGI/mr/${trait}-{1}-{2}
    '
  done
}

function _outcome_r6()
{
  export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20210415/results/20210614/
  export TMPDIR=${HPC_WORK}/work
  export b38tob37=~/hpc-work/bin/hg38ToHg19.over.chain.gz
  for trait in A2 B1 B2 C2
  do
    export trait=${trait}
    export src=${HGI}/COVID19_HGI_${trait}_ALL_leave_23andme_20210607.txt.gz
    if [ ${trait} == "A2" ]; then export fields=82-91
    elif [ ${trait} == "B1" ]; then export fields=76-85
    elif [ ${trait} == "B2" ]; then export fields=136-145
    elif [ ${trait} == "C2" ]; then export fields=229-238
    fi
    (
      echo -e "#CHROM\tstart\tend\tSNP"
      gunzip -c ${src} | \
      sed '1d' | \
      cut -f1,2,5 | \
      awk -vFS="\t" -vOFS="\t" '{$1="chr" $1 "\t" $2-1};1'
    ) > ${INF}/HGI/r6-b38-${trait}.bed
    liftOver ${INF}/HGI/r6-b38-${trait}.bed ${b38tob37} ${INF}/HGI/r6-b37-${trait}.bed ${INF}/HGI/r6-b37-${trait}.unlifted.bed
    join -13 <(cut -f1,3,4 ${INF}/HGI/r6-b37-${trait}.bed | sed 's/chr//' | sort -k3,3) \
             <(gunzip -c ${src} | sed '1d' | cut -f5,${fields} | sort -k1,1) \
    > ${INF}/HGI/r6-b37-${trait}.dat
    cut -f5,6,8,9 --output-delimiter=' ' ${INF}/work/INF1.merge-rsid | \
    sed '1d' | \
    parallel --env trait -C' ' '
    (
       awk -v snp="SNP" \
           -v chr="chr" \
           -v pos="pos" \
           -v c1="all_inv_var_meta_beta" \
           -v c2="all_inv_var_meta_sebeta" \
           -v c3="all_inv_var_meta_p" \
           -v c4="all_inv_var_meta_cases" \
           -v c5="all_inv_var_meta_controls" \
           -v c6="all_inv_var_meta_effective" \
           -v c7="all_inv_var_het_p" \
           -v c8="all_meta_sample_N" \
           -v c9="all_meta_AF" \
           -v c10="rsid" \
           "BEGIN{print snp,chr,pos,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10}"
       awk -v M=1e6 -v chr={3} -v pos={4} "\$2==chr && \$3>=pos-M && \$3<=pos+M" ${INF}/HGI/r6-b37-${trait}.dat | \
       sort -k2,2n -k3,3n
    ) > ${INF}/HGI/mr/r6-${trait}-{1}-{2}
    '
  done
}
_outcome_r6

function _exposure_vcf()
  for trait in A2 B2 C2
  do
    export trait=${trait}
    cut -f5,6,8,9 --output-delimiter=' ' ${INF}/work/INF1.merge-rsid | \
    sed '1d' | \
    awk -vM=1e6 '{start=$4-M;if(start<0) start=1;end=$4+M;print $1,$2,$3":"start"-"end}' | \
    parallel -j15 --env INF --env trait -C' ' '
    (
      echo -e "prot\trsid\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tlogP\tN\tChromosome\tPosition"
      bcftools query -f "{1}\t%ID\t%ALT\t%REF\t%AF\t[%ES]\t[%SE]\t[%LP]\t[%SS]\t%CHROM\t%POS\n" \
                     -r {3} ${INF}/METAL/gwas2vcf/{1}.vcf.gz | \
      awk "a[\$2]++==0"
    ) | \
    gzip -f > ${INF}/HGI/mr/${trait}-{1}-{2}.gz
    '
  done
_exposure_vcf

function _exposure_METAL()
{
  for trait in A2 B2 C2
  do
    export trait=${trait}
    cut -f5,6,8,9 --output-delimiter=' ' ${INF}/work/INF1.merge-rsid | \
    sed '1d' | \
    parallel -j15 --env trait -C' ' '
      gunzip -c ${INF}/METAL/{1}-1.tbl.gz | \
      cut -f1-6,10-12,18 | \
      awk -vchr={3} -vpos={4} -vM=1e6 "\$1==chr && \$2>=pos-M && \$2 <= pos+M" > ${INF}/HGI/mr/${trait}-{1}-{2}.mri
      (
        echo -e "prot\trsid\tChromosome\tPosition\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tlogP\tN"
        awk "{\$4=toupper(\$4);\$5=toupper(\$5);print}" ${INF}/HGI/mr/${trait}-{1}-{2}.mri | \
        sort -k3,3 | \
        join -23 ${INF}/work/INTERVAL.rsid - | \
        awk -v prot={1} "a[\$2]++==0{\$1=prot;print}" | \
        tr " " "\t"
      ) | gzip -f > ${INF}/HGI/mr/${trait}-{1}-{2}.mrx
    '
  done
}
_exposure_METAL
