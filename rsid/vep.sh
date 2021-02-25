#!/usr/bin/bash

function prune()
{
  export sq="'"
  sed 's/\t/ /g' ${INF}/work/INF1.merge.alleles | parallel -C' ' --env INF '
    echo {1}-{2}
    export rt=${INF}/annotate/{1}-{2}
    plink --bfile ${INF}/INTERVAL/cardio/INTERVAL --snp {2} --window 1000 \
          --geno 0.1 --mind 0.1 --maf 0.005 --indep-pairwise 1000kb 1 0.8 --out ${rt}
    export incl=$(grep -w {2} ${rt}.prune.in)
    if [ "${incl}" == "" ]; then
       export i=$(grep -w -f ${rt}.prune.in ${INF}/INTERVAL/cardio/INTERVAL.bim | \
                  awk -vpos={4} "function abs(x) {if (x<0) return -x; else return x;} {d=abs(\$4-pos);print \$1, \$2, \$4, d}" | \
                  sort -r -k4,4n | \
                  awk "NR==1 {print \$2}" \
                )
       sed -i "s/${i}/{2}/;s/${sq}//g" ${rt}.prune.in
    fi
    (
      echo "##fileformat=VCFv4.0"
      echo "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO" | \
      tr "," "\t"
      sed "s/chr//" ${rt}.prune.in | \
      awk -vOFS="\t" "{snpid=\$1;gsub(/:|_/,\" \",snpid);split(snpid,a,\" \");print a[1],a[2],\"chr\"\$1,a[3],a[4],\".\",\".\",\".\"}"
      sort -k1,1n -k2,2n
    ) > ${rt}.vcf
    rm ${rt}.prune.* ${rt}.log ${rt}.nosex
  '
}

cd ${HPC_WORK}/loftee
for f in $(ls ${INF}/annotate/*vcf | sed 's/.vcf//g')
do
  echo ${f}
# snpid --> rsid
  cp ${f}.vcf ${f}.vepinput
  (
    sed '1,2d' ${f}.vcf | cut -f3 | sort | join - ${INF}/work/INTERVAL.rsid | \
    parallel --dry-run -C' ' "
      export s={1};
      export r={2};
      sed -i 's/'\"\${s}\"'/'\"\${r}\"'/g' ${f}.vepinput
    "
  ) | bash
# VEP annotation
  vep --input_file ${f}.vepinput --output_file ${f}.tab --force_overwrite \
      --cache --dir_cache ${HPC_WORK}/ensembl-vep/.vep --dir_plugins ${HPC_WORK}/loftee --offline \
      --species homo_sapiens --assembly GRCh37 --pick --nearest symbol --symbol --plugin TSSDistance \
      --plugin LoF,loftee_path:.,human_ancestor_fa:human_ancestor.fa.gz,conservation_file:phylocsf_gerp.sql.gz \
      --tab
# per_gene
  vep -i ${f}.vcf -o ${f}.pergene --per_gene --check_existing --force_overwrite --offline \
      --species homo_sapiens --everything --assembly GRCh37 --nearest symbol \
      --symbol --pubmed --uniprot --protein --sift b --polyphen b --tab
done
cd -
# Problematic with loftee
# CCL11-chr7:75495667_A_G
# CCL19-chr6:32444093_C_G
# CCL25-chr12:578100_A_G
# CCL25-chr19:49206674_A_G
# CDCP1-chr6:32602396_C_T
# CST5-chr12:11058117_C_T
# CST5-chr19:49206145_C_G
# CX3CL1-chr6:32424882_C_T
# CXCL1-chr4:74739076_G_T
# CXCL5-chr4:74858051_A_G
# CXCL6-chr4:74703999_C_T
# FGF.19-chr19:49206172_C_T
# FGF.21-chr19:49260677_A_C
# IL.10-chr1:206954566_A_G
# IL.10-chr6:32434716_A_C
# IL.12B-chr6:31154493_A_G
# IL.1.alpha-chr6:32586222_A_G
# IL.6-chr1:154426970_A_C
# IL.8-chr4:74574265_A_G
# LAP.TGF.beta.1-chr19:41847860_A_G
# MMP.10-chr11:102649482_C_T
# MMP.10-chr19:49206145_C_G
# MMP.1-chr11:102697731_A_G
# SCF-chr19:54793830_C_G
# SCF-chr9:107661742_A_C
# SLAMF1-chr17:7106378_A_G
# SLAMF1-chr5:95263427_A_G
# ST1A1-chr16:28561581_C_T
# TNFB-chr6:31540757_A_C
# TRAIL-chr1:196710916_C_T
# TRANCE-chr8:23085868_A_G
# uPA-chr17:7063667_C_T
