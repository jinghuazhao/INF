#!/usr/bin/bash

function prune()
{
  sed 's/\t/ /g' ${INF}/work/INF1.merge.alleles | parallel -C' ' --env INF '
    echo {1}-{2}
    export rt=${INF}/annotate/{1}-{2}
    plink --bfile ${INF}/INTERVAL/cardio/INTERVAL --snp {2} --window 1000 \
          --geno 0.1 --mind 0.1 --maf 0.005 --indep-pairwise 1000kb 1 0.8 --out ${rt}
    if [ $(grep -w {2} ${rt}.prune.in | wc -l) -eq 0 ]; then
       export i=$(grep -w -f ${rt}.prune.in ${INF}/INTERVAL/cardio/INTERVAL.bim | \
                  awk -vpos={4} "function abs(x) {if (x<0) return -x; else return x;} {d=abs(\$4-pos);print \$1, \$2, \$4, d}" | \
                  sort -r -k4,4n | \
                  awk "NR==1 {print \$2}" \
                )
       sed -i 's/'"$i"'/'"{2}"'/g' ${rt}.prune.in
    fi
    (
      echo "##fileformat=VCFv4.0"
      echo "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO" | \
      tr "," "\t"
      sed "s/chr//" ${rt}.prune.in | \
      awk -vOFS="\t" "{snpid=\$1;gsub(/:|_/,\" \",snpid);split(snpid,a,\" \");print a[1],a[2],\$1,a[3],a[4],\".\",\".\",\".\"}"
      sort -k1,1n -k2,2n
    ) > ${rt}.vcf
    rm ${rt}.prune.* ${rt}.log ${rt}.nosex
  '
}

export LOFTEE=${HPC_WORK}/loftee
export VEP=${HPC_WORK}/ensembl-vep
export TMPDIR=/rds/user/jhz22/hpc-work/work

cd ${LOFTEE}
for f in $(ls ${INF}/annotate/*vcf | sed 's/.vcf//g')
do
# snpid --> rsid
  cp ${f}.vcf ${f}-rsid.vcf
  (
    sed '1,2d' ${f}.vcf | cut -f3 | sort | join - ${INF}/work/INTERVAL.rsid | \
    parallel --dry-run -C' ' "
      export s={1};
      export r={2};
      sed -i 's/'\"\${s}\"'/'\"\${r}\"'/g' ${f}-rsid.vcf
    "
  ) | bash
# VEP annotation
  vep --input_file ${INF}/work/${f}-rsid.vepinput --output_file ${INF}/annotate/${f}-rsid.tab --force_overwrite \
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
