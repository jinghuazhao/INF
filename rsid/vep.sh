#!/usr/bin/bash

export LOFTEE=${HPC_WORK}/loftee
export VEP=${HPC_WORK}/ensembl-vep
export TMPDIR=/rds/user/jhz22/hpc-work/work
export dir=${INF}

for f in INF1.proxy INF1.proxy.trans
do
# snpid --> rsid
  cd ${INF}/work
  cp ${f}.vepinput ${f}-rsid.vepinput
  (
  sed '1,2d' ${f}.vepinput | cut -f3 | sort | join - INTERVAL.rsid | \
  parallel --dry-run -C' ' "
    export s={1};
    export r={2};
    sed -i 's/'\"\${s}\"'/'\"\${r}\"'/g' ${f}-rsid.vepinput
  "
  ) | bash
# VEP annotation
  cd ${LOFTEE}
  vep --input_file ${INF}/work/${f}-rsid.vepinput --output_file ${INF}/work/${f}-rsid.tab --force_overwrite \
      --cache --dir_cache ${HPC_WORK}/ensembl-vep/.vep --dir_plugins ${HPC_WORK}/loftee --offline \
      --species homo_sapiens --assembly GRCh37 --pick --nearest symbol --symbol --plugin TSSDistance \
      --plugin LoF,loftee_path:.,human_ancestor_fa:human_ancestor.fa.gz,conservation_file:phylocsf_gerp.sql.gz \
      --tab
done
# per_gene
sed 's/chr//'  INF1.merge.cis | \
grep -f - INF1.merge-rsid.vepinput > INF1.merge.cis.vcf
export s=INF1.merge.cis
vep -i ${s}.vcf -o ${s}.vepoutput --per_gene --check_existing --force_overwrite --offline \
    --species homo_sapiens --everything --assembly GRCh37 --nearest symbol \
    --symbol --pubmed --uniprot --protein --sift b --polyphen b --tab
## cross-check
grep -v '#' ${s}.vepoutput | \
grep -e missense_variant -e stop_lost -e splice_acceptor_variant | cut -f1 | grep -f - INF1.merge.cis
cd ${dir}

### clumping

cat work/INF1.merge.alleles | \
tr '\t' ' ' | \
parallel -C' ' --env INF '
  echo {2}
  plink --bfile ${INF}/INTERVAL/cardio/INTERVAL \
        --snp {2} --window 2000 \
        --clump METAL/gwas2vcf/{1}.tsv.gz \
        --clump-snp-field snpid \
        --clump-p1 1 --clump-p2 1 --clump-r2 0.8 --clump-field p \
        --out ${INF}/annotate/{1}-{2}
  awk "NR>1{print \$3}" {1}-{2}.clumped > ${INF}/annotate/{1}-{2}.snpid
  export exist=$(grep {2} ${INF}/annotate/{1}-{2}.snpid)
  if [ ${exist} =="" ]; then echo {1} >> ${INF}/annotate/{1}-{2}.snpid; fi
  (
    echo "##fileformat=VCFv4.0"
    echo "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO" | \
    tr "," "\t"
    awk "{print \$3}" {1}-{2}.snpid | \
    awk "{snpid=\$1;gsub(/chr|:|_/,\" \",snpid);split(snpid,a,\" \");print a[1],a[2],\$1,a[3],a[4],\".\",\".\",\".\"}" | \
    sort -k1,1n -k2,2n | \
    tr " " "\t"
  ) > ${INF}/annotate/{1}-{2}.vcf
'
