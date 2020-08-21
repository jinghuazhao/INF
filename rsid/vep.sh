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
      --species homo_sapiens --assembly GRCh37 --pick --symbol --plugin TSSDistance \
      --plugin LoF,loftee_path:.,human_ancestor_fa:human_ancestor.fa.gz,conservation_file:phylocsf_gerp.sql.gz \
      --tab
done
cd ${dir}

sed 's/chr//'  work/INF1.merge.cis | \
grep -f - work/INF1.merge-rsid.vepinput > INF1.merge.cis.vcf
export s=INF1.merge.cis
vep -i ${s}.vcf -o ${s}.vepoutput --per_gene --check_existing --distance 500000 --force_overwrite --offline \
    --species homo_sapiens --everything --assembly GRCh37 \
    --symbol --pubmed --uniprot --protein --sift b --polyphen b --tab
