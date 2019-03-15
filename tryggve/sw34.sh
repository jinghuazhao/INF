#!/usr/bin/bash

function plink19()
{
  export module=plink2/1.90beta5.4
  export c=plink
  export d=import-dosage
  export p=fam
  export m=map
}
function plink20()
{
  export module=plink2/2.00alpha20181028
  export c=plink2
  export d=import-dosage
  export p=psam
  export m=map
}

plink20
module load $module

export rt=/data/andmala/STANLEY
ls $rt/ | \
grep .gz | \
grep -v test | \
sed 's/dos_bip_sw34_eur_rk-qc.hg19.ch.fl.chr//g;s/.out.dosage.fam//g;s/.out.dosage.gz//g;s/.out.dosage.map//g' > sw34.list
cat sw34.list | \
parallel -j5 --env rt --env c --env d --env p --env m -C' ' '
  $c --${d} $rt/dos_bip_sw34_eur_rk-qc.hg19.ch.fl.chr{}.out.dosage.gz \
     --${p} $rt/dos_bip_sw34_eur_rk-qc.hg19.ch.fl.chr{}.out.dosage.fam \
     --${m} $rt/dos_bip_sw34_eur_rk-qc.hg19.ch.fl.chr{}.out.dosage.map \
     --export bgen-1.1 --out {}
'

module unload $module
module load bgen/20180807
cat-bgen -g $(awk '{$1=$1 ".bgen";printf $1 " "}' sw34.list) -og sw34
