#!/usr/bin/bash

function plink1()
{
  export module=plink2/1.90beta5.4
  export c=plink
  export d=import-dosage
  export p=fam
  export m=map
}
function plink2()
{
  export module=plink2/2.00alpha20181028
  export c=plink2
  export d=import-dosage
  export p=psam
  export m=map
}

module load bgen/20180807
module load $module
export rt=/data/andmala/STANLEY

plink2

ls $rt/*gz | \
sed 's/dos_bip_sw34_eur_rk-qc.hg19.ch.fl.chr//g;s/.out.dosage.fam//g;s/.out.dosage.gz//g;s/.out.dosage.map//g' | \
head -2 | \
parallel -j5 --env rt --env c --env d --env p --env m -C' ' '
  $c --${d} $rt/dos_bip_sw34_eur_rk-qc.hg19.ch.fl.chr{}.out.dosage.gz \
     --${p} $rt/dos_bip_sw34_eur_rk-qc.hg19.ch.fl.chr{}.out.dosage.fam \
     --${m} $rt/dos_bip_sw34_eur_rk-qc.hg19.ch.fl.chr{}.out.dosage.map \
     --bgen --log {} --out {}
'
module unload $module
