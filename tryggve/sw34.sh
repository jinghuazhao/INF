#!/usr/bin/bash

module load plink2/1.90beta5.4

export rt=/data/andmala/STANLEY 
ls $rt | \
sed 's/dos_bip_sw34_eur_rk-qc.hg19.ch.fl.chr//g;s/.out.dosage.fam//g;s/.out.dosage.gz//g;s/.out.dosage.map//g' | \
parallel -j5 --env rt -C' ' '
  plink --dosage $rt/dos_bip_sw34_eur_rk-qc.hg19.ch.fl.chr{}.out.dosage.gz \
        --fam $rt/dos_bip_sw34_eur_rk-qc.hg19.ch.fl.chr{}.out.dosage.fam \
        --map $rt/dos_bip_sw34_eur_rk-qc.hg19.ch.fl.chr{}.out.dosage.map
        --freq --out work/sw34.{}.freq
'
