#!/usr/bin/bash

for trait in CAD FEV1
do
  echo ${trait}
  if [ -f ${INF}/gsmr/INF1_${trait}.gsmr ]; then rm ${INF}/gsmr/INF1_${trait}.gsmr; fi
  (
    cat ${INF}/gsmr/gsmr_${trait}*.gsmr | \
    head -1
    ls ${INF}/gsmr/gsmr_${trait}*gsmr | \
    parallel -j1 -C' ' '
      if [ -f {} ]; then
         awk "NR>1" {}
      fi
    '
  ) | \
  grep -v nan > ${INF}/gsmr/INF1_${trait}.gsmr
done

join -a1 -e "NA" -o1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8 \
     <(gunzip -c ${INF}/gsmr/gsmr_FEV1_LIF.R.eff_plot.gz | \
       awk '/effect_begin/,/effect_end/' | \
       grep -v effect | \
       sort -k1,1) \
     <(gunzip -c ${INF}/gsmr/gsmr_CAD_LIF.R.eff_plot.gz | \
       awk '/effect_begin/,/effect_end/' | \
       grep -v effect | \
       sort -k1,1) | \
join work/INTERVAL.rsid - > ${INF}/gsmr/INF1_CAD-FEV1.snp_effects

(
  echo SNP b.LIF.R SE.LIF.R b.FEV1 SE.FEV1 b.CAD SE.CAD
  cut -d' ' -f2,6,7,8,9,16,17 ${INF}/gsmr/INF1_CAD-FEV1.snp_effects
) > ${INF}/gsmr/INF1_CAD-FEV1.gsmr

for trait in A2 B1 B2 C2
do
  if [ -f ${INF}/HGI/INF1_${trait}.gsmr ]; then rm ${INF}/HGI/INF1_${trait}.gsmr; fi
  (
    cat ${INF}/HGI/gsmr_${trait}*.gsmr | \
    head -1
    ls ${INF}/HGI/gsmr_${trait}*gsmr | \
    parallel -j1 -C' ' '
      if [ -f {} ]; then
         awk "NR>1" {}
      fi
    '
  ) | \
  grep -v nan > ${INF}/HGI/INF1_${trait}.gsmr
done

(
  echo prot uniprot A2 b_A2 se_A2 p_A2 n_A2 B2 b_B2 se_B2 p_B2 n_B2 C2 b_C2 se_C2 p_C2 n_C2
  join -a1 -e "NA"   -o1.1,1.2,2.2,2.3,2.4,2.5,2.6 <(sort -k1,1 ${INF}/work/inf1.tmp) <(sed '1d' ${INF}/HGI/INF1_A2.gsmr | cut -f1-6 | sort -k1,1) | \
  join -a1 -e "NA" - -o1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6 <(sed '1d' ${INF}/HGI/INF1_B2.gsmr | cut -f1-6 | sort -k1,1 ) | \
  join -a1 -e "NA" - -o1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.2,2.3,2.4,2.5,2.6 \
                     <(sed '1d' ${INF}/HGI/INF1_C2.gsmr | cut -f1-6 | sort -k1,1) | \
  awk '$3!="NA"'
) > ${INF}/HGI/A2-B2-C2.txt

for trait in A2 B2 C2
do
  gunzip -c ${INF}/HGI/gsmr_${trait}_LIF.R.eff_plot.gz | \
  awk '/effect_begin/,/effect_end/' | \
  grep -v effect > ${INF}/HGI/gsmr_${trait}_LIF.R.dat
done

R --no-save -q <<END
  options(width=200)
  INF <- Sys.getenv("INF")
  A2 <- read.table(file.path(INF,"HGI","gsmr_A2_LIF.R.dat"),col.names=c("snpid","A2.a1","A2.a2",paste0("x",1:5)))
  B2 <- read.table(file.path(INF,"HGI","gsmr_B2_LIF.R.dat"),col.names=c("snpid","B2.a1","B2.a2",paste0("y",1:5)))
  C2 <- read.table(file.path(INF,"HGI","gsmr_C2_LIF.R.dat"),col.names=c("snpid","C2.a1","C2.a2",paste0("z",1:5)))
  A2_B2_C2 <- merge(merge(A2,B2,by="snpid",all=TRUE),C2,by="snpid",all=TRUE)
  vars <- c("snpid",paste0("x",1:5),paste0("y",4:5),paste0("z",4:5))
  snp_effects <- A2_B2_C2[vars]
  names(snp_effects) <- c("SNP", "MAF.LIF.R", "b.LIF.R", "SE.LIF.R", "b.A2", "SE.A2", "b.B2", "SE.B2", "b.C2", "SE.C2")
  write.table(snp_effects,file=file.path(INF,"HGI","INF1_A2-B2-C2.gsmr"),row.names=FALSE,quote=FALSE)
END

for ext in png tif
do
  convert ${INF}/gsmr/INF1_CAD-FEV1.pdf ${INF}/gsmr/INF1_CAD-FEV1.${ext}
  convert ${INF}/HGI/INF1_A2-B2-C2.pdf ${INF}/HGI/INF1_A2-B2-C2.${ext}
done
