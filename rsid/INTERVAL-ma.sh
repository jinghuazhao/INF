#!/usr/bin/bash

module load ceuadmin/stata/15
ls $SCALLOP/jp549/olink-merged-output/INTERVAL_inf1*gz | grep inf1 | xargs -l basename | \
sed 's|INTERVAL_inf1_||g;s|___| |g;s|_chr_merged.gz||g;s|\*||g' | \
cut -d' ' -f1,2 | head -2 | \
parallel --env INF --env SCALLOP -C' ' '
   echo {1}-{2}
   export f=${SCALLOP}/jp549/olink-merged-output/INTERVAL_inf1_{1}___{2}_chr_merged.gz | \
   (
     zcat ${f} | \
     head -1 | awk -vOFS="\t" "{\$1=\"Chrom\";\$2=\"Start\" \"\t\" \"End\" \"\t\" \"MarkerName\";print}"
     zcat ${f} | \
     awk -vOFS="\t" "(NR>1 && \$24!=\"NA\" && \$25 != \"NA\" && \$22 < 1-e3) {
       chr=\$3+0;a1=toupper(\$6);a2=toupper(\$5);
       if(a1<a2) snpid=\"chr\" chr \":\" \$4 \"_\" a1 \"_\" a2; else snpid=\"chr\" chr \":\" \$4 \"_\" a2 \"_\" a1;
       if (\$2!=\".\") MarkerName=\$2; else MarkerName=snpid;
       \$1="chr" chr; \$2=\$4-1 \"\t\" \$4 \"\t\" MarkerName;print}" | \
     awk "!(\$1 == \"chr6\" && \$3 >= 25392021 && \$3 < 33392022)"
   ) > {1}.a
   export p0=$(awk "NR>1&&\$22==0" {1}.a | wc -l | cut -d" " -f1)
   (
      echo SNP A1 A2 freq b se p N
      if [ ${p0} -eq 0 ]; then
         awk "NR>1{print \$4,toupper(\$6),toupper(\$5),\$19,\$24,\$25,\$22,\$18}" {1}.a
         rm {1}.a
      else
         awk -vOFS="\t" "{if(NR==1) print \$0, \"N\"; else print \$0,NR-1}" ${INF}/work/{1}.merged > {1}.b
         bedtools intersect -wb -a {1}.a -b {1}.b > {1}.ab
         stata -b -q ${INF}/csd3/ma.do 
         awk -f ${INF}/csd3/ma.awk {1}.txt
         rm {1}.a {1}.b {1}.ab {1}.txt {1}.dta {1}0.dta
      fi
   ) > INTERVAL/{1}.ma
    export lines=$(wc -l INTERVAL/{1}.ma | cut -d" " -f1)
    if [ $lines -eq 1 ]; then
      echo removing {1}.ma with $lines lines
      rm INTERVAL/{1}.ma
    fi
'
