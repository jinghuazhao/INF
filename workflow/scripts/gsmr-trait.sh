function efo_update()
# switch to gsmr_trait in mr.sb apeared problematic
{
  export EFO_UPDATE=${INF}/OpenGWAS/efo-update.txt
  sed '1d' ${EFO_UPDATE} | grep -e ebi -e ieu | awk -vFS="\t" '{print $1,2/(1/$3+1/$4)}' | \
  while read efo N
  do
    export efo=${efo}
    export N=${N}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | awk '!/TNFB/' | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 "{print \$1, \$2, \$3\":\"\$4-M\"-\"\$5+M}" | \
    parallel -C' ' -j15 --env INF --env efo --env N --env suffix '
      (
        echo -e "SNP A1 A2 freq b se p N"
        bcftools query -f "%CHROM %POS %ID %ALT %REF [%AF] [%ES] [%SE] [%LP] [%SS] \n" -r {3} ${INF}/OpenGWAS/${efo}.vcf.gz | \
        awk -vN=${N} "{
               if(\$4<\$5) snpid=\"chr\"\$1\":\"\$2\"_\"\$4\"_\"\$5;
               else snpid=\"chr\"\$1\":\"\$2\"_\"\$5\"_\"\$4
               \$9=10^-\$9
               if (\$10==\".\") \$10=N
               print snpid, \$1, \$2, \$4, \$5, \$6, \$7, \$8, \$9, \$10
             }" | sort -k2,2n -k3,3n -k9,9gr | cut -d" " -f2,3 --complement | awk "a[\$1]++==0"
      ) | \
      gzip -f > ${INF}/mr/gsmr/trait/{2}-${efo}.gz
    '
  done
  sed '1d' ${EFO_UPDATE} | grep -v -e ebi -e ieu -e finn | awk -vFS="\t" '{print $1,2/(1/$3+1/$4)}' | \
  while read efo N
  do
    export efo=${efo}
    export N=${N}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | awk '!/TNFB/' | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 "{print \$1, \$2, \$3\":\"\$4-M\"-\"\$5+M}" | \
    parallel -C' ' -j15 --env INF --env efo --env cases --env controls '
      echo ${efo} {2}
      (
        echo SNP A1 A2 freq b se p N
        case ${efo} in
        GCST90019016)
           tabix ${INF}/OpenGWAS/${efo}/harmonised/${efo}.tsv.gz {3} | \
           awk -vefo=${efo} -vN=${N} "
           {
               chr=\$3; pos=\$4; a1=\$6; a2=\$5
               if (a1<a2) snpid=\"chr\" chr \":\" pos \"_\" a1 \"_\" a2; else snpid=\"chr\" chr \":\" pos \"_\" a2 \"_\" a1
               print snpid, a1, a2, \$8, \$7, \$9, \$10, N
           }" | \
           sort -k1,1 -k7,7g | \
           awk "a[\$1]++==0"
           ;;
        GCST90010715)
           tabix ${INF}/OpenGWAS/${efo}/harmonised/${efo}.tsv.gz {3} | \
           awk -vefo=${efo} -vN=${N} "
           {
               chr=\$3; pos=\$4; a1=\$5; a2=\$6
               if (a1<a2) snpid=\"chr\" chr \":\" pos \"_\" a1 \"_\" a2; else snpid=\"chr\" chr \":\" pos \"_\" a2 \"_\" a1
               print snpid, a1, a2, \$7, \$8, \$9, \$10, \$11
           }" | \
           sort -k1,1 -k7,7g | \
           awk "a[\$1]++==0"
           ;;
        GCST90061440 | GCST90014325)
           tabix ${INF}/OpenGWAS/${efo}/harmonised/${efo}.tsv.gz {3} | \
           awk -vefo=${efo} -vN=${N} "
           {
               chr=\$3; pos=\$4; a1=\$6; a2=\$5
               if (a1<a2) snpid=\"chr\" chr \":\" pos \"_\" a1 \"_\" a2; else snpid=\"chr\" chr \":\" pos \"_\" a2 \"_\" a1
               print snpid, a1, a2, \$8, \$7, \$10, \$9, N
           }" | \
           sort -k1,1 -k7,7g | \
           awk "a[\$1]++==0"
           ;;
        GCST90014023)
           tabix ${INF}/OpenGWAS/${efo}/harmonised/${efo}.tsv.gz {3} | \
           awk -vefo=${efo} -vN=${N} "
           {
               chr=\$1; pos=\$2; a1=\$7; a2=\$6
               if (a1<a2) snpid=\"chr\" chr \":\" pos \"_\" a1 \"_\" a2; else snpid=\"chr\" chr \":\" pos \"_\" a2 \"_\" a1
               print snpid, \$7, \$6, \$9, \$8, \$11, \$10, N
           }" | \
           sort -k1,1 -k7,7g | \
           awk "a[\$1]++==0"
           ;;
        GCST90134602)
           tabix ${INF}/OpenGWAS/${efo}/${efo}.tsv.gz {3} | \
           awk "{
               chr=\$1;pos=\$2;a1=\$4;a2=\$3;
               if (a1<a2) snpid=\"chr\" chr \":\" pos \"_\" a1 \"_\" a2; else snpid=\"chr\" chr \":\" pos \"_\" a2 \"_\" a1
               print snpid, a1, a2, \$11, \$5, \$6, \$7, \$10
           }"
        ;;
      # chr pos A2 A1 b se p cases controls N freq rsid
        *)
        ;;
        esac
      ) | \
      gzip -f > ${INF}/mr/gsmr/trait/{2}-${efo}.gz
    '
  done
  sed '1d' ${EFO_UPDATE} | grep finn | awk -vFS="\t" '{print $1,2/(1/$3+1/$4)}' | \
  while read efo N
  do
    export efo=${efo}
    export N=${N}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 "{print \$1, \$2, \$3\":\"\$4-M\"-\"\$5+M}" | \
    parallel -C' ' -j15 --env INF --env efo --env N '
      echo ${efo} {2} {3}
      cat <(echo SNP A1 A2 freq b se p N) \
          <(tabix ${INF}/OpenGWAS/${efo}.gz {3} | \
            awk -vN=${N} "{
               chr=\$1;pos=\$2;a1=toupper(\$4);a2=toupper(\$3)
               if (a1<a2) snpid=\"chr\" chr \":\" pos \"_\" a1 \"_\" a2;
                     else snpid=\"chr\" chr \":\" pos \"_\" a2 \"_\" a1
               print snpid,a1,a2,\$11,\$9,\$10,\$7,N
               }" | \
            sort -k1,1 -k7,7g | \
            awk "a[\$1]++==0"
           ) | \
      gzip -f > ${INF}/mr/gsmr/trait/{2}-${efo}.gz
    '
  done
# chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls
  sed '1d' ${EFO_UPDATE} | awk -vFS="\t" '{print $1,2/(1/$3+1/$4)}' | \
  while read efo N
  do
    export efo=${efo}
    export N=${N}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    parallel -C' ' -j15 --env INF --env efo '
      echo ${efo} ${INF}/mr/gsmr/trait/{2}-${efo}.gz > ${INF}/mr/gsmr/trait/gsmr_{2}-${efo}
    '
  done
}

efo_update
