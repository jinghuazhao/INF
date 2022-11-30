function mr_rsid()
# data with RSid for TwoSampleMR
{
  for r in {1..59}
  do
    export region=$(awk -vr=${r} 'NR==r{print $4":"$5"-"$6}' ${INF}/TNFB/cis.dat)
    export prot=$(awk -vr=${r} 'NR==r{print $3}' ${INF}/TNFB/cis.dat)
    echo ${prot} ${region}
    if [ ${prot} == "TNFB" ]; then continue; fi
    for OpenGWAS in $(sed '1d' ${INF}/OpenGWAS/efo-update.txt | cut -f1 | awk '/ieu|ebi|bbj/')
    do
      export N=$(grep -w ${OpenGWAS} ${INF}/OpenGWAS/efo-update.txt | awk -vFS='\t' '{print 2/(1/$3+1/$4)}')
      (
        echo -e "SNP A1 A2 freq b se p N"
        bcftools query -f "%ID %ALT %REF [%AF] [%ES] [%SE] [%LP] [%SS]\n" -r ${region} ${INF}/OpenGWAS/${OpenGWAS}.vcf.gz | \
        awk -vN=${N} '{$7=10^-$7;if ($8==".") $8=N;print}'
      ) > ${INF}/mr/gsmr/trait/${prot}-${OpenGWAS}-rsid.txt
    done
    for GCST in $(sed '1d' ${INF}/OpenGWAS/efo-update.txt | cut -f1 | awk '/^GCST/')
    do
      export N=$(grep -w ${GCST} ${INF}/OpenGWAS/efo-update.txt | awk -vFS='\t' '{print 2/(1/$3+1/$4)}')
      (
        echo -e "SNP A1 A2 freq b se p N"
        case ${GCST} in
        GCST90010715)
          tabix ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.tsv.gz ${region} | \
          awk '{print $2,$5,$6,$7,$8,$9,$10,$11}'
          ;;
        GCST90014325 | GCST90061440)
          tabix ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.tsv.gz ${region} | \
          awk -vN=${N} '{print $2,$6,$5,$8,$7,$10,$9,N}'
          ;;
        GCST90019016)
          tabix ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.tsv.gz ${region} | \
          awk -vN=${N} '{print $2,$6,$5,$8,$7,$9,$10,N}'
          ;;
        GCST90014023)
          tabix ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.tsv.gz ${region} | \
          awk -vN=${N} '{print $3,$7,$6,$9,$8,$11,$10,N}'
          ;;
        GCST90134602)
          tabix ${INF}/OpenGWAS/${GCST}/${GCST}.tsv.gz ${region} | \
          awk '{print $12,$4,$3,$11,$5,$6,$7,$10}'
          ;;
        *)
          ;;
        esac
      ) > ${INF}/mr/gsmr/trait/${prot}-${GCST}-rsid.txt
    done
    for finngen in $(sed '1d' ${INF}/OpenGWAS/efo-update.txt | cut -f1 | awk -vFS='\t' '/^finngen/')
    do
      export N=$(grep -w ${finngen} ${INF}/OpenGWAS/efo-update.txt | awk -vFS='\t' '{print 2/(1/$3+1/$4)}')
      (
         echo -e "SNP A1 A2 freq b se p N"
         tabix ${INF}/OpenGWAS/${finngen}.gz ${region} | \
         awk -vN=${N} '{print $5,$4,$3,$11,$9,$10,$7,N}'
      ) > ${INF}/mr/gsmr/trait/${prot}-${finngen}-rsid.txt
    done
  done
}

mr_rsid
