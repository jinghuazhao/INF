#!/usr/bin/bash

export job=${1}
export list=${INF}/work/INF1.merge-rsid
export p=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $5}')
export r=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $6}')
export pr=${p}-${r}
export chr=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $8}')
export pos=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $9}')
export flanking=1e6
export start=$(sed '1d' ${list} | awk -v job=${job} -v d=${flanking} 'NR==job{start=$9-d;if(start<0) start=0; print start}')
export end=$(sed '1d' ${list} | awk -v job=${job} -v d=${flanking} 'NR==job{print $9+d}')
export sumstats=$SCALLOP/jp549/olink-merged-output/$(ls $SCALLOP/jp549/olink-merged-output | grep -v cvd | grep ${p}___)
export sample=${INF}/INTERVAL/o5000-inf1-outlier_in-r2.sample
export study=INTERVAL
export N=4994
export n=1000
export TMPDIR=${HPC_WORK}/work
export dir=${INF}/sentinels
# scaled results
export list=${INF}/work/INF1.merge
export p=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $5}')
export r=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $6}')
export sumstats=${INF}/finemapping/${p}-${r}.out
export sample=${INF}/finemapping/${p}-${r}.sample

cd ${INF}/cs
# z0 - INTERVAL data only excluding NLRP2 region if overlapped and employs P value cutoff to avoid nonsense results
(
  echo MarkerName Chromosome Position Allele2 Allele1 MAF Effect StdErr
# zcat ${sumstats} | \
# awk 'NR > 1' | \
  cat ${sumstats} | \
  awk 'NR > 17' | \
  awk -vchr=${chr} -vstart=${start} -vend=${end} '
  {
    if (($3==chr && $4 >= start && $4 < end && $22 < 1e-3) && !(chr==19 && $4 >= 53296855 && $4 <= 54500000)) {
       chr=$3+0;a1=toupper($5);a2=toupper($6);
       if(a1<a2) snpid=chr ":" $4 "_" a1 "_" a2; else snpid=chr ":" $4 "_" a2 "_" a1;
       if ($2!=".") MarkerName=$2; else MarkerName=snpid;
       if ($24 != "NA" && $25 != "NA") print MarkerName, chr, $4, a2, a1, $19, $24, $25
    }
  } ' | \
  sort -k1,1 | \
  join - <(sort -k1,1 ${INF}/prune/${pr}.prune) | sort -k2,2n -k3,3n
) > ${pr}.z

export list=${INF}/work/INF1.merge-rsid
export p=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $5}')
export r=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $6}')
R --no-save < ${INF}/rsid/cs.R > /dev/null
cd -
