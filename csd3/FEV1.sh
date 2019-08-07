# 7-8-2019 JHZ

(
  echo "SNP A1 A2 freq b se p N"
  gunzip -c UKB_FEV1_results.txt.gz | awk '
  {
    CHR=$2
    POS=$3
    a1=$4
    a2=$5
    if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
    else snpid="chr" CHR ":" POS "_" a1 "_" a2
    if (NR>1) print snpid, a1, a2, $7, $8, $9, $10, 321047
  }'
) > gsmr_FEV1.txt

# gunzip -c UKB_FEV1_results.txt.gz | head -1 | sed 's/\t/\n/g' | awk '{print "#" NR, $1}'
#1 #SNP
#2 Chromosome
#3 Position_b37
#4 Coded
#5 Non_coded
#6 INFO
#7 Coded_freq
#8 beta
#9 SE_GC
#10 P_GC

# gunzip -c SpiroMeta_FEV1_results.txt.gz | head -1 | sed 's/\t/\n/g' | awk '{print "#" NR, $1}'
#1 #SNP
#2 Chromosome
#3 Position_b37
#4 Coded
#5 Non_coded
#6 N
#7 Neff
#8 Coded_freq
#9 beta
#10 SE
#11 P
