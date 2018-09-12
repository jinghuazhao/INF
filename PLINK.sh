# 3-9-2018 JHZ

module load plink2/1.90beta5.4

plink-1.9 --bfile $bfile --clump $rt.tab \
             --clump-snp-field MarkerName --clump-field P.value \
             --clump-kb 500 \
             --clump-p1 5e-08 --clump-p2 0.01 \
             --clump-r2 0.1 \
             --clump-snp-field snpid \
--out $rt.clump
