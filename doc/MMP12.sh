# 31/8/2018 JHZ

# MMP12 (UCSC uc001phk.3) at chr11:102733464-102745764
# MMP12 (RefSeq) at chr11:102733460-102745764

# One should use whole-genome counterpart to replicate analysis in the paper
gunzip -c MMP12.4496.60.2/MMP12.4496.60.2_chrom_11_meta_final_v1.tsv.gz | awk 'NR==1||($3>=102733464 && $3<=102745764)' > MMP12.txt
awk 'NR==1||($2==11 && $3>=102733464 && $3<=102745764)' cad.add.160614.website.txt > CAD.txt

/c/Users/jhz22/R-3.5.1/bin/R --no-save  < MMP12.R

# PhenoScanner

awk '(NR>1) {print "chr" $2 ":" $3}' MMP12.txt > MMP12.ps.txt
awk '(NR>1) {print $1}' CAD.txt > CAD.ps.txt
