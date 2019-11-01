## Notes on UniProt IDs

The list of proteins on inflammation is within one of the [OLINK](https://www.olink.com/products/) panels (as in [Olink validation data all panels.xlsx](doc/Olink%20validation%20data%20all%20panels.xlsx)) each containing 92 proteins. Information contained in these panels can be retrieved into R via [OLINK.R](../cardio/OLINK.R), which also attempts to compromise earlier version and annotations.

The use of [UniProt](https://www.uniprot.org/) IDs is rationalised in two aspects,

1. The protein list in [inf1.csv](doc/inf1.csv) found O43508 and Q4ACW9 are associated with the TNFSF12 and TWEAK genes, respectively.

2. Q8NF90 and Q8WWJ7 were not listed at the UCSC, their availability on UniProt seem to be for backward compatibility as on query they 
point to P12034 and P30203 (Q8WWJ7_HUMAN should have been CD6_HUMAN). [hgTables.tsv](hgTables.tsv) is based on UCSC, checked over
UniProt IDs as follows,
```bash
grep doc/inf1 olink.prot.list.txt | \
sed 's/inf1_//g;s/___/\t/g' > inf1.list
sort -k2,2 inf1.list > 1
awk '{
   FS=OFS="\t"; 
   split($4,f,"-");
   $4=f[1];
   if(!index($1,"_")) print
}' doc/hgTables.tsv | \
sort -k4,4 > 2
join -t$'\t' -12 -24 1 2 > 12
# 90 lines
wc -l 12
# Q8NF90 (FGF.5), Q8WWJ7 (CD6) are missing
join -v2 -22 12 1
rm 1 2 12
```
A UniProt ID may be associated with multiple chromosomes, e.g., Q6IEY1 with chromosomes 1 and 5. While [inf1.csv](inf1.csv) 
edits Q4ACW9, [inf2.csv](inf2.csv) is inline with UCSC with respect to P12034 and P30203.

A version by Jimmy is [olink.inf.panel.annot.tsv](olink.inf.panel.annot.tsv) from [olink.annotation.R](olink.annotation.R).

BDNF has recently been removed from the assay and replaced with CD8A, https://www.olink.com/bdnf-info/, and there are also changes on TNF and IFN.gamma, https://www.olink.com/inflammation-upgrade/.
