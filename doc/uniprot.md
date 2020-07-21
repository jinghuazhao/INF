## Notes on UniProt IDs

The list of proteins on inflammation is within one of the [Olink](https://www.olink.com/products/) panels (as in [Olink validation data all panels.xlsx](doc/Olink%20validation%20data%20all%20panels.xlsx)) each containing 92 proteins. Information contained in these panels can be retrieved into R via [Olink.R](Olink.R), which also attempts to compromise earlier version and annotations. Nevertheless O43508 is replaced with Q4ACW9 for TWEAK.

BDNF has recently been removed from the assay and replaced with CD8A, https://www.olink.com/bdnf-info/, and there are also changes on TNF and IFN.gamma, https://www.olink.com/inflammation-upgrade/.

A [UniProt](https://www.uniprot.org/) ID may be associated with multiple chromosomes, e.g., Q6IEY1 with chromosomes 1 and 5. While [inf1.csv](inf1.csv) 
edits Q4ACW9, [inf2.csv](inf2.csv) is inline with UCSC with respect to P12034 and P30203.

The use of uniprot IDs is noted in two aspects,

1. The protein list in [inf1.csv](inf1.csv) notes both O43508 and Q4ACW9.

2. Q8NF90 and Q8WWJ7 were not listed at the UCSC, their availability on UniProt seem to be for backward compatibility as on query they 
point to P12034 and P30203 (Q8WWJ7_HUMAN should have been CD6_HUMAN). [hgTables.tsv](hgTables.tsv) is based on UCSC, checked over
UniProt IDs as follows,
```bash
grep inf1 doc/olink.prot.list.txt | \
sed 's/inf1_//g;s/___/\t/g' > inf1.list
join -t$'\t' -12 -24 <(sort -k2,2 inf1.list) <(awk '{split($4,f,"-"); $4=f[1]; if(!index($1,"_")) print}' OFS='\t' doc/hgTables.tsv | sort -k4,4) > 12
# 90 lines
wc -l 12
# Q8NF90 (FGF.5), Q8WWJ7 (CD6) are missing
join -v2 -22 12 <(sort -k2,2 inf1.list)
rm 12
```
Likewise, [olink.inf.panel.annot.tsv](olink.inf.panel.annot.tsv) from [olink.annotation.R](olink.annotation.R) also has the following two entries

"target" | "target.short" | "uniprot" | "panel" | "prot.on.multiple.panel" | "panels.with.prot" | "hgnc_symbol" | "chromosome_name" | "start_position" | "end_position" | "olink.id" | "alternate.uniprot"
-------------------------------------|----------------|-----------|---------|--------------------------|--------------------|---------------|-------------------|------------------|----------------|------------|--------------------
"Fibroblast growth factor 5 (FGF-5)" | "FGF-5" | "Q8NF90" | "inf" | FALSE | NA | NA | "4" | 81187753 | 81257834 | "141_FGF-5" | "P12034"
"T-cell surface glycoprotein CD6 isoform (CD6)" | "CD6" | "Q8WWJ7" | "inf" | FALSE | NA | NA | "11" | 60739115 | 60787849 | "131_CD6" | "P30203"

whose hgnc_symbol can be amended as follows
```bash
awk '!/BDNF/ && NR > 1 {gsub(/"/,"",$7);if ($3=="Q8NF90") $7="FGF5"; else if ($3=="Q8WWJ7") $7="CD6"};1' FS='\t' OFS='\t' doc/olink.inf.panel.annot.tsv
```
The overlap with SomaLogic panel is characterised with [Olink.R](Olink.R) which also gives a Venn diagram.
<img src="Olink-SomaLogic-Venn-diagram.png" width="300" height="300" align="right">
