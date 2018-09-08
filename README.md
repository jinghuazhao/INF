# INF
SCALLOP-INF analysis

The list of proteins on inflammation is within one of the panels from [OLINK](https://www.olink.com/products/), each containing 92 proteins as in [Olink validation data all panels.xlsx](doc/Olink%20validation%20data%20all%20panels.xlsx). Information contained in these panels can be retrieved into R via [OLINK.R](doc/OLINK.R), which also attempts to compromise earlier version and annotations.

File / root  | Description
---------|---------------------------------------------------------------------------------------------------------------------
[inf1.csv](doc/inf1.csv) | UniProt ID, protein, target for the INF panel
[log10p.md](doc/log10p.md) | On a very competitive log10(p) calculator
MMP12 | Bash ([MMP12.sh](doc/MMP12.sh)) and R ([MMP12.R](doc/MMP12.R)) scripts for MMP12-CHD analysis
SERPINA1 | [SERPINA1.R](doc/SERPINA1.R) and [SERPINA1.pdf](doc/SERPINA1.pdf) which emulate [original version](doc/fig2.R) of the SOMAscan paper (Sun et al. 2018)
SOMAscan | Methods ([SOMAscan.md](doc/SOMAscan.md), [SOMAscan.pdf](doc/SOMAscan.pdf), [SOMAscan.tif](doc/SOMAscan.tif)) and supplementary tables ([SOMAscan.xlsx](doc/SOMAscan.xlsx)) for the SOMAscan paper
SecureCloud.md | Information for SecureCloud (unavailable here yet)

As noted in MMP12.sh, the MMP12 case could have been done genomewide; the CAD summary statistics is described 
[here](https://github.com/jinghuazhao/Omics-analysis/tree/master/CAD). A colocalisation analysis on simulated 
data can be found in the association analysis section of [software-notes](https://github.com/jinghuazhao/software-notes) as 
well as the [BMI example](https://github.com/jinghuazhao/Omics-analysis/tree/master/BMI).

## Notes on UniProt IDs

The use of UniProt IDs is rationalised in two aspects,

1. The protein list in [inf1.csv](doc/inf1.csv) found O43508 and Q4ACW9 are associated with the TNFSF12 and TWEAK gene, respectively.

2. Q8NF90 and Q8WWJ7 were not listed at the UCSC, their availability on UniProt seem to be for backward compatibility as on query they 
point to P12034 and P30203 (Q8WWJ7_HUMAN should have been CD6_HUMAN). [hgTables.txt](doc/hgTables.txt) is based on UCSC, checked over
UniProt IDs as follows,
```bash
sort -k2,2 inf1.list > 1
awk '{FS=OFS="\t"; split($4,f,"-");$4=f[1];if(!index($1,"_")) print}' hgTables.txt | sort -k4,4 > 2
join -t$'\t' -12 -24 1 2 > 12
# 90 lines
wc -l 12
# Q8NF90 (FGF.5), Q8WWJ7 (CD6) are missing
join -v2 -22 12 1
rm 1 2 12
```
A UniProt ID may be associated with multiple chromosomes, e.g., Q6IEY1 with chromosomes 1 and 5. While [inf1.csv](doc/inf1.csv) 
edits Q4ACW9, [inf2.csv](doc/inf2.csv) is inline with UCSC with respect to P12034 and P30203.

## References

Stacey D, et al. (2017), [ProGeM](https://github.com/ds763/ProGeM): A framework for the prioritisation of candidate causal genes at molecular 
quantitative trait loci, http://dx.doi.org/10.1101/230094.

Sun BB, et al. (2018). Genomic atlas of the human plasma proteome. *Nature* 558: 73–79
