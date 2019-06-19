<img src="doc/INF1.circlize.png" width="300" height="300" align="right">

---

# SCALLOP Olink/INFlammation analysis

---

## TABLE OF CONTENTS

* [Analysis](https://github.com/jinghuazhao/INF#analysis)
* [References](https://github.com/jinghuazhao/INF#references)
* [Use of approximately independent LD blocks](https://github.com/jinghuazhao/INF#use-of-approximately-independent-ld-blocks)
* [Notes on UniProt IDs](https://github.com/jinghuazhao/INF#notes-on-uniprot-ids)
* [Files, URLs, and downloading](https://github.com/jinghuazhao/INF#files-urls-and-downloading)

<img src="doc/OPG-qmlf.png" width="400" height="800" align="right">

## Analysis

1. To implement the analysis plan, we started with analysis on INTERVAL as with [INTERVAL.sh](tryggve/INTERVAL.sh) and [cardio.sh](cardio/cardio.sh), which includes cis/trans classification via Jimmy's [cis.vs.trans.classification.R](cardio/cis.vs.trans.classification.R) and a Bash version via bedtools, [cistrans.sh](cardio/cistrans.sh).
Jimmy's competitive log10(p) calculator is also documented in [log10p.md](doc/log10p.md) in comparison with R/Rmpfr. 

2. Data preprocessig was done with [list.sh](tryggve/list.sh) and [format.sh](tryggve/format.sh). The meta-analysis was according to [metal.sh](tryggve/metal.sh) using METAL whose results were cross-examined with [QCGWAS.sh](tryggve/QCGWAS.sh) together with addtional investigation.

3. The main analysis follows [analysis.sh](tryggve/analysis.sh), which contains codes for Q-Q/Manhattan/LocusZoom/Forest plots (see the figure on the right for the OPG example, which replicated results of Kwan et al. (2014) as identified by PhenoScanner), clumping using PLINK and conditional analysis using GCTA. The clumping results were classified into cis/trans signals. As the meta-analysis stabilised especially with INTERVAL reference, analysis has increasingly been done locally with [cardio](cardio). In particular, approximately independent LD blocks are utilised as descrbed below.

4. Further downstream analysis will be considered. The CAD summary statistics used for MAGMA and MR is described [here](https://github.com/jinghuazhao/Omics-analysis/tree/master/CAD)
-- as noted in MMP12.sh, the MMP12 case could have been done genomewide. A colocalisation analysis on simulated data can be found in the
association analysis section of [software-notes](https://github.com/jinghuazhao/software-notes)
as well as the [BMI example](https://github.com/jinghuazhao/Omics-analysis/tree/master/BMI).

5. TRYGGVE-specific issues were noted in [tryggve.md](tryggve.md). The `cis.vs.trans.classification`, `circos.cis.vs.trans.plot` (for the circular plot at the top of the page); together with `log10p`, `gc.lambda`, `invnormal`, `METAL_forestplot`, `mhtplot.trunc` functions they are now with gap at CRAN and latest updates at [R/gap](https://github.com/jinghuazhao/R/tree/master/gap) repository.

## References

Folkersen L, et al. (2017). Mapping of 79 loci for 83 plasma protein biomarkers in cardiovascular disease. *PLoS Genetics* 13(4), doi.org/10.1371/journal.pgen.1006706.

Kwan JSH, et al. (2014). Meta-analysis of genome-wide association studies identiﬁes two loci associated with circulating osteoprotegerin levels. *Hum Mol Genet* 23(24): 6684–6693.

Niewczas MA, et al. (2019). A signature of circulating inflammatory proteins and development of end-stage renal disease in diabetes. *Nat Med*. https://doi.org/10.1038/s41591-019-0415-5

Sun BB, et al. (2018). Genomic atlas of the human plasma proteome. *Nature* 558: 73–79.

---

## Use of approximately independent LD blocks

The implementation involves several steps:

1. Set up 1703 autosomal regions as defined in [EURLD.bed](tryggve/EURLD.bed).
2. Extract variants outside the [regions in high LD](tryggve/high-LD-regions-hg19.txt) to [1672 regions](tryggve/EURLD-no-high-LD-regions-hg19.bed) by [EURLD.sh](tryggve/EURLD.sh).
3. Overlap regions and GWAS sumstats:
   * Tag GWAS sumstats with regions through [aild-rma.sb](cardio/aild-rma.sb).
   * Pair protein-region which contains genomewide significant signals by [aild-list.sb](cardio/aild-list.sb).
   * Independently, list variants by region in the reference panel by [aild-snplist.sb](cardio/aild-snplist.sb).
4. clump via [aild-clump.sb](cardio/aild-clump.sb).
5. cojo via [aild-cojo.sb](cardio/aild-cojo.sb).
6. Downstream analyses with PhenoScanner (preferably v2) as in [ps.sh](cardio/ps.sh) and forest plots on TRYGGVE (with study-specific sumstats).

The regions are predefined. As shown in [EURLD.tsv](tryggve/EURLD.tsv) by [EURLD.R](tryggve/EURLD.R), the LD patterns across the genome are more variable than the norm in a typical genomewide association analysis therefore slide windows such as 250kb (36), 500kb (300), or even 10Mb (1071), seeing that the sentinel variant may not necessarily lie right in the middle of a window. The number of signals in our case were close to GCTA but overestimated (53 by PLINK) as in [clump-cojo.md](cardio/clump-cojo.md). For instance, it is often seen from the PLINK --clump-range outputs that sliding windows can give results in two neighbouring LD blocks.

Note that pairing regions of interest would reduce the burden of genomewide analysis, and also that region-specific reference will not affect results from steps 4 and 5 regarding use of variants from GWAS sumstats.

Steps 4 and 5 both use `INF1.aild`, which contains all the protein-region pairs. The results are classified as in [analysis.sh](tryggve/analysis.sh). In particular, for step 5 this is done with [aild.sh](cardio/aild.sh).

## Notes on UniProt IDs

The list of proteins on inflammation is within one of the [OLINK](https://www.olink.com/products/) panels (as in [Olink validation data all panels.xlsx](doc/Olink%20validation%20data%20all%20panels.xlsx)) each containing 92 proteins. Information contained in these panels can be retrieved into R via [OLINK.R](cardio/OLINK.R), which also attempts to compromise earlier version and annotations.

The use of [UniProt](https://www.uniprot.org/) IDs is rationalised in two aspects,

1. The protein list in [inf1.csv](doc/inf1.csv) found O43508 and Q4ACW9 are associated with the TNFSF12 and TWEAK genes, respectively.

2. Q8NF90 and Q8WWJ7 were not listed at the UCSC, their availability on UniProt seem to be for backward compatibility as on query they 
point to P12034 and P30203 (Q8WWJ7_HUMAN should have been CD6_HUMAN). [hgTables.tsv](doc/hgTables.tsv) is based on UCSC, checked over
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
A UniProt ID may be associated with multiple chromosomes, e.g., Q6IEY1 with chromosomes 1 and 5. While [inf1.csv](doc/inf1.csv) 
edits Q4ACW9, [inf2.csv](doc/inf2.csv) is inline with UCSC with respect to P12034 and P30203.

A version by Jimmy is [olink.inf.panel.annot.tsv](doc/olink.inf.panel.annot.tsv) from [olink.annotation.R](doc/olink.annotation.R).

BDNF has recently been removed from the assay and replaced with CD8A, https://www.olink.com/bdnf-info/, and there are also changes on TNF and IFN.gamma, https://www.olink.com/inflammation-upgrade/.

## Files, URLs and downloading

### Files

File     | Description
---------|--------------------------------------------------------
[doc/](doc) | Some documents and auxiliary files
[cardio/](cardio) | Work on CEU's cardio
[tryggve/](tryggve) | Analysis programs on TRYGGVE
[SCALLOP_INF1_analysis_plan.md](SCALLOP_INF1_analysis_plan.md) | Analysis plan
[SCALLOP_INF1_analysis_plan.docx](SCALLOP_INF1_analysis_plan.docx) |
[tryggve.md](tryggve.md) | TRYGGVE notes

### URLs

SCALLOP [consortium](https://www.olink.com/scallop/), [GitHub repository](https://github.com/lassefolkersen/scallop), [securecloud](https://secureremote.dtu.dk/vpn/index.html).

SomaLogic plasma protein GWAS summary statistics, http://www.phpc.cam.ac.uk/ceu/proteins.

[turbo.sh](cardio/turbo.sh) shows tweak of [turboman](https://github.com/bpprins/turboman) and [turboqq](https://github.com/bpprins/turboqq); see also aild() in tryggve/analysis.sh.

[Olink publications](https://www.olink.com/data-you-can-trust/publications/).

### downloading

``` {.bash}
git clone https://github.com/jinghua/INF
```
