# Cardio

Cardio is the Linux clusters once used by the Cardiovascular Epidemiology Unit at University of Camrbridge.

## INTERVAL reference

/scratch/jhz22/data/INTERVAL

## INTERVAL sumstats

location: /scratch/jhz22/INF/sumstats/INTERVAL

## Manhattan plots by cohort

location: /scratch/jhz22/INF/plots

## Meta-analysis sumstats + Manhattan/Q-Q/LocusZoom plots

location: /scratch/jhz22/INF/METAL

---

# Results based on aproximately independent LD blocks (AILD)

*See https://bitbucket.org/nygcresearch/ldetect-data*

## SNPs+indels

location, /scratch/jhz22/INF/aild

* cojo/, **GCTA** --cojo-slct results for all (INF1*) and individual proteins
  * ps/, **PhenoScanner** v2 results. coded as in /scratch/jhz22/INF/cardio/ps.sh
* clump/, **PLINK** --clump results for all (INF1*) and individual proteins

## SNPs only

location, /scratch/jhz22/INF/snps

* cojo/, **GCTA** --cojo-slct results for all (INF1*) and individual proteins
  * ps/, **PhenoScanner** v2 results.
  * lz/, **LocusZoom** plots
* clump/, **PLINK** --clump results for all (INF1*) and individual proteins
* clump0/ **PLINK** --clump-r2 0 for independent signals

---

# Genomewide results based on slide windows

## Clumping results

location: /scratch/jhz22/INF/clumping

Files: study.reference.LD-cutoff.clumped_reults, e.g., INTERVAL.UK10K+1KG.r2-0.1.clumped

## Conditional/joint analysis

location: /scratch/jhz22/INF/cojo

## Conditional/joint analysis with regions in high LD

location: /scratch/jhz22/INF/aild-ld

* aild-indel/, all variants
  * ps, PhenoScanner v2 and v1.1/ 
* aild-snp/,  SNPs only

Results from PhenoScanner which contains protein-specific and whole (INF1_* and INF1r_*, see below) results.
```bash
# those INF1_*
# All results, no LD, p=1e-5
# phenoscanner -c All -l No -p 0.00001 -i INF1.ps -o INF1

# revised
# All results, no LD, p=1e-7
phenoscanner -c All -l No -p 0.0000001 -i INF1.ps -o INF1r
```
The results on GWAS with LD r2=0.6 in INF1r_* are available from R as follows,
```bash
R --no-save -q <<END
  library(MendelianRandomization)
# SNPID not possible
  INF1 <- read.table("INF1.ps",as.is=TRUE)
# batches in 100 only
  INF1 <- with(read.csv("INF1_PhenoScanner_SNP_Info.csv"),rsID)
  r1 <- phenoscanner(snpquery=  INF1[1:100], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  r2 <- phenoscanner(snpquery=INF1[101:200], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  r3 <- phenoscanner(snpquery=INF1[201:300], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  r4 <- phenoscanner(snpquery=INF1[301:376], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  save(r1,r2,r3,r4,file="INF1r.rda",version=2)
END
```

## Approximately independent LD blocks:

1. Set up 1703 autosomal regions as defined in [EURLD.bed](../tryggve/EURLD.bed).
2. Extract variants outside the [regions in high LD](../tryggve/high-LD-regions-hg19.txt) to [1672 regions](../tryggve/EURLD-no-high-LD-regions-hg19.bed) by [EURLD.sh](../tryggve/EURLD.sh).
3. Overlap regions and GWAS sumstats:
   * Tag GWAS sumstats with regions through [aild-rma.sb](../cardio/aild-rma.sb).
   * Pair protein-region which contains genomewide significant signals by [aild-list.sb](../cardio/aild-list.sb).
   * Independently, list variants by region in the reference panel by [aild-snplist.sb](../cardio/aild-snplist.sb).
4. clump via [aild-clump.sb](../cardio/aild-clump.sb).
5. cojo via [aild-cojo.sb](../cardio/aild-cojo.sb).
6. Downstream analyses with PhenoScanner (preferably v2) as in [ps.sh](../cardio/ps.sh) and forest plots on TRYGGVE (with study-specific sumstats).

The regions are predefined. As shown in [EURLD.tsv](../tryggve/EURLD.tsv) by [EURLD.R](../tryggve/EURLD.R), the LD patterns across the genome are more variable than the norm in a typical genomewide association analysis therefore slide windows such as 250kb (36), 500kb (300), or even 10Mb (1071), seeing that the sentinel variant may not necessarily lie right in the middle of a window. The number of signals in our case were close to GCTA but overestimated (53 by PLINK) as in the following table. For instance, it is often seen from the PLINK --clump-range outputs that sliding windows can give results in two neighbouring LD blocks.

Note that pairing regions of interest would reduce the burden of genomewide analysis, and also that region-specific reference will not affect results from steps 4 and 5 regarding use of variants from GWAS sumstats.

Steps 4 and 5 both use `INF1.aild`, which contains all the protein-region pairs. The results are classified as in [analysis.sh](../tryggve/analysis.sh). In particular, for step 5 this is done with [aild.sh](../cardio/aild.sh).

## Near-independent signals from cojo and clumping

**Run** | **Option** | **cis** | **trans** | **total** | **Comments/location<sup>\+</sup>**
-----------|----------|--------------|-----------|------------|----------------------------------------------
**GCTA** |
1 | LD blocks | 210 | 147 | 357 | only SNPs, cojo/aild-snp/INF1.jma.*, also doc/INF1.paper.xlsx
\+ indels | LD blocks | 220 | 155 | 375 | SNPs+indels, cojo/aild-indel/INF1.jma.*
2 | default | 234 | 173 | 407 | --cojo-collinear 0.9 --cojo-wind 10000, doc/SCALLOP_INF1-260419.xlsx
3 | small R2 & window | 189 | 186 | 375 | --cojo-collinear 0.1 --cojo-wind 500, doc/SCALLOP_INF1-260419.xlsx
**PLINK** |
4 | LD blocks | 594 | 252 | 846 | only SNPs, clumping/aild-snp/INF1.jma.*, also doc/INF1.paper.xlsx
\+ indels | LD blocks | 621 | 258 | 879 | SNPs+indels, clumping/aild-indel/INF1.jma.*
5 | INTERVAL LD panel | 657 | 275 | 932 | --clump-r2 0.1 --clump-kb 500, doc/SCALLOP_INF1-120419.xlsx
6 | 1000G LD panel | 405 | 229 | 634 | --clump-r2 0.1 --clump-kb 500, clumping/INF1.1KG.r2-0.1.clumped.*
7 | INTERVAL data | 424 | 188 | 612 | --clump-r2 0.1 --clump-kb 500, doc/SCALLOP_INF1-120419.xlsx
8 | 1000G LD panel | 402 | 226 | 628 | --clump-r2 0.1 --clump-kb 1000, on tryggve

<sup>\+</sup>The directories are relative to /scratch/jhz22/INF, i.e., doc/, cojo/ and clumping/, Results in 2 and 3 include regions in high LD excluded in other analyses.

A few observations can be made,

* indels lead to more signals in cojo (1) and clumping (4) analyses.
* **default GCTA --cojo-collinear and --cojo-wind parameters did quite well in numbers, esp. taking ~30 away regions in LD** (1, 2).
* Although it looks close, GCTA --cojo-collinear 0.1 produces considerably less signals compared to --cojo-collinear 0.9 (3).
* the number of signals increase with the values of GCTA parameters (2, 3), yet moderate changes in LD window have less impact than the reference panel (5, 8).
* PLINK --clump gives more signals than GCTA --cojo (1, 4 and 2, 5).
* **Specification of sliding LD windows disregarding AILD patterns in clumping gives 53 additional signals** (4, 5).
* Thanks to the larger sample size and perhaps greater variant number, INTERVAL as LD reference leads to more signals than 1000Genomes (5, 6).
* Summary statistics from larger sample size gives more signals (5, 7).
* Unpruned results are likely to give more cis signals but this is subject to scrutiny perhaps on individual cases.

It can be concluded that it is desirable to employ approximately independent LD blocks for both GCTA (1) and PLINK (4), and also that reference such as UK10K+1KG would be desirable with respect to both sample size and variant number.

*Date last changed:* **23/8/2021**

