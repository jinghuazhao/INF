# GWAS sumstats and plots

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

* cojo/, GCTA --cojo-slct results for all (INF1*) and individual proteins
  * ps/, PhenoScanner v2 results. coded as in /scratch/jhz22/INF/cardio/ps.sh
* clump/, PLINK --clump results for all (INF1*) and individual proteins

## SNPs only

location, /scratch/jhz22/INF/snps

* cojo/, GCTA --cojo-slct results for all (INF1*) and individual proteins
  * ps/, PhenoScanner v2 results.
  * lz/, LocusZoom plots
* clump/, PLINK --clump results for all (INF1*) and individual proteins
* clump0/ PLINK --clump-r2 0 for independent signals

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
