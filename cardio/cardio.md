# INTERVAL reference

/scratch/jhz22/data/INTERVAL

# INTERVAL sumstats

location: /scratch/jhz22/INF/sumstats/INTERVAL

# Manhattan plots by cohort

location: /scratch/jhz22/INF/plots

# Clumping results

location: /scratch/jhz22/INF/clumping

Files: study.reference.LD-cutoff.clumped_reults, e.g., INTERVAL.UK10K+1KG.r2-0.1.clumped, while AILD is approximately independent linkage disequilibrium as described in FM-pipeline.

# Conditional/joint analysis

location: /scratch/jhz22/INF/cojo

where aild-indel/ and aild-snp/ contain results for all variants and SNPs, respectively.

# PhenoScanner outputs

Location, /scratch/jhz22/INF/ps

which contains protein-specific and whole (INF1_* and INF1r_*, see below) results.
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

# Meta-analysis sumstats + Manhattan/Q-Q/LocusZoom plots

location: /scratch/jhz22/INF/METAL

---

# GCTA (1-3) and Clumping (4-7) results

The directories are relative to /scratch/jhz22/INF,

 No | option | [cis, trans, total] | Comments and locations (cojo|clumping|doc)
---------- | ------- | ------------------- | -----------------------------------------------------------
1 | LD blocks | [228, 182, 410] | only SNPs, cojo/aild-snp/INF1.jma.\*, also doc/INF1.paper.xlsx
~ | ~ | [254,  191, 445] | SNPs+indels, cojo/aild-indel/INF1.jma.\*
2 | default  | [234, 173, 407] | --cojo-collinear 0.9 --cojo-wind 10000, doc/SCALLOP_INF1-260419.xlsx
3 | small R2 & window | [189, 186, 375] | --cojo-collinear 0.1 --cojo-wind 500, doc/SCALLOP_INF1-260419.xlsx
4 | INTERVAL LD panel | [659, 275, 934] | --clump-r2 0.1 --clump-kb 500, doc/SCALLOP_INF1-120419.xlsx
5 | 1000G LD panel | [405, 229, 634] | --clump-r2 0.1 --clump-kb 500, clumping/INF1.1KG.r2-0.1.clumped*
6 | INTERVAL data  | [424, 188, 612] | --clump-r2 0.1 --clump-kb 500, doc/SCALLOP_INF1-120419.xlsx
7 | 1000G LD panel | [402, 226, 628] | --clump-r2 0.1 --clump-kb 1000, on tryggve

Factors on number of signals,

* indels lead to more signals (1, 2).
* increase with value of --cojo-collinear and --cojo-wind (2, 3), yet moderate changes in LD window have less impact than panel (5, 7).
* PLINK --clump gives more signals than GCTA --cojo (3, 4).
* INTERVAL as LD reference leads to more signals than 1000Genomes (4,  5).
* sample size. INF1 has more than INTERVAL alone (5, 6).
