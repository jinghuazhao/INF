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
# default all results, no LD, p=1e-5
# phenoscanner -c All -l No -p 0.00001 -i INF1.ps -o INF1

# revised
# GWAS only, LD r=0.6, p=1e-7
phenoscanner -c GWAS -l 1000G -p 0.0000001 -i INF1.ps -o INF1r
```
The results in INF1r_* are also available from R as follows,
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
