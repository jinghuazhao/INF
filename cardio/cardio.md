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

# /scratch/jhz22/INF/ps

which contains protein-specific and whole (INF1_* and INF1r_*, see below) results.

```bash
# those INF1_*
# default all results, no LD, p=1e-5
# phenoscanner -c All -l No -p 0.00001 -i INF1.ps -o INF1

# revised
# GWAS only, LD r=0.6, p=1e-7
phenoscanner -c All -l 1000G -p 0.0000001 -i INF1.ps -o INF1r
```

# Meta-analysis sumstats + Manhattan/Q-Q/LocusZoom plots

location: /scratch/jhz22/INF/METAL
