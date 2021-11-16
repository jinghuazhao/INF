# CSD3

CSD3 is the HPC Linux clusters used at the University of Cambridge, with which this document is associated.

**NOTE**: [merge.sh](merge.sh) and [INTERVAL-fm.ini](INTERVAL-fm.ini) now work directly with the original INTERVAL bgen files. The alphabetical list is as follows.

## 3d scatter plot

Furnished with [js.R](js.R), the JSON output could be used as data in the Supplementary Figure 1 of Sun et al. (2018).

## Annotation

[annotate.sh](annotate.sh) involves **ANNOVAR**, **PolyPhen 2**, **VEP** and some R packages.

## Known variants

[lookup.sh](lookup.sh), which uses PhenoScanner output from [pqtl.R](pqtl.R), can be used to check if pQTLs were known.

File [cvd1.sh](cvd1.sh) gives sentinels overlapped with CVD I as well as LD information corresponding to each protein.

## LD reference panels

The .sh versions below extract data from INTERVAL and ukb, each calling a .sb to generate `binary_ped` with SNPIDs (bgen/*map)
* by **qctool** -- it also creates `bgen` with SNPIDs (ukb/nodup) to avoid hard-called genotypes.
* by **PLINK** -- it has a `_snpid` tag.

In both cases duplicates (bgen/*rmdup.list) are excluded.

1. The utility [INTERVAL.sh](INTERVAL.sh) extracts data from  
```
/DO-NOT-MODIFY-SCRATCH/bp406/data_sets/interval_subset_olink/genotype_files/unrelated_4994_pihat_0.1875_autosomal_typed_only
/DO-NOT-MODIFY-SCRATCH/bp406/data_sets/interval_subset_olink/genotype_files/unrelated_4994_pihat_0.1875_autosomal_imputed_info_0.4_phwe_1e-4_filtered/per_chr
```
on Cardio into region-specific data in `bgen` format according to work/INF1.merge.

2. The utility [ukb.sh](ukb.sh) extracts data from ukb_imp_chr[x-xx]_v3.bgen as in
```
/DO-NOT-MODIFY-SCRATCH/uk_biobank/500k/imputed_v3
/DO-NOT-MODIFY-SCRATCH/curated_genetic_data/uk_biobank/reference_files/full_release
```
on Cardio into region-specific data in `bgen` format (ukb/bgen) according to work/INF1.merge.

## SNPID-rsid mappings, conditional analysis and finemapping

This is furnished with [snpid_rsid.sb](snpid_rsid.sb), whose results will be attached to **GCTA**/**finemap** results via `finemap.R`, `slct.R` and `jam.R`.
After that, it is ready to use script `slct.sb` followed by `slct.sh`. Optionally, the results are fed into `finemap.sb` via `--n-causal-snps`.

* [finemap.sb](finemap.sb) uses the unpruned version.
* [fm.sb](fm.sb) and [INTERVAL-fm.sb](INTERVAL-fm.sb) use pruned variants to compromise `JAM`. The `.z` file is also appropriate for both `finemap` and `JAM`.
* [INTERVAL-fm.sh](INTERVAL-fm.sb) and [INTERVAL-fm.ini](INTERVAL-fm.ini) works on INTERVAL data.

```bash
# ldstore v1.1

export pr=IL.6-chr1:154426970_A_C

ldstore --bcor ${pr} --bgen ${pr}.bgen
ldstore --bcor ${pr} --merge 1
ldstore --bcor ${pr} --matrix ${pr}.ld

# ldstore 2 (and finemap-1.4)

ldstore --in-files ${pr}.master2 --write-bcor --write-bdose

finemap --sss --in-files ${pr}.master2 --n-causal-snps 10
```

where ${pr}.master2 contains two lines,
```
z;bgen;bgi;bcor;bdose;snp;config;cred;log;n_samples
IL.6-chr1:154426970_A_C.z;IL.6-chr1:154426970_A_C.bgen;IL.6-chr1:154426970_A_C.bgi;IL.6-chr1:154426970_A_C.bcor;IL.6-chr1:154426970_A_C.bdose;IL.6-chr1:154426970_A_C.snp;IL.6-chr1:154426970_A_C.config;IL.6-chr1:154426970_A_C.cred;IL.6-chr1:154426970_A_C.log;4994
```
The environmental variable `LDREF` provides an option to use either INTERVAL or UKB data.

## Low abundance check

This is down with [qc.sh](qc.sh).

## PGS

These are [prsice.sh](prsice.sh) and [pgs.sh](pgs.sh) for the CRP example.

## PhenoScanner

This is furnished with [ps.sh](ps.sh).

## QQ and Manhattan plots

The version for meta-analysis was part of `analysis.sh` at tryggve/. [qqman.sh](qqman.sh) calls turboqq and turboman by Bram Prins.

## Sentinel identification

This is now furnished with [merge.sh](merge.sh).

Clumping by fixed distance is superseded with its failure to handle long LD regions.

There are scripts for heritability analysis and proportion of variance explained.

## Trans-pQTL hotspots and proteins as polygenic targets

This is illustrated with circos plots and by default this works on `work/INF1.merge.cis.vs.trans` and requires [glist-hg19](glist-hg19).
Because circos plots are gene-centric, in both cases, the protein-coding gene is handled in mind.

Try
```bash
./hotspot.sh chr1:159175354_A_G
./hotspot.sh chr12:111884608_C_T
```
giving results linking *ACKR1* and *SH2B3*, respectively, while
```bash
./polygene.sh IL12B
./polygene.sh TNFSF10
```
linking IL.12B and TRAIL, respectively.

To generate all possible plots, wo do
```bash
for h in $(cut -f6 ${INF}/work/INF1.merge | sed '1d' | sort -k1,1 | uniq); do echo $h; hotspot.sh $h; done
for g in $(cat ${INF}/work/INF1.merge.gene); do echo $g; polygene.sh $g; done
```

A brick-and-tile operation is illustrated with the following code,

```bash
# long format

bedtools intersect -a <(cut -d, -f3,4,5,7 ${INF}/work/INF1.merge.cis.vs.trans | awk -vFS="," -vOFS="\t" 'NR>1{print "chr" $1,$2-1, $2, $3, $4}') \
                   -b <(sort -k1,1n -k2,2n ${INF}/csd3/glist-hg19 | grep -v -w -e X -e Y -e XY | awk '{$1="chr" $1;print}' | tr ' ' '\t') \
                   -wa -wb -loj > sentinel-glist.long

# compact format
Rscript -e '
  options(width=1000)
  d <- read.table("sentinel-glist.long",as.is=TRUE)
  a <- aggregate(d,by=list(with(d,V1),with(d,V2),with(d,V3)),FUN="c")
  sink("sentinel-glist.wide")
  a
  sink()
  for (v in c(1:4,6:9))
  {
      s <- lapply(a[,3+v],unique,"[[")
      a[,3+v] <- unlist(lapply(s,paste,collapse=","))
  }
  sink("sentinel-glist.short")
  a[paste0("V",1:9)]
  sink()
'
```

*Date laste changed:* **16/11/2021**
