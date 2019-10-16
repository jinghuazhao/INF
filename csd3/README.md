# Scripts on CSD3

## Sentinel identification

This is furnished with `merge.sh`.

## LD reference panels

The .sh versions below extract data from INTERVAL and ukb, each calling a .sb to generate `binary_ped` with SNPIDs (bgen/*map)
* by qctool -- it also creates `bgen` with SNPIDs (ukb/nodup) to avoid hard-called genotypes.
* by PLINK -- it has a `_snpid` tag.

In both cases the duplicates (bgen/*rmdup.list) are excluded.

### INTERVAL.sh (INTERVAL.sb)

The utility `INTERVAL.sh` extracts data from  
```
/DO-NOT-MODIFY-SCRATCH/bp406/data_sets/interval_subset_olink/genotype_files/unrelated_4994_pihat_0.1875_autosomal_typed_only
/DO-NOT-MODIFY-SCRATCH/bp406/data_sets/interval_subset_olink/genotype_files/unrelated_4994_pihat_0.1875_autosomal_imputed_info_0.4_phwe_1e-4_filtered/per_chr
```
on Cardio into region-specific data in `bgen` format according to work/INF1.merge.

### ukb (ukb.sh)

The utility `ukb.sh` extracts data from ukb_imp_chr[x-xx]_v3.bgen as in
```
/DO-NOT-MODIFY-SCRATCH/uk_biobank/500k/imputed_v3
/DO-NOT-MODIFY-SCRATCH/curated_genetic_data/uk_biobank/reference_files/full_release
```
on Cardio into region-specific data in `bgen` format (ukb/bgen) according to work/INF1.merge.

## SNPID-rsid mappings

This is furnished with `snpid_rsid.sb`, whose results will be attached to GCTA/finemap results via `finemap.R`, `gcta.R` and `jam.R`.

* `finemap.sb` and `slct.sb` use the unpruned version.
* `fm.sb` and `INTERVAL-fm.sb` use a version which only contains pruned variants to comproise `JAM`. The `.z` file is also appropriate for both `finemap` and `JAM`.
```bash
# ldstore v1.1

export pr=IL.6-chr1:154426970_A_C

ldstore --bcor ${pr} --bgen ${pr}.bgen
ldstore --bcor ${pr} --merge 1
ldstore --bcor ${pr} --matrix ${pr}.ld

# ldstore 2.0b (and finemap-1.4)

ldstore-2 --in-files ${pr}.master2 --write-bcor --write-bdose

finemap-1.4 --sss --in-files ${pr}.master2 --n-causal-snps 10
```
where ${pr}.master2 contains two lines,
```
z;bgen;bgi;bcor;bdose;snp;config;cred;log;n_samples
IL.6-chr1:154426970_A_C.z;IL.6-chr1:154426970_A_C.bgen;IL.6-chr1:154426970_A_C.bgi;IL.6-chr1:154426970_A_C.bcor;IL.6-chr1:154426970_A_C.bdose;IL.6-chr1:154426970_A_C.snp;IL.6-chr1:154426970_A_C.config;IL.6-chr1:154426970_A_C.cred;IL.6-chr1:154426970_A_C.log;4994
```

## Conditional analysis and finemapping

After snpid-rsid.sb is called, it is ready to use script `slct.sb` followed by `slct.sh`.

Optionally, the results are fed into `finemap.sb` via `--n-causal-snps`.

## Clumping by fixed distance

This is available as `sentinels_nold.sh` but is superseded with its failure to handle long LD regions.

There are scripts for heritability analysis and proportion of variance explained.

## Protein with polygenic effects and trans-pQTL hotspots

This is illustrated with circos plots and by default this works on `work/INF1.merge.cis.vs.trans` and requires [glist-hg19](glist-hg19).
Because circos plots are gene-centric, in both cases, the protein-coding gene is handled in mind.

Try
```bash
sh hotspot.sh chr1:159175354_A_G
sh hotspot.sh chr12:111884608_C_T
```
giving *ACKR1* and *SH2B3*, respectively, and
```bash
sh polygene.sh TNFSF10
```
linking TRAIL.

## PhenoScanner

This is furnished with `ps.sh`.
