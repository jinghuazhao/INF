# scripts on CSD3

`rsid.sb` generates two versions of SNPID-rsid mapping in `.rsid` files to be used by `finemap.R`, `gcta.R` and `jam.R`.

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

# ukb

The utility `ukb.sh` extracts from ukb_imp_chr[x-xx]_v3.bgen as in
```
/DO-NOT-MODIFY-SCRATCH/uk_biobank/500k/imputed_v3
/DO-NOT-MODIFY-SCRATCH/curated_genetic_data/uk_biobank/reference_files/full_release
```
on Cardio into region-specific data in `bgen` format (ukb/bgen) according to work/INF1.merge.

It involves a utility called `ukb.sb` to generate `binary_ped` with SNPIDs (bgen/*map)
* by qctool -- it also creates `bgen` with SNPIDs (ukb/nodup) to avoid hard-called genotypes.
* by PLINK -- it has a `_snpid` tag.

In both cases the duplicates (bgen/*rmdup.list) are excluded.

# SNPID-rsid mappings

This is furnished with `snpid_rsid.sb`, whose results will be attached to GCTA/finemap results.
