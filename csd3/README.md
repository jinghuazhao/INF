# scripts on CSD3

`rsid.sb` generates two versions of SNPID-rsid mapping in `.rsid` files to be used by `finemap.R`, `gcta.R` and `jam.R`.

* `finemap.sb` and `slct.sb` use the unpruned version.
* `fm.sb` and `INTERVAL-fm.sb` use a version which only contains pruned variants to comproise `JAM`. The `.z` file is also appropriate for both `finemap` and `JAM`.
```bash
# ldstore-1

export pr=IL.6-chr1:154426970_A_C

ldstore --bcor ${pr} --bgen ${pr}.bgen
ldstore --bcor ${pr} --merge 1
ldstore --bcor ${pr} --matrix ${pr}.ld

# ldstore-2 (and finemap-1.4)

ldstore-2 --in-files ${pr}.master2 --write-bcor --write-bdose

finemap-1.4 --sss --in-files ${pr}.master2 --n-causal-snps 10
```
where ${pr}.master2 contains two lines,
```
z;bgen;bgi;bcor;bdose;snp;config;cred;log;n_samples
IL.6-chr1:154426970_A_C.z;IL.6-chr1:154426970_A_C.bgen;IL.6-chr1:154426970_A_C.bgi;IL.6-chr1:154426970_A_C.bcor;IL.6-chr1:154426970_A_C.bdose;IL.6-chr1:154426970_A_C.snp;IL.6-chr1:154426970_A_C.config;IL.6-chr1:154426970_A_C.cred;IL.6-chr1:154426970_A_C.log;4994
```
