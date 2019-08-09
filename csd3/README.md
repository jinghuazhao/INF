# scripts on CSD3

`rsid.sb` generates two versions of SNPID-rsid mapping in `.rsid` files to be used by `finemap.R`, `GCTA.R` and `JAM.R`.

* `finemap.sb` and `slct.sb` use the unpruned version.
* `fm.sb` and `INTERVAL-fm.sb` use a version which only contains pruned variants to comproise `JAM`. The `.z` file is also appropriate for both finemap and `JAM`.

