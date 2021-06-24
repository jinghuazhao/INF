# Tryggve

The repository covers various aspects such as software and data.

## Modules

Instances are exemplified as follows,

``` {.bash}
module load anaconda2/4.4.0
module load bedtools/2.28.0
module load bgen/20180807
module load perl/5.24.0 annovar/2019oct24
module load bcftools/1.9
module load emacs/26.1
module load gcta/1.91.0beta
module load intel/redist/2019 intel/perflibs/64/2019 gcc/5.4.0 lapack/3.8.0 R/3.5.3-ICC-MKL rstudio/1.1.453
module load libreoffice/6.0.5.2
module load locuszoom/1.4
module load metal/20180828
module load pandoc/2.1
module load parallel/20190122
module load plink2/1.90beta5.4
module load vcftools/0.1.16
module load xpdf/3.04
export threads=1
```

e.g., 

```bash
pandoc tryggve.md -o tryggve.docx
soffice tryggve.docx
xterm -fa arial -fs 12 -bg black -fg white
```

## Parallel computing

[METAL.qsub](METAL.qsub) was intended to run batch jobs under [Terascale Open-source
Resource and QUEue Manager (TORQUE)](https://en.wikipedia.org/wiki/TORQUE) controls batch jobs and
distributed compute nodes, enabling integration with [Moab cluster
suite](https://en.wikipedia.org/wiki/Moab_Cluster_Suite) and extending
[Portable Batch System (PBS)](https://en.wikipedia.org/wiki/Portable_Batch_System) with respect
to extend scalability, fault tolerance, and functionality. Information on
job scheduling on Computerome is availble from https://www.computerome.dk/display/CW/Batch+System,

However, the SCALLOP securitycloud is a 1-node cloud and we resort to
parallel instead, e.g,

```bash
module load metal/20110325 parallel/20190122
ls METAL/*.run | parallel --dry-run --env HOME -j8 -C' ' 'metal $HOME/INF/{}'
```
NB METAL add -1 to the filenames.

## SOFTWARE ISSUES

### Bash functions

There might be problems with R/3.5.3-ICC-MKL so functions are defined in `.bashrc` which can be invoked at the start of a session,
then one invokes R.3.5.3 when necessary. 
```bash
# This could be very slow
function R.3.5.3()
{
  export R_LIBS=/data/$USER/R:/services/tools/R/3.5.3-ICC-MKL/lib64/R/library
  module load intel/redist/2019 intel/perflibs/64/2019 gcc/5.4.0 lapack/3.8.0 R/3.5.3-ICC-MKL lapack/3.8.0
  source /data/jinhua/parallel-20190222/bin/env_parallel.bash
  alias R='/services/tools/R/3.5.3-ICC-MKL/bin/R -q $@'
}
# An alternative
function R.3.3.1()
{
  module load intel/redist/2019 intel/perflibs/64/2019 gcc/5.4.0 R/3.3.1-ICC-MKL
}
```
We also use R 3.3.x below.

## Code extraction

This is possible with codedown, as noted [here](https://github.com/jinghuazhao/physalia/blob/master/NOTES.md).

### LocusZoom 1.4

To use LocusZoom 1.4, one only needs to start with
```bash
module load gcc/5.4.0
module load R/3.2.5
module load anaconda2/4.4.0
module load locuszoom/1.4
```

### finemap/ldstore

finemap 1.3.1 and ldstore 1.1 are available from /data/jinhua and can be enabled as follows,
```bash
ln -fs /data/jinhua/ldstore_v1.1_x86_64/ldstore $HOME/bin/ldstore_v1.1
ln -fs /data/jinhua/finemap_v1.3.1_x86_64/finemap_v1.3.1_x86_64 $HOME/bin/finemap_v1.3.1
````
finemap 1.4 and ldstore 2.0b are available from /data/jinhua/finemap-1.4.

### qctool

TRYGGVE now fixed issue with qctool/2.0.1 for lack of lapack shared libraries as in /data/jinhua/lapack-3.8.0/ and its installation described on GitHub repository, https://github.com/jinghuazhao/Computational-Statistics.

It is then possible to issue `module load gcc/5.4.0 lapack/3.8.0 qctool/2.0.1`.

## NEW SOFTWARE

There are a number of software updates/additions which are worthy of note.

### GCTA

A more recenve version is available from /data/jinhua/gcta_1.91.7beta/.
This version can handle chi-squared statistics instead of p values in the joint/conditional (COJO) analysis; it also allows for `--grm file --pca --out file`, i.e., same file root.

### GNU parallel

The latest version parallel-20190222 has new features and use them directly on TRYGGVE without invoking modules as follows,
```bash
export src=/data/jinhua/parallel-20190222/bin
for i in $(ls $src); do ln -fs $src/$i $HOME/bin/$i; done
export MANPATH=/data/jinhua/parallel-20190222/share/man:$MANPATH
```
The last line enables `man parallel` and `info parallel`.

### ImageMagick

This is version 7.0.8-22, made available due to the inability to use the imagemagick/7.0.8-16 module, e.g.,
```bash
export PATH=/data/jinhua/ImageMagick-7.0.8-22/bin:$PATH
convert OPG.lz-1.png -resize 130% OPG.lz-3.png
convert \( OPG.qq.png -append OPG.manhattan.png -append OPG.lz-3.png -append \) +append OPG-qml.png
```
used to generate the figure in the front page. Another very useful utility is its `display`.

### METAL

The version as in /data/jinhua/METAL-2018-08-28 contains modification which allows for CUSTOMVARIABLE to use integer position rather than scientific format as in [software-notes](https://github.com/jinghuazhao/software-notes/blob/master/AA.md).

An extension was made with the Direction column so that if a variant with '+' and '-' effects becomes 'p' and 'n', respectively with P <= 0.05.

To avoid loading the default /usr/bin/metal, one can add
```bash
ln -sf /data/jinhua/METAL-2018-08-28/metal $HOME/bin/metal
export PATH=$HOME/bin:/data/jinhua/ImageMagick-7.0.8-22/bin:$PATH
```
into $HOME/.bashrc.

### OpenVPN

In case you see message "All TAP-Windows aapters on this system are currently in use", read [this post](https://vpn.ac/knowledgebase/64/OpenVPN-Error-All-TAP-Windows-adapters-on-this-system-are-currently-in-use.html), esp. install [tap-windows-8.21.2](https://swupdate.openvpn.org/community/releases/tap-windows-9.21.2.exe).

The log is usually located at "%USERPROFILE%/OpenVPN/log".

### R

**EasyQC/EasyStrata**

The version is 18.1 rather than 9.2 and 8.6 currently online, https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/software/.

**gap**

First, run R.3.5.3 as defined above.

This version contains functions cis.vs.trans.classification, gc.lambda, invnormal, and here is the way to go
```bash
module load intel/redist/2019 intel/perflibs/64/2019 gcc/5.4.0 R/3.5.3-ICC-MKL
tar xvfz gap_1.2.2.tar.gz
cd gap/src
gcc -I/services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/include -L/services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/lib/release -L/services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/lib -Xlinker --enable-new-dtags -Xlinker -rpath -Xlinker /services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/lib/release -Xlinker -rpath -Xlinker /services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/lib -Xlinker -rpath -Xlinker /opt/intel/mpi-rt/2017.0.0/intel64/lib/release -Xlinker -rpath -Xlinker /opt/intel/mpi-rt/2017.0.0/intel64/lib -lmpifort -lmpi -ldl -lrt -lpthread -L/services/tools/intel/perflibs/2019//compilers_and_libraries_2019.0.117/linux/mpi/intel64/libfabric/lib -fPIC -c *.c *.f
gcc -shared -L/services/tools/R/3.5.3-ICC-MKL/lib64/R/lib -L/usr/local/lib64 -o gap.so 2k.o 2ld.o cline.o gcontrol_c.o gcx.o gif_c.o hap_c.o hwe.hardy.o kin.morgan.o makeped_c.o mia.o muvar.o package_native_routine_registration_skeleton.o pfc.o pfc.sim.o pgc_c.o whscore.o -L/usr/lib/gcc/x86_64-redhat-linux/4.8.2 -lgfortran -lm -lquadmath -L/services/tools/R/3.5.3-ICC-MKL/lib64/R/lib -lR
cd -
R CMD INSTALL gap -l /data/jinhua/R
```
The rather imtimidating compiler flags are derived from `mpiicc -show` as with initial pass requested for libfabric.so.1, and the idea is to get around check for icc which does not exist on the system. To faciliate compiling, [gap.sh](gap.sh) contains these lines for use.

To take advantage of the `circos.mhtplot` function, package R/gap.datasets is also made available. The following is an example for cross-check,
```bash
R --no-save -q <<END
  gz <- "sumstats/STABILITY/STABILITY.IFN.gamma.gz"
  STABILITY <-  read.table(gz,as.is=TRUE,header=TRUE,sep="\t")
  summary(STABILITY)
  library("GenABEL")
  estlambda(with(STABILITY,PVAL), method="median")
END
```
Somehow R/gap thus compiled was very slow; so attempt was also made with R 3.3.1 via [gap.3.3.1.sh](gap.3.3.1.sh) which requires registration of routine `family` (defined in `pfc.f`) altered to `family_` and the change is made as default.

The `fp()` function in `analysis.sh` requires `R.3.3.1` to be called.

**QCGWAS**

The version contains fix to the use of HapMap reference as in software-notes above. The HapMap data as with code from the packages's quick guide is /data/jinhua/data/QCGWAS.

## A SUMMARY OF FILES

File specification | Descrription
-------------------|---------------------------------------------------------------
doc/               | Oringal documents
 -- [KORA.prot.preproc.belowlod.v2.R*](KORA.prot.preproc.belowlod.v2.R) | R code for data preprocessing
 -- kora.normalised.prot.txt* | sex, age, normalised proteins
 -- [KORA.pc.below.llod.pdf*](KORA.pc.below.llod.pdf) | llod check
METAL/             | METAL/output scripts by protein
sumstats/          | File lists and study directories
tryggve/           | Auxiliary files
 -- [list.sh](list.sh)         | Generation of file list and directory
 -- [format.sh](format.sh)     | Code for format GWAS summary statistics
 -- [lz14.sh](lz14.sh)         | Code to extract LocusZoom 1.4 databases for [analysis.sh](analysis.sh)
 -- [analysis.sh](analysis.sh) | Bash code for analysis calling [analysis.ini](analysis.ini)
 -- [metal.sh](metal.sh)       | Generation/execution of METAL scripts
 -- [METAL.qsub](METAL.qsub)   | TORQUE qsub script for METAL
 -- [INTERVAL.sh](INTERVAL.sh) | INTERVAL analysis
 -- [KORA.sh](KORA.sh) | Bash/R scripts to hand KORA data
 -- [KORA.txt](KORA.txt) | obsolete version to use simulated data
 -- [KORA.R](KORA.R) | R code to simulate phenotypes
 -- [EURLD.bed](EURLD.bed) | approximately independent LD blocks for Europeans.
 -- [QCGWAS.sh](QCGWAS.sh) | QCGWAS for specific proteins, calling [QCGWAS.R](QCGWAS.R)
 -- [snpstats.sh](snpstats.sh) | sumstats/qctool -snp-stats summary, see also [qctool.txt](qctool.txt)
tryggve.md         | This document

\* from Jimmy

In total, 92 proteins are expected as given in [olink.prot.list.txt](../doc/olink.prot.list.txt).

-   **BioFinder**. 91 (no BDNF) proteins. sumstats file named after genes and converted to protein names.
-   **NSPHS**. 91 (no BDNF) proteins, originally a tar.gz file is unpacked into \$HOME/INF/work leading to 10 proteins
-   **EGCUT**. 91 (no BDNF) proteins, orginally only 18 proteins though stratified by chromsomes
-   **INTERVAL**. 92 proteins. raw SNPTEST output with information such as info/chip SNPs to be added
-   **KORA**. 91 (no BDNF) proteins, age, sex and individual level imputed genotypes
-   **LifeLinesDeep**. Only 1/25 proteins (unused)
-   **MadCam**. 91 (no IL.6) proteins
-   **ULSAM**. 25 proteins (unused)
-   **PIVUS**. 23 proteins (unused)
-   **ORCADES**. 91 (no BDNF) protein results are available but adding CCL3 which overlaps with MMP.1
-   **RECOMBINE**. 91 (no BDNF) protein results are available with information as described in the analysis plan
-   **VIS**. 91 (no BDNF) protein restults as with ORCADES
-   **STABILITY**. 90 (no BDNF, IL.2) protein.
-   **STANLEY**. 91 (no BDNF) largely complete protein results for lah1 and swe6


## NOTEs on results

The directory contains several versions of results,

Directory | Version of results
----------|----------------------------------
CEU | Results excluding 18 regions in high LD
HLA | Similar to CEU/ bove but keeps HLA
nold| No exclusion of regions
INTERVAL_nold | nold/ result for INTERVAL

The PLINK --clump-r2 0 results are included here for
comparison. Comprehensive information is contained
in .cis.vs.trans and cis/trans classificaiton in .out.
