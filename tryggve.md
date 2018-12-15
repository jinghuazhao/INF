# SCALLOP/INFlammation anlaysis on TRYGGVE

The repository contains the following files,

File specification | Descrription
-------------------|---------------------------------------------------------------
doc/               | Oringal documents
 -- [KORA.prot.preproc.R*](doc/KORA.prot.preproc.R) | R code for data preprocessing
 -- kora.normalised.prot.txt* | sex, age, normalised proteins
 -- [KORA.pc.below.llod.pdf*](doc/KORA.pc.below.llod.pdf) | llod check
METAL/             | METAL/output scripts by protein
sumstats/          | File lists and study directories
tryggve/           | Auxiliary files
 -- [list.sh](tryggve/list.sh)         | Generation of file list and directory
 -- [format.sh](tryggve/format.sh)     | Code for format GWAS summary statistics
 -- [lz14.sh](tryggve/lz14.sh)         | Code to extract LocusZoom 1.4 databases for [analysis.sh](tryggve/analysis.sh)
 -- [analysis.sh](tryggve/analysis.sh) | Bash code for analysis calling [analysis.ini](tryggve/analysis.ini)
 -- [metal.sh](tryggve/metal.sh)       | Generation/execution of METAL scripts
 -- [METAL.qsub](tryggve/METAL.qsub)   | TORQUE qsub script for METAL
 -- [INTERVAL.sh](tryggve/INTERVAL.sh) | INTERVAL analysis
 -- [KORA.sh](tryggve/KORA.sh) | Bash/R scripts to hand KORA data
 -- [KORA.R](tryggve/KORA.R) | R code to simulate phenotypes
 -- [EURLD.bed](tryggve/EURLD.bed) | approximately independent LD blocks for Europeans.
tryggve.md         | This document

\* from Jimmy

In total, 92 proteins are expected as given in
[olink.prot.list.txt](doc/olink.prot.list.txt). Data processing so far
is done for all studies altogether; it might be appropriate to furnish
individually though this has been done by cut and paste interactively.

## Modules and related applications

Instances are exemplified as follows,

``` {.bash}
module load anaconda2/4.4.0
module load perl/5.24.0 annovar/2018apr16
module load bcftools/1.8
module load gcta/1.91.0beta
module load intel/redist/2019 intel/perflibs/64/2019 gcc/5.4.0 R/3.5.0-ICC-MKL rstudio/1.1.453
module load libreoffice/6.0.5.2
module load locuszoom/1.4
module load metal/20110325
module load pandoc/2.1
module load parallel/20170822
module load plink2/1.90beta5.4
module load vcftools/0.1.15
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

[METAL.qsub](tryggve/METAL.qsub) was intended to run batch jobs under [Terascale Open-source
Resource and QUEue Manager (TORQUE)](https://en.wikipedia.org/wiki/TORQUE) controls batch jobs and
distributed compute nodes, enabling integration with [Moab cluster
suite](https://en.wikipedia.org/wiki/Moab_Cluster_Suite) and extending
[Portable Batch System (PBS)](https://en.wikipedia.org/wiki/Portable_Batch_System) with respect
to extend scalability, fault tolerance, and functionality. Information on
job scheduling on Computerome is availble from https://www.computerome.dk/display/CW/Batch+System,

However, the SCALLOP securitycloud is a 1-node cloud and we resort to
parallel instead, e.g,

```bash
module load metal/20110325 parallel/20170822
ls METAL/*.run | parallel --dry-run --env HOME -j8 -C' ' 'metal $HOME/INF/{}'
```
NB METAL add -1 to the filenames.

## INSTALLATIONS

There are a number of software updates/additions which are worthy of note.

### GCTA

The version that can handle z statistics instead of p values is at /data/jinhua/gcta_1.91.7beta/

### METAL

The version contains modification which allows for CUSTOMVARIABLE to use integer position rather than scientific format as in [software-notes](https://github.com/jinghuazhao/software-notes/blob/master/AA.md).

### GNU parallel

The version has some improvement and is installed at /data/jinhua/parallel-20181022/.

### qctool

The system has module qctool/2.0.1 installed but nevertheless not functional for lack of lapack shared libraries, so lapack-3.8.0/ is available from /data/jinhua as composed to module lapack/3.8.0.

### R/gap

The version at TRYGGVE and CRAN are not the latest, which contains functions cis.vs.trans.classification, gc.lambda, invnormal, and here is the way to go
```bash
module load intel/redist/2019 intel/perflibs/64/2019 gcc/5.4.0 R/3.5.0-ICC-MKL
tar xvfz gap_1.1-23.tar.gz
cd gap/src
gcc -I/services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/include -L/services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/lib/release -L/services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/lib -Xlinker --enable-new-dtags -Xlinker -rpath -Xlinker /services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/lib/release -Xlinker -rpath -Xlinker /services/tools/intel/perflibs/2019/compilers_and_libraries/linux/mpi/intel64/lib -Xlinker -rpath -Xlinker /opt/intel/mpi-rt/2017.0.0/intel64/lib/release -Xlinker -rpath -Xlinker /opt/intel/mpi-rt/2017.0.0/intel64/lib -lmpifort -lmpi -ldl -lrt -lpthread -L/services/tools/intel/perflibs/2019//compilers_and_libraries_2019.0.117/linux/mpi/intel64/libfabric/lib -fPIC -c *.c *.f
gcc -shared -L/services/tools/R/3.5.0-ICC-MKL/lib64/R/lib -L/usr/local/lib64 -o gap.so 2k.o 2ld.o cline.o gcontrol_c.o gcx.o gif_c.o hap_c.o hwe.hardy.o kin.morgan.o makeped_c.o mia.o muvar.o package_native_routine_registration_skeleton.o pfc.o pfc.sim.o pgc_c.o whscore.o -L/usr/lib/gcc/x86_64-redhat-linux/4.8.2 -lgfortran -lm -lquadmath -L/services/tools/R/3.5.0-ICC-MKL/lib64/R/lib -lR
cd -
R CMD INSTALL gap
```
The rather imtimidating compiler flags are derived from `mpiicc -show` as with initial pass requested for libfabric.so.1, and the idea is to get around check for icc which does not exist on the system.

### R/QCGWAS

The version contains fix to the use of HapMap reference as in software-notes above. The HapMap data as with code from the packages's quick guide is /data/jinhua/data/QCGWAS.

### R/EasyQC and R/EasyStrata

The version is 18.1 rather than 9.2 and 8.6 currenly online, https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/software/.

## URLs

SCALLOP consortium, https://www.olink.com/scallop/

SCALLOP GitHub repository, https://github.com/lassefolkersen/scallop

SCALLOP securecloud, https://secureremote.dtu.dk/vpn/index.html

SomaLogic plasma protein GWAS summary statistics, http://www.phpc.cam.ac.uk/ceu/proteins

## Notes on studies

-   **NSPHS**. 91 proteins, originally a tar.gz file is unpacked into \$HOME/INF/work leading to 10 proteins
-   **EGCUT**. 91 proteins, orginally only 18 proteins though stratified by chromsomes
-   **INTERVAL**. raw SNPTEST output with information such as strand/chip SNPs to be added
-   **LifeLinesDeep**. Only 1/25 proteins
-   **ULSAM**. 25 proteins
-   **PIVUS**. 23 proteins
-   **ORCADES**. 91 protein results are available but adding CCL3 which overlaps with MMP.1
-   **VIS**. 91 protein restults as with ORCADES
-   **STABILITY**. 90 protein.
-   **STANLEY**. 92 largely complete protein results for lah1 and swe6

It is necessary to generate an annotation file similar to [cvd1\_annotation.tsv](doc/cvd1_annotation.tsv).

96, MMP.12 (CVD II)
97, MMP.9 (CVD III)
153, TNFRSF12A (UnitProt ID Q9NP84) (reported in Sun et al. but nonoverlapped with CVD/INF)
     ~ MMP.9 connection on PLoS One by Yang et al. (2018)
