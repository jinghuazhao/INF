# SCAndinavian coLLaboration for Olink plasma Protein (SCALLOP)

INFlmmation anlaysis

The working directory is \$HOME/INF which can be made available from
/data or GitHub -- when the Internet connection is enabled one can deploy
with

``` {.bash}
git clone https://github.com/jinghua/INF
```

The repository contains the following files,

Filespec     | Descrription
-------------|-----------------------------------------
tryggve.md   | This document
tryggve.sh   | git batch file
list.sh      | Generation of file list and directory
format.sh    | Code for format GWAS summary statistics
analysis.sh  | Bash code for analysis calling analysis.ini
METAL.qsub   | TORQUE qsub script for METAL
metal.sh     | Generation/execution of METAL scripts
doc/         | Oringal documents
sumstats/    | File lists and study directories
METAL/       | METAL/output scripts by protein

In total, 92 proteins are expected as given in
[olink.prot.list.txt](doc/olink.prot.list.txt). Data processing so far
is done for all studies altogether; it might be appropriate to furnish
individually though this has been done by cut and paste interactively.

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

## KORA and INTERVAL

File | Description
-----|-------------------------------------------------------------
[KORA.sh](files/KORA.sh) | Bash/R scripts to hand KORA data
[KORA.R](files/KORA.R) | R code to simulate phenotypes
[KORA.prot.preproc.R*](doc/KORA.prot.preproc.R) | R code
[kora.normalised.prot.txt*](doc/kora.normalised.prot.txt) | sex, age, normalised proteins
[KORA.pc.below.llod.pdf*](doc/KORA.pc.below.llod.pdf) | llod check
[INTERVAL.sh](files/INTERVAL.sh) | INTERVAL analysis

\* from Jimmy

## Technical notes

The system is managed with modules, e.g.,

``` {.bash}
module load anaconda2/4.4.0
module load perl/5.24.0 annovar/2018apr16
module load bcftools/1.8
module load gcta/1.91.0beta
module load intel/redist/2018 intel/perflibs/64/2018 gcc/5.4.0 R/3.5.0-ICC-MKL rstudio/1.1.453
module load libreoffice/6.0.5.2
module load locuszoom/1.4
module load metal/20110325
module load pandoc/2.1
module load parallel/20170822
module load plink2/1.90beta5.4
module load qctool/1.4
module load vcftools/0.1.15
module load xpdf/3.04
export threads=1

xterm -fa arial -fs 12 -bg black -fg white
```

For instance, once libreoffice and pandoc are loaded, we can view this file with

``` {.bash}
pandoc SecureCloud.md -o SecureCloud.docx
soffice SecureCloud.docx
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
module load metal/20110325 parallel/20170822
ls METAL/*.run | parallel --dry-run --env HOME -j8 -C' ' 'metal $HOME/INF/{}'
```
NB METAL add -1 to the filenames.

## URLs

SCALLOP consortium, https://www.olink.com/scallop/

SCALLOP GitHub repository, https://github.com/lassefolkersen/scallop

SCALLOP securecloud, https://secureremote.dtu.dk/vpn/index.html

SomaLogic plasma protein GWAS summary statistics, http://www.phpc.cam.ac.uk/ceu/proteins

96, MMP.12 (CVD II)
97, MMP.9 (CVD III)
153, TNFRSF12A (UnitProt ID Q9NP84) (reported in Sun et al. but nonoverlapped with CVD/INF)
     ~ MMP.9 connection on PLoS One by Yang et al. (2018)
