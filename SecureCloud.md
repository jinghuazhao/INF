# SCAndinavian coLLaboration for Olink plasma Protein (SCALLOP)

INFlmmation anlaysis

The working directory is \$HOME/INF which can be made available from
/data or GitHub -- once the Internet connection is enabled one can deloy
with

``` {.bash}
git clone https://github.com/jinghua/INF
```

The repository contains the following files,

Filespec     | Descrription
-------------|-----------------------------------------
SecureCloud.md    | This document
[list.sh](list.sh)      | Generation of file list and directory
[format.sh](format.sh)    | Code for format GWAS summary statistics
doc/         | Orignal documents
sumstats/    | Directory for file list and study
METAL/       | METAL/output scripts by protein
[METAL.qsub](METAL.qsub)   | TORQUE qsub script for METAL

In total, 92 proteins are expected as given in
[olink.prot.list.txt](doc/olink.prot.list.txt). Data processing so far
is done for all studies altogether; it might be appropriate to furnish
individually though this has been done by cut and paste.interactively.

## Notes on studies

-   **NSPHS**. A tar.gz file is unpacked into \$HOME/INF/work leading to 10 proteins
-   **EGCUT**. Only 18 proteins though stratified by chromsomes
-   **Estonian)biobank**. 18 proteins as with EGCUT
-   **ULSAM**. 25 proteins
-   **PIVUS**. 23 proteins
-   **INTERVAL**. raw SNPTEST output with information such as strand/chip SNPs to be added
-   **LifeLinesDeep**. Only 1/25 proteins
-   **ORCADES**. 91 protein results are available but adding CCL3 which overlaps with MMP.1
-   **VIS**. 91 protein restults as with ORCADES

It is necessary to generate an annotation file similar to [cvd1\_annotation](doc/cvd1_annotation).

## Technical notes

The system is managed with modules, e.g.,

``` {.bash}
# R
module avail R
module unload R
module load gcc/5.4.0
module load R/3.2.5
module load rstudio/1.1.453

module avail gcta
module load gcta/1.91.0beta
module load libreoffice/6.0.5.2
module load metal/20110325
module load pandoc/2.1
module load parallel/20170822
module load plink2/1.90beta5.4
module load xpdf/3.04

xterm -fa arial -fs 12 -bg black -fg white
```

For instance, once libreoffice and pandoc are loaded, we can view this
file with

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
ls METAL/*.run | parallel --dry-run --env HOME -j3 -C' ' 'metal $HOME/INF/{}'
```

## URLs

SCALLOP consortium, https://www.olink.com/scallop/

SCALLOP GitHub repository, https://github.com/lassefolkersen/scallop

SCALLOP securecloud, https://secureremote.dtu.dk/vpn/index.html

SomaLogic plasma protein GWAS summary statistics, http://www.phpc.cam.ac.uk/ceu/proteins

96, MMP.12 (CVD II)
97, MMP.9 (CVD III)
153, TNFRSF12A (UnitProt ID Q9NP84) (reported in Sun et al. but nonoverlapped with CVD/INF)
     ~ MMP.9 connection on PLoS One by Yang et al. (2018)
