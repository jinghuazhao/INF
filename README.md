<img src="doc/INF1.circlize.png" width="300" height="300" align="right">

---

# SCALLOP Olink/INFlammation analysis

---

## TABLE OF CONTENTS ([Web site](https://jinghuazhao.github.io/INF/))

* [Analysis plan](doc/SCALLOP_INF1_analysis_plan.md) ([docx](doc/SCALLOP_INF1_analysis_plan.docx))
* Analysis
* References
* Additional information
  * [Approximately independent LD blocks](doc/aild.md)
  * [Competitive logp/log10p functions](doc/logplog10p.md)
  * [Notes on UniProt IDs](doc/uniprot.md)
  * [TRYGGVE](https://neic.no/tryggve/)-[specific issues](tryggve/tryggve.md)
  * URLs and downloading

<img src="doc/OPG-qmlf.png" width="400" height="800" align="right">

## Analysis

1. To implement the analysis plan, we started with analysis on INTERVAL as with [INTERVAL.sh](tryggve/INTERVAL.sh) and [cardio.sh](cardio/cardio.sh).

2. Data preprocessig was done with [list.sh](tryggve/list.sh) and [format.sh](tryggve/format.sh). The meta-analysis was according to [metal.sh](tryggve/metal.sh) using METAL whose results were cross-examined with [QCGWAS.sh](tryggve/QCGWAS.sh) together with addtional investigation.

3. The main analysis follows [analysis.sh](tryggve/analysis.sh), initially from [tryggve](tryggve), which contains codes for Q-Q/Manhattan/LocusZoom/Forest plots (see the figure on the right for the OPG example, which replicated results of Kwan et al. (2014) as identified by PhenoScanner), clumping using PLINK and conditional analysis using GCTA. The clumping results were classified into cis/trans signals. As the meta-analysis stabilised especially with INTERVAL reference, analysis has been intensively done locally with [cardio](cardio) and [CSD3](csd3). cis/trans classificaiton has been done via [cis.vs.trans.classification.R](cardio/cis.vs.trans.classification.R) as validated by [cistrans.sh](cardio/cistrans.sh).

4. The `cis.vs.trans.classification`, `circos.cis.vs.trans.plot` (as in this page) as with `log10p`, `gc.lambda`, `invnormal`, `METAL_forestplot`, `mhtplot.trunc`, `mhtplot2d` functions are now part of R/gap at [CRAN](https://CRAN.R-project.org/package=gap) and updates are made at [GitHub](https://github.com/jinghuazhao/R/tree/master/gap).

5. Downstream analysis links examples such as CAD summary statistics for MR and MAGMA [here](https://jinghuazhao.github.io/Omics-analysis/CAD/), colocalisation analysis of simulated data in [software notes on association analysis](https://jinghuazhao.github.io/software-notes/AA.html) as well as the [BMI example](https://jinghuazhao.github.io/Omics-analysis/BMI/).

## References

Choi WM, Mak TSH, O’Reilly PF (2018). A guide to performing Polygenic Risk Score analyses. [a tutorial](https://choishingwan.github.io/PRS-Tutorial/) ([GitHub](https://github.com/choishingwan/PRSice)),
https://www.biorxiv.org/content/10.1101/416545v1.

Folkersen L, et al. (2017). Mapping of 79 loci for 83 plasma protein biomarkers in cardiovascular disease. *PLoS Genetics* 13(4), doi.org/10.1371/journal.pgen.1006706.

Kwan JSH, et al. (2014). Meta-analysis of genome-wide association studies identiﬁes two loci associated with circulating osteoprotegerin levels. *Hum Mol Genet* 23(24): 6684–6693.

Niewczas MA, et al. (2019). A signature of circulating inflammatory proteins and development of end-stage renal disease in diabetes. *Nat Med*. https://doi.org/10.1038/s41591-019-0415-5

Sun BB, et al. (2018). Genomic atlas of the human plasma proteome. *Nature* 558: 73–79.

## URLs and downloading

### URLs

* [SCALLOP consortium](http://www.scallop-consortium.com/).
* [Olink location](https://www.olink.com/scallop/), [What is NPX](https://www.olink.com/question/what-is-npx/), [F2F London meeting](https://www.olink.com/scallop-f2f-2019/), [Data Science Workshop 2019](https://www.olink.com/data-science-workshop-2019/).
* [GitHub repository](https://github.com/lassefolkersen/scallop) for the 2017 *PLoS Genetics* paper above.
* [Tryggve](https://neic.no/tryggve/) and [securecloud](https://secureremote.dtu.dk/vpn/index.html).
* [Olink publications](https://www.olink.com/data-you-can-trust/publications/).
* [SomaLogic plasma protein GWAS summary statistics](http://www.phpc.cam.ac.uk/ceu/proteins).
* [Aging Plasma Proteome](https://twc-stanford.shinyapps.io/aging_plasma_proteome/) ([DEswan](https://github.com/lehallib/DEswan)).
* [ImmunoBase](https://genetics.opentargets.org/immunobase).

### downloading

```{bash}
git clone https://github.com/jinghua/INF
```

