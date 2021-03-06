# SCALLOP-INF meta-analysis

## Flow of analysis

(When the diagram is not rendered, view it from [Mermaid live editor](https://mermaid-js.github.io/mermaid-live-editor/))

```mermaid
graph TB;
  tryggve ==> cardio;
  cardio ==> csd3;
  csd3 --> csd3Analysis[Conditional analysis,finemapping, etc];
  csd3 --> software[R Packages at CRAN/GitHub]; 
  tryggveAnalysis[Meta analysis: list.sh, format.sh,metal.sh, QCGWAS.sh, analysis.sh] --> GWAS[pQTL selection and Characterization];
  GWAS --> Prototyping[Prototyping: INTERVAL.sh, cardio.sh, ...];
  Prototyping --> Multi-omics-analysis;
```

## Comments

The [tryggve](tryggve), [cardio](cardio) and [csd3](csd3) directories here are associated with the named Linux cluster(s) used for the analysis over time.

1. Data pre-processing was done initially from tryggve with [list.sh](tryggve/list.sh) and [format.sh](tryggve/format.sh), followed by meta-analysis according to [metal.sh](tryggve/metal.sh) using METAL whose results were cross-examined with [QCGWAS.sh](tryggve/QCGWAS.sh) together with additional investigation.

2. The main analysis followed with [analysis.sh](tryggve/analysis.sh) containing codes for Q-Q/Manhattan/LocusZoom/Forest plots such as the OPG example (see the diagram below), which replicated results of Kwan et al. (2014) as identified by PhenoScanner), clumping using PLINK and conditional analysis using GCTA. The clumping results were classified into cis/trans signals. As the meta-analysis stabilised especially with INTERVAL reference, analysis has been intensively done locally with cardio and csd3. cis/trans classification has been done via [cis.vs.trans.classification.R](cardio/cis.vs.trans.classification.R) as validated by [cistrans.sh](cardio/cistrans.sh).

3. We prototyped our analysis on cardio with INTERVAL such as [INTERVAL.sh](tryggve/INTERVAL.sh) and [cardio.sh](cardio/cardio.sh) as well as individual level data analysis for the KORA study.

4. The `cis.vs.trans.classification`, `circos.cis.vs.trans.plot` as with `cs`, `log10p`, `logp`, `gc.lambda`, `invnormal`, `METAL_forestplot`, `mhtplot.trunc`, `mhtplot2d`, `pvalue` functions are now part of R/gap at [CRAN](https://CRAN.R-project.org/package=gap) and updates such as `pqtl2dplot/pqtl2dplotly/pqtl3dplotly` are made at [GitHub](https://github.com/jinghuazhao/R/).

5. Downstream analyses links colocalisation and Mendelian randomisation with CAD, FEV1 and the meta-analysis summary statistics are now described in [pQTLtools articles](https://jinghuazhao.github.io/pQTLtools/articles/index.html).

6. A nested predictive model based on genotype data G, which link with proteins P1, P2, ..., Pn as predictors for outcome y. Alternative routes are T(P)WAS with [fusion_twas](http://gusevlab.org/projects/fusion/) and EWAS with [EWAS-fusion](https://jinghuazhao.github.io/EWAS-fusion/).
```mermaid
graph TD;
  G --> P1
  G --> P2
  G --> ...
  G --> Pn
  P1 --> y
  P2 --> y
  ... --> y
  Pn --> y

  SNP["LD reference panel (bed,bim,fam)"] --> |"EWAS reference panel(top1, blup, lasso, enet, bslmm)"| Methylation;
  Methylation --> Protein;
  SNP --> |"GWAS summary statistics (SNP, A1, A2, Z)"| Protein;
```

## The OPG example

This proves to be a positive control. The stacked image below shows Q-Q, Manhattan, LocusZoom and forest plots.

<p align="center"><img src="doc/OPG-qmlf.png"></p>

## EWAS with IL-12B

(EWAS, joint/conditional) Q-Q and Manhattan plots from `ewas-plot.R`.

<p align="center"><img src="doc/ewas-plot.png"></p>

## Summary statistics

The link will be added here when available.

## Related links

* [OlinkAnalyze](https://github.com/Olink-Proteomics/OlinkRPackage)
* [Olink Insights Stat Analysis](https://tinyurl.com/shj46ukj)
* [SCALLOP consortium](http://www.scallop-consortium.com/).
* [Olink location](https://www.olink.com/scallop/), [What is NPX](https://www.olink.com/question/what-is-npx/), [F2F London meeting](https://www.olink.com/scallop-f2f-2019/), [Data Science Workshop 2019](https://www.olink.com/data-science-workshop-2019/).
* [GitHub repository](https://github.com/lassefolkersen/scallop) for the 2017 *PLoS Genetics* paper above.
* [securecloud](https://secureremote.dtu.dk/vpn/index.html).
* [Olink publications](https://www.olink.com/data-you-can-trust/publications/).
* [SomaLogic plasma protein GWAS summary statistics](http://www.phpc.cam.ac.uk/ceu/proteins).
* [Aging Plasma Proteome](https://twc-stanford.shinyapps.io/aging_plasma_proteome/) ([DEswan](https://github.com/lehallib/DEswan)).
* [ImmunoBase](https://genetics.opentargets.org/immunobase).
* [Worldwide PDB](http://www.wwpdb.org/)

## References

Choi SW, Mak TS, O'Reilly PF Tutorial: a guide to performing polygenic risk score analyses. *Nat Protoc* 15, 2759-2772 (2020), [GitHub](https://github.com/choishingwan/PRSice) [documentation](https://choishingwan.github.io/PRS-Tutorial/).

Folkersen L, et al. (2017). Mapping of 79 loci for 83 plasma protein biomarkers in cardiovascular disease. *PLoS Genetics* 13(4), doi.org/10.1371/journal.pgen.1006706.

Kwan JSH, et al. (2014). Meta-analysis of genome-wide association studies identiﬁes two loci associated with circulating osteoprotegerin levels. *Hum Mol Genet* 23(24): 6684–6693.

Niewczas MA, et al. (2019). A signature of circulating inflammatory proteins and development of end-stage renal disease in diabetes. *Nat Med*. https://doi.org/10.1038/s41591-019-0415-5

Sun BB, et al. (2018). Genomic atlas of the human plasma proteome. *Nature* 558: 73–79.
