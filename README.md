# SCALLOP-INF meta-analysis

## Flow of analysis (view diagram by pasting script to [mermaid live editor](https://mermaid-js.github.io/mermaid-live-editor/))

```mermaid
graph TB;
tryggve ==> cardio;
cardio ==> csd3;
tryggveAnalysis[Meta analysis: list.sh, format.sh,metal.sh, QCGWAS.sh, analysis.sh] --> GWAS[pQTL selection and Characterisation];
GWAS --> Prototyping[Prototyping: INTERVAL.sh, cardio.sh, ...];
Prototyping --> Multi-omics-analysis;
cardio --> cardioAnalysis[Prototyping and KORA data analysis];
csd3 --> csd3Analysis[Conditional analysis,finemapping, etc];
csd3 --> software[R Packages at CRAN/GitHub]; 
```

### Comments

1. Data pre-processing was done initially from [tryggve](tryggve) with [list.sh](tryggve/list.sh) and [format.sh](tryggve/format.sh), followed by meta-analysis according to [metal.sh](tryggve/metal.sh) using METAL whose results were cross-examined with [QCGWAS.sh](tryggve/QCGWAS.sh) together with additional investigation.
2. The main analysis followed with [analysis.sh](tryggve/analysis.sh) containing codes for Q-Q/Manhattan/LocusZoom/Forest plots such as the OPG example (see the diagram below), which replicated results of Kwan et al. (2014) as identified by PhenoScanner), clumping using PLINK and conditional analysis using GCTA. The clumping results were classified into cis/trans signals. As the meta-analysis stabilised especially with INTERVAL reference, analysis has been intensively done locally with [cardio](cardio) and [CSD3](csd3). cis/trans classification has been done via [cis.vs.trans.classification.R](cardio/cis.vs.trans.classification.R) as validated by [cistrans.sh](cardio/cistrans.sh).
3. We prototyped our analysis on cardio with INTERVAL such as [INTERVAL.sh](tryggve/INTERVAL.sh) and [cardio.sh](cardio/cardio.sh) as well as individual level data analysis for the KORA study.
4. The `cis.vs.trans.classification`, `circos.cis.vs.trans.plot` as with `cs`, `log10p`, `logp`, `gc.lambda`, `invnormal`, `METAL_forestplot`, `mhtplot.trunc`, `mhtplot2d`, `pvalue` functions are now part of R/gap at [CRAN](https://CRAN.R-project.org/package=gap) and updates such as `pqtl2dplot/pqtl2dplotly/pqtl3dplotly` are made at [GitHub](https://github.com/jinghuazhao/R/).
5. Downstream analyses links colocalisation and Mendelian randomisation with CAD, FEV1 and the meta-analysis summary statistics are now described in [pQTLtools articles](https://jinghuazhao.github.io/pQTLtools/articles/index.html).
6. A nested predictive model based on genotype data G1, G2, G3, which link with proteins P1, P2, P3 as predictors for outcome y is sketched as follows,

![](rsid/grViz.png)

## References

Choi SW, Mak TS, O'Reilly PF Tutorial: a guide to performing polygenic risk score analyses. *Nat Protoc* 15, 2759-2772 (2020), [GitHub](https://github.com/choishingwan/PRSice) [documentation](https://choishingwan.github.io/PRS-Tutorial/).

Folkersen L, et al. (2017). Mapping of 79 loci for 83 plasma protein biomarkers in cardiovascular disease. *PLoS Genetics* 13(4), doi.org/10.1371/journal.pgen.1006706.

Kwan JSH, et al. (2014). Meta-analysis of genome-wide association studies identiﬁes two loci associated with circulating osteoprotegerin levels. *Hum Mol Genet* 23(24): 6684–6693.

Niewczas MA, et al. (2019). A signature of circulating inflammatory proteins and development of end-stage renal disease in diabetes. *Nat Med*. https://doi.org/10.1038/s41591-019-0415-5

Sun BB, et al. (2018). Genomic atlas of the human plasma proteome. *Nature* 558: 73–79.

## Download

```{bash}
git clone https://github.com/jinghua/INF
```

## Study Information

* [Analysis plan](doc/SCALLOP_INF1_analysis_plan.md) ([docx](doc/SCALLOP_INF1_analysis_plan.docx))
* [Approximately independent LD blocks](doc/aild.md)
* [Joint/conditional analysis and fine-mapping](rsid/README.md)
* [Notes on UniProt IDs](doc/uniprot.md)
* [TRYGGVE](https://neic.no/tryggve/)-[specific notes](tryggve/tryggve.md)

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

The OPG example,<img src="doc/OPG-qmlf.png">
