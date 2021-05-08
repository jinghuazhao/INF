<img src="doc/INF1.circlize.png" width="300" height="300" align="right">
<script>
  valine: 'https://unpkg.com/valine/dist/Valine.min.js',
  mathjax: 'https://cdn.bootcss.com/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML',
  mermaid: 'https://cdn.bootcss.com/mermaid/8.0.0-rc.8/mermaid.min.js'
</script>
# SCALLOP-INF meta-analysis

## Flow of analysis

<pre><code class="language-mermaid">
graph TB;
tryggve ==> cardio;
cardio ==> csd3;
tryggveAnalysis[Meta analysis: list.sh, format.sh,metal.sh, QCGWAS.sh, analysis.sh] --> GWAS[pQTL selection and Characterisation];
GWAS --> Prototyping[Prototyping: INTERVAL.sh, cardio.sh, ...];
Prototyping --> Multi-omics-analysis;
cardio --> cardioAnalysis[Prototyping and KORA data analysis];
csd3 --> csd3Analysis[Conditional analysis,finemapping, etc];
csd3 --> software[R Packages at CRAN/GitHub]; 
</code></pre>

### Comments

1. Data pre-processing was done initially from [tryggve](tryggve) with [list.sh](tryggve/list.sh) and [format.sh](tryggve/format.sh), followed by meta-analysis according to [metal.sh](tryggve/metal.sh) using METAL whose results were cross-examined with [QCGWAS.sh](tryggve/QCGWAS.sh) together with additional investigation.
2. The main analysis followed with [analysis.sh](tryggve/analysis.sh) containing codes for Q-Q/Manhattan/LocusZoom/Forest plots such as the OPG example (see the diagram below), which replicated results of Kwan et al. (2014) as identified by PhenoScanner), clumping using PLINK and conditional analysis using GCTA. The clumping results were classified into cis/trans signals. As the meta-analysis stabilised especially with INTERVAL reference, analysis has been intensively done locally with [cardio](cardio) and [CSD3](csd3). cis/trans classification has been done via [cis.vs.trans.classification.R](cardio/cis.vs.trans.classification.R) as validated by [cistrans.sh](cardio/cistrans.sh).
3. We prototyped our analysis on cardio with INTERVAL such as [INTERVAL.sh](tryggve/INTERVAL.sh) and [cardio.sh](cardio/cardio.sh) as well as individual level data analysis for the KORA study.
4. The `cis.vs.trans.classification`, `circos.cis.vs.trans.plot` (as in this page) as with `cs`, `log10p`, `logp`, `gc.lambda`, `invnormal`, `METAL_forestplot`, `mhtplot.trunc`, `mhtplot2d`, `pvalue` functions are now part of R/gap at [CRAN](https://CRAN.R-project.org/package=gap) and updates are made at [GitHub](https://github.com/jinghuazhao/R/).
5. Downstream analysis links examples such as CAD summary statistics for MR and MAGMA [here](https://jinghuazhao.github.io/Omics-analysis/CAD/), colocalisation analysis of simulated data in [software notes on association analysis](https://jinghuazhao.github.io/software-notes/AA.html) as well as the [BMI example](https://jinghuazhao.github.io/Omics-analysis/BMI/).

## References

Choi WM, Mak TSH, O’Reilly PF (2018). A guide to performing Polygenic Risk Score analyses. [a tutorial](https://choishingwan.github.io/PRS-Tutorial/) ([GitHub](https://github.com/choishingwan/PRSice)),
https://www.biorxiv.org/content/10.1101/416545v1.

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
* [pQTLtools articles](https://jinghuazhao.github.io/pQTLtools/articles/index.html)
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
