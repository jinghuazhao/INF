# SCALLOP-INF meta-analysis

A companion web site for this paper, 

Jing Hua Zhao, David Stacey, Niclas Eriksson, Erin Macdonald-Dunlop, Asa H Hedman, Anette Kalnapenkis, Stefan Enroth, Domenico Cozzetto, Jonathan Digby-Bell, Jonanthan Marten, Lasse Folkersen, Christian Herder, Lina Jonsson, Sarah E. Bergen, Christian Gieger, Elise J Needham, Praveen Surendran, Estonia Biobank Research Team, Dirk S Paul, Ozren Polasek, Barbara Thorand, Harald Grallert, Michael Roden, Urmo Vosa, Tonu Esko, Caroline Hayward, Asa Johansson, Ulf Gyllensten, Nicholas Powell, Oskar Hansson, Niklas Mattsson-Carlgren, Peter K Joshi, John Danesh, Leonid Padyukov, Lars Klareskog, Mikael Landen, James F Wilson, Agneta Siegbahn, Lars Wallentin, Anders Malarstig, Adam S Butterworth, James E. Peters
medRxiv 2023.03.24.23287680; doi: https://doi.org/10.1101/2023.03.24.23287680, <https://www.medrxiv.org/content/10.1101/2023.03.24.23287680v1>.
<img src="https://connect.medrxiv.org/qr/qr_img.php?id=2023.03.24.23287680" width=60 height=60>.

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

The [tryggve](tryggve), [cardio](cardio) and [csd3](csd3) directories here are associated with the named Linux cluster(s) used for the analysis over time. Most recent implementations are documented in `Supplementary notes` ([rsid](rsid)). <font color="blue"><b>To view the code inside the browser, select the `GitHub` button from the menu.</b></font>

1. Data pre-processing was done initially from tryggve with [list.sh](tryggve/list.sh) and [format.sh](tryggve/format.sh), followed by meta-analysis according to [metal.sh](tryggve/metal.sh) using METAL whose results were cross-examined with [QCGWAS.sh](tryggve/QCGWAS.sh) together with additional investigation.

2. The main analysis followed with [analysis.sh](tryggve/analysis.sh) containing codes for Manhattan/Q-Q/forest/LocusZoom plots such as the OPG example (see the diagram below), which replicated results of Kwan et al. (2014) as identified by PhenoScanner), clumping using PLINK and conditional analysis using GCTA. The clumping results were classified into cis/trans signals. As the meta-analysis stabilised especially with INTERVAL reference, analysis has been intensively done locally with cardio and csd3. cis/trans classification has been done via [cis.vs.trans.classification.R](cardio/cis.vs.trans.classification.R) as validated by [cistrans.sh](cardio/cistrans.sh).

3. We prototyped our analysis on cardio with INTERVAL such as [INTERVAL.sh](tryggve/INTERVAL.sh) and [cardio.sh](cardio/cardio.sh) as well as individual level data analysis for the KORA study. Most analyses were done locally on CSD3.

4. In alphabetical order, the `cis.vs.trans.classification`, `circos.cis.vs.trans.plot`, `cs`, `log10p`, `logp`, `gc.lambda`, `get_b_se`, `get_sdy`, `get_pve_se`, `invnormal`, `METAL_forestplot`, `mhtplot.trunc`, `mhtplot2d`, `pqtl2dplot/pqtl2dplotly/pqtl3dplotly`, `pvalue` functions are now part of R/gap at [CRAN](https://CRAN.R-project.org/package=gap) or [GitHub](https://github.com/jinghuazhao/R/). Some aspects of the downstream analyses links colocalisation and Mendelian randomisation are also available from [gap vignette](https://jinghuazhao.github.io/R/vignettes/gap.html) and [pQTLtools articles](https://jinghuazhao.github.io/pQTLtools/articles/index.html).

## The OPG example

This proves to be a positive control. The stacked image below shows Manhattan, Q-Q, forest and LocusZoom plots.

<p align="center"><img src="doc/OPG.png"></p>

## Summary statistics

The link will be added here when made available.

<p align="center"><img src="doc/circos.png"></p>
