# SCALLOP-INF meta-analysis

A companion web site for this paper[^medRxiv], 

Zhao, J.H., Stacey, D., Eriksson, N., Macdonald-Dunlop, E., Hedman, Å.K., Kalnapenkis, A., Enroth, S., Cozzetto, D., Digby-Bell, J., Marten, J., Folkersen, L., Herder, C., Jonsson, L., Bergen, S.E., Gieger, C., Needham, E.J., Surendran, P., Metspalu, A., Milani, L., Mägi, R., Nelis, M., Hudjašov, G., Paul, D.S., Polasek, O., Thorand, B., Grallert, H., Roden, M., Võsa, U., Esko, T., Hayward, C., Johansson, Å., Gyllensten, U., Powell, N., Hansson, O., Mattsson-Carlgren, N., Joshi, P.K., Danesh, J., Padyukov, L., Klareskog, L., Landén, M., Wilson, J.F., Siegbahn, A., Wallentin, L., Mälarstig, A., Butterworth, A.S., Peters, J.E., and Estonian Biobank Research Team (2023). **Genetics of circulating inflammatory proteins identifies drivers of immune-mediated disease risk and therapeutic targets**. *Nature Immunology*,
URL <https://www.nature.com/articles/s41590-023-01588-w>.

Quick links to codes for figures,

- Figure 1, [circos2.R](https://github.com/jinghuazhao/INF/blob/master/rsid/circos2.R)
- Figure 2, [hotspot.sh](https://github.com/jinghuazhao/INF/blob/master/csd3/hotspot.sh), [utils.sh](https://github.com/jinghuazhao/INF/blob/master/rsid/utils.sh), [IL.12B.sh](https://github.com/jinghuazhao/INF/blob/master/csd3/IL.12B.sh), [TRAIL.sh](https://github.com/jinghuazhao/INF/blob/master/csd3/TRAIL.sh)
- Figure 3, [IL.18-rs385076.sh](https://github.com/jinghuazhao/INF/blob/master/rsid/IL.18-rs385076.sh)
- Figure 4,
- Figure 5, [gsmr.r](https://github.com/jinghuazhao/INF/blob/master/workflow/scripts/gsmr.r)
- Figure 6, [utils.sh](https://github.com/jinghuazhao/INF/blob/master/rsid/utils.sh)

- Extended Data Figure 1,
- Extended Data Figure 2, [utils.sh](https://github.com/jinghuazhao/INF/blob/master/rsid/utils.sh), [IL.17C.R](https://github.com/jinghuazhao/INF/blob/master/rsid/IL.17C.R)
- Extended Data Figure 3, [aristotle.sh](https://github.com/jinghuazhao/INF/blob/master/csd3/aristotle.sh)
- Extended Data Figure 4, [h2pve.R](https://github.com/jinghuazhao/INF/blob/master/rsid/h2pve.R)
- Extended Data Figure 5, [rs12075.R](https://github.com/jinghuazhao/INF/blob/master/rsid/rs12075.R)
- Extended Data Figure 6, [utils.sh](https://github.com/jinghuazhao/INF/blob/master/rsid/utils.sh)
- Extended Data Figure 7,
- Extended Data Figure 8, [pqtlGWAS.R](https://github.com/jinghuazhao/INF/blob/master/rsid/pqtlGWAS.R)
- Extended Data Figure 9, [pqtlGWAS.R](https://github.com/jinghuazhao/INF/blob/master/rsid/pqtlGWAS.R)
- Extended Data Figure 10, [utils.sh](https://github.com/jinghuazhao/INF/blob/master/rsid/utils.sh)

- [Supplementary Figure 1](doc/manhattan-qq.pdf) (450dpi), [qqmanhattanlz.sb](https://github.com/jinghuazhao/INF/blob/master/rsid/qqmanhattanlz.sb), [utils.sh](https://github.com/jinghuazhao/INF/blob/master/rsid/utils.sh)
- [Supplementary Figure 2](doc/fp-lz.pdf), [utils.sh](https://github.com/jinghuazhao/INF/blob/master/rsid/utils.sh)
- [Supplementary Figure 3](doc/eQTLGen-INF.pdf), [eQTLGen.sh](https://github.com/jinghuazhao/INF/blob/master/rsid/eQTLGen.sh)
- [Supplementary Figure 4](doc/GWAS-INF.pdf), [coloc-disease.sh](https://github.com/jinghuazhao/INF/blob/master/rsid/coloc-disease.sh)

- Supplementary item, [js.R](https://github.com/jinghuazhao/INF/blob/master/csd3/js.R), [merge.sh](https://github.com/jinghuazhao/INF/blob/master/csd3/merge.sh)

- Supplementary Tables, [tables.R](https://github.com/jinghuazhao/INF/blob/master/rsid/tables.R)


## Flow of analysis

The diagram can also be rendered via [Mermaid live editor](https://mermaid-js.github.io/mermaid-live-editor/).

```mermaid
graph TB;
  tryggve ==> cardio;
  cardio ==> csd3;
  csd3 --> csd3Analysis[Conditional analysis, finemapping, etc];
  csd3 --> software[R Packages at CRAN/GitHub]; 
  tryggveAnalysis[Meta analysis: list.sh, format.sh,metal.sh, QCGWAS.sh, analysis.sh] --> GWAS[pQTL selection and Characterization];
  GWAS --> Prototyping[Prototyping: INTERVAL.sh, cardio.sh, ...];
  Prototyping --> Multi-omics-analysis;
```

## Comments

<font color="blue"><b>To view the code inside the browser, select the `GitHub` button from the menu</b></font>.

The [tryggve](tryggve), [cardio](cardio) and [csd3](csd3) directories here are associated with the named Linux cluster(s) used for the analysis over time. Early implementation involves the following aspects,

1. Data pre-processing from tryggve with [list.sh](tryggve/list.sh) and [format.sh](tryggve/format.sh), followed by meta-analysis according to [metal.sh](tryggve/metal.sh) using METAL whose results were cross-examined with [QCGWAS.sh](tryggve/QCGWAS.sh) together with additional investigation.
2. The main analysis with [analysis.sh](tryggve/analysis.sh) containing codes for Manhattan/Q-Q/forest/LocusZoom plots, clumping using PLINK and conditional analysis using GCTA. The clumping results were classified into cis/trans signals. As the meta-analysis stabilised especially with INTERVAL reference, analysis has been intensively done locally with cardio and csd3. cis/trans classification has been done via [cis.vs.trans.classification.R](cardio/cis.vs.trans.classification.R) as validated by [cistrans.sh](cardio/cistrans.sh).
3. Prototyping analysis on cardio with INTERVAL such as [INTERVAL.sh](tryggve/INTERVAL.sh) and [cardio.sh](cardio/cardio.sh) as well as individual level data analysis for the KORA study. Most analyses were done locally on CSD3.

Most recent implementations are documented in `Supplementary notes` ([rsid](rsid)). Over time, many functions become part of two R packages, gap ([CRAN](https://CRAN.R-project.org/package=gap), [GitHub](https://github.com/jinghuazhao/R/), [vignette](https://jinghuazhao.github.io/R/vignettes/gap.html)) and pQTLtools ([Web page](https://jinghuazhao.github.io/pQTLtools/)).

## A benchmark

As revealed by PhenoScanner[^phenoscanner], Osteoprotegerin (OPG) proves to be a positive control[^OPG] involving both a cis and a trans pQTLs, with the cis-pQTL showing stronger association -- see a stacked image containg forest+LocusZoom and Manhattan+Q-Q plots, [OPG.pdf](doc/OPG.pdf).

## Summary statistics

They will be available from 

* Cardiovascular Epidemiology Unit (CEU), <https://www.phpc.cam.ac.uk/ceu/proteins>
* GWAS catalog, <https://www.ebi.ac.uk/gwas/publications/37563310> (accession GCST90274758-GCST90274848 [[gcst_list.xlsx](doc/gcst_list.xlsx)], [protein-target-gene mapping](https://github.com/jinghuazhao/INF/blob/master/doc/prot_target_gene.tsv))
    - <https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/>
    - <https://www.ebi.ac.uk/gwas/studies/GCST90274758><pre>&#x2026;</pre><https://www.ebi.ac.uk/gwas/studies/GCST90274848>

<p align="center"><img src="doc/circos.svg"></p>

The diagram is based on [circos](http://circos.ca/) in the named directory [circos](https://github.com/jinghuazhao/INF/tree/master/circos) which highlights pQTLs [causal genes]; the significance levels of association can be seen from the inner scatter plot whose ceiling for -log10(P) is set to be 150.

## References

<img src="https://connect.medrxiv.org/qr/qr_img.php?id=2023.03.24.23287680" align="right" width=80 height=80>.

[^medRxiv]: The SCALLOP consortium. Jing Hua Zhao, David Stacey, Niclas Eriksson, Erin Macdonald-Dunlop, Asa H Hedman, Anette Kalnapenkis, Stefan Enroth, Domenico Cozzetto, Jonathan Digby-Bell, Jonanthan Marten, Lasse Folkersen, Christian Herder, Lina Jonsson, Sarah E. Bergen, Christian Gieger, Elise J Needham, Praveen Surendran, Estonia Biobank Research Team, Dirk S Paul, Ozren Polasek, Barbara Thorand, Harald Grallert, Michael Roden, Urmo Vosa, Tonu Esko, Caroline Hayward, Asa Johansson, Ulf Gyllensten, Nicholas Powell, Oskar Hansson, Niklas Mattsson-Carlgren, Peter K Joshi, John Danesh, Leonid Padyukov, Lars Klareskog, Mikael Landen, James F Wilson, Agneta Siegbahn, Lars Wallentin, Anders Malarstig, Adam S Butterworth, James E. Peters.
**Mapping pQTLs of circulating inflammatory proteins identifies drivers of immune-related disease risk and novel therapeutic targets**.
medRxiv 2023.03.24.23287680; doi: <https://doi.org/10.1101/2023.03.24.23287680>.

[^phenoscanner]: Kamat, M.A. et al. PhenoScanner V2: an expanded tool for searching human genotype-phenotype associations. Bioinformatics 35, 4851-4853 (2019).

[^OPG]: Kwan, J.S. et al. Meta-analysis of genome-wide association studies identifies two loci associated with circulating osteoprotegerin levels. Hum Mol Genet 23, 6684-6693 (2014).
