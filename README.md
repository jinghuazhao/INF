# SCALLOP-INF meta-analysis

## Flow of analysis

(view diagram by pasting script to [mermaid live editor](https://mermaid-js.github.io/mermaid-live-editor/))

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

Further details is available from [https://jinghuazhao.github.io/INF/](https://jinghuazhao.github.io/INF/).
