## Predictive modeling

Given the high level of proportions of variance explained here and empirical data{Ganz, 2016 #4}{Williams, 2019 #15}{Brandes, 2020 #234}. it is 
appropriate to construct scores on significant variants (polygenic score, PGS) or genomewide variants based on total number of changes (polygenic 
risk score, PRS). The most straightforward approach is to consider effects of individual proteins through implementations such as PRSice{Choi, 2020 
#238}, LDpred{Vilhjálmsson, 2015 #203}{Privé, 2020 #242} and lassosum{Mak, 2017 #202} and use them as predictors for disease outcomes (e.g., through 
UKBB data). With individual level data, one can attempt to construct scores on multiple proteins{Benson Mark, 2018 #38} jointly with genotype data. 
In general, let y be our trait of interest, X be a matrix of M protein abundances measured in N individuals. The contribution of proteins to trait 
as a whole is captured through a linear predictor=Xwhererepresents the regression coefficients or protein scores{Ganz, 2016 #4}{Mosley, 2018 #168} 
in a generalized linear model (GLM) framework such that E(y)==g-1(), Var(y)=Var(), where g(.) is a link function used to accommodate various types 
of outcomes such as continuous (e.g., FEV1) and binary outcomes (e.g. CAD). In this case unlike in the usual GLM the proteins themselves are further 
modeled by relevant pQTLs, i.e., Xi=jijZij, with Zij being the dosage of the genetic variants for the i-th protein and j-th variant as in PGS or 
PRS. Note that ij is the associated weight involving a variable number of variants across proteins. This is reminiscent of a nested model in path 
analysis and the model only using individual protein scores can be considered marginal. One may also wonder how effect sizes were mirrored in an 
established biomarker such as C-reactive protein (CRP). Effect sizes of the 162 unique variants from the meta-analysis for invnormal(CRP) in the 
UKBB (http://www.nealelab.is/blog/2019/9/16/biomarkers-gwas-results) showed small correlation (r=0.049) with those in the meta-analysis; there was 
weak correlation (r=-0.22, p=0.27) from those variants from invnormal(CRP) with CVD but it was not optimal compared to those based on PRSice (data 
not shown).

