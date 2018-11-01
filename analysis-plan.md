The cvd1 plan is quoted here since it largely overlaps with INFlammation analysis.

## cvd1 analysis plan

## Olink protein assays

The aim is to keep as many proteins as possible for the meta-analysis and we therefore suggest to dichotomize proteins with more than 20% missing values. The expectation is that these proteins will largely be the same across cohorts. Please include all proteins for which you have usable data, and further filtering can be done at the meta-analysis stage. At a later stage in the collaboration, imputation of missing values will be discussed.

*	Use standard linear regression for assays with 80 % of samples above the lower detection limit. 
*	Dichotomize proteins with more than 20 % of samples below the lower detection limit and code values below the detection limit as 0 and those above as 1. 
*	Rank-based inverse normal transformation, e.g., invnormal function from https://github.com/jinghuazhao/R/tree/master/gap

## List of proteins

## SNPs

*	1000 genomes imputation, any version 
*	SNPs will be filtered for imputation quality at time of meta-analysis, but please filter out SNPs with IMPUTE INFO quality less than 0.2
*	Standard QC, including call rate <95 % or failed Illumina genotype calling, gender mismatch, abnormal inbreeding coefficient, failed cryptic relatedness test, ancestry outlier, sample call rate <95 %, Bonferroni corrected Hardy-Weinberg Equilibrium test.

## Association analysis

*	Linear regression with adjustment for study-specific covariates. These should always include age at time of sample collection, gender and adjustment for population structure / geography if applicable (e.g across countries). Sample storage time and season of collection if applicable. 
*	Use imputation-dosages
*	Additive genetic model
*	Separate the analyses for men and women for X chromosome SNPs (exception for cohorts that have already performed analyses)

## Stratification

*	Analyse patients and controls separately –results will be merged at meta-analysis stage

6.	Descriptive statistics

Please fill out the attached descriptive statistics spreadsheet and use the naming convention: 

*	STUDY.descriptives.DATE.xls
*	Where, STUDY is a short (14 characters or less) identifier for the population studied, which is the same for all files provided by your study.
*	DATE is the date on which the file was prepared, in the format “YYYYMMDD”.

7.	GWAS results submission and file formats

SNP table for association results. Please include the following columns. Missing values are coded as “NA”.

V# | Variable name | Description
---|---------------|------------
V1 | SNPID | SNP ID as rs number
V2 | CHR | Chromosome number (1-22)
V3 | POS | Physical position for the reference sequence (please indicate NCBI build in descriptive file)
V4 | STRAND | Indicator of strand direction. Please specify “+” if positive or forward strand and “-” if negative or reverse strand. 
V5 | N | Number of non-missing observations
V6 | EFFECT_ALLELE | Allele for which the effect (beta coefficient) is reported. For example, in an A/G SNP in which AA = 0, AG=1, and GG=2, the coded allele is G.
V7 | REFERENCE_ALLELE | Second allele at the SNP (the other allele). In the example above, the non-coded allele is A. 
V8 | CODE_ALL_FQ | Allele frequency for the coded allele – “NA” if not available
V9 | BETA | Effect size for the coded allele, beta estimate from the genotype-phenotype association, with at least 5 decimal places. Note: if not available, please report “NA” for this variable.
V10 | SE | Standard error of the beta estimate, to at least 5 decimal places - “NA” if not available. 
V11 | PVAL | p-value of Wald test statistic – “NA” if not available
V12 | RSQ | Residual phenotypic variance explained by SNP. “NA” if not available
V13 | RSQ_IMP| Observed divided by expected variance for imputed allele dosage.
V14 | IMP | Please specify whether the SNP was imputed or genotyped: 1: imputed SNP, 0: directly genotyped SNP

8.	Meta-analysis

Meta-analysis will be performed using the inverse-N weighted analysis of p-values, as implemented in the software METAL (www.sph.umich.edu/csg/abecasis/metal). METAL implements a weighted Z-score method using the following formula: where the weight wi = square root of the sample size of the ith study, zi= -1(1-(pi/2))*(effect direction for study i), and pi is the P-value for the ith study.  
 
We will apply genomic control and the appropriate marker filters at this stage (i.e. please provide unfiltered results). 

a.	Marker exclusion filters: we will apply imputation quality filters at the meta-analysis stage. Please do not apply these filters yourself and provide unfiltered results. 
b.	Genomic control (GC): genomic control will be applied to each study at the meta-analysis stage (single GC). Please do not apply GC to GWAS results and provide uncorrected standard errors, as (double) GC will be applied at the meta-analysis stage. 
c.	Significance: the threshold for the genome-wide analyses will be set at 5 x 10-8. The results will be replicated in independent cohorts so no need for additional correction.
