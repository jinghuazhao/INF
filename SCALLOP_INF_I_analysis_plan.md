# SCAndinavian coLLaboration for Olink plasma Protein genetics – INF panel proteins

## 1. Background

The SCALLOP consortium was created to work collaboratively on discovery and follow-up of pQTLs for proteins measured using Olink technology. A meta-analysis has been conducted on data from participating cohort Olink CVD panel and consequent contributions have been made on Olink INF panel. This document therefore follows closely the analysis plan for the analysis while highlighting relevant information which will facilitate the meta-analysis

## 2. Aims

As with the CVD I meta-analysis, the tasks will include

* Identification of pQTLs in SCALLOP discovery cohorts
* Replication of pQTLs in SCALLOP replication cohorts
* Investigation of the mechanistic basis of identified cis- and trans-pQTL by functional annotation
* Investigation of pleiotropic effects of the pQTLs
* Evaluate whether the CVD I proteins are causal in e.g. CHD and stroke

## 3. Data analysis

* Use standard linear regression for assays with 80% of samples above the lower detection limit. 
* Dichotomize proteins with more than 20% of samples below the lower detection limit and code values below the detection limit as 0 and those above as 1. 
* Rank-based inverse normal transformation, e.g., invnormal function from https://github.com/jinghuazhao/R/tree/master/gap

### SNPs

* 1000 genomes imputation
* SNPs will be filtered for imputation quality at time of meta-analysis, but please filter out SNPs with IMPUTE INFO quality less than 0.2
* Standard QC, including call rate < 95% or failed Illumina genotype calling, gender mismatch, abnormal inbreeding coefficient, failed cryptic relatedness test, ancestry outlier, sample call rate < 95%, Bonferroni corrected Hardy-Weinberg Equilibrium test.

### Association analysis

* Linear regression with adjustment for study-specific covariates. These should always include age at time of sample collection, gender and adjustment for population structure / geography if applicable (e.g across countries). Sample storage time and season of collection if applicable. 
* Use imputation-dosages
* Additive genetic model
* Separate the analyses for men and women for X chromosome SNPs (exception for cohorts that have already performed analyses)

### Stratification

* Analyse patients and controls separately –results will be merged at meta-analysis stage

## 4. Descriptive statistics

Please fill out the attached descriptive statistics spreadsheet and use the naming convention: 

* STUDY.descriptives.DATE.xls
* Where, STUDY is a short (14 characters or less) identifier for the population studied, which is the same for all files provided by your study.
* DATE is the date on which the file was prepared, in the format “YYYYMMDD”.

## 5. GWAS results submission and file formats

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

Note that due possibly to the large number of proteins for GWAS, some cohorts employed PLINK to expedite analysis in which case one sees the following information: 

No | Name | Description | Comment
--|----------|--------------------
1 | BP | Position in basepairs	
2 | CHR | Chromosome	
3 | SNP | SNP name/chr:pos_a1_a2	
4 | HWE | Hardy-Weinberg equilibrium P	
5 | MAF | Minor allele frequency	
6 | A1 | Allele 1	
7 | A2 | Allele 2	
8 | N |	Sample size	
9 | BETA | Regression coefficient	
10 | CHI2 | Regression statistics	
11 | P | P value	

In this case, if is preferable to provide strand, effect allele, effect allele frequency, and the information measures.

## 6. Meta-analysis

Meta-analysis will be performed using the inverse-N weighted analysis of p-values, as implemented in the software METAL (www.sph.umich.edu/csg/abecasis/metal). METAL implements a weighted Z-score method using the following formula: where the weight wi = square root of the sample size of the ith study, zi= -1(1-(pi/2))*(effect direction for study i), and pi is the P-value for the ith study.  
 
We will apply genomic control and the appropriate marker filters at this stage (i.e. please provide unfiltered results). 

*. Marker exclusion filters: we will apply imputation quality filters at the meta-analysis stage. Please do not apply these filters yourself and provide unfiltered results. 
*. Genomic control (GC): genomic control will be applied to each study at the meta-analysis stage (single GC). Please do not apply GC to GWAS results and provide uncorrected standard errors, as (double) GC will be applied at the meta-analysis stage. 
*. Significance: the threshold for the genome-wide analyses will be set at 5 x 10<sup>-8</sup>. The results will be replicated in independent cohorts so no need for additional correction.

## 7. Uploading of results data to TRYGGVE server

See CVD I analysis plan.

## 8. SCALLOP consortium

### Information is available from https://www.olink.com/scallop/.
