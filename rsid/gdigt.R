library(MRInstruments)
library(pQTLtools)

options(width=200)
dim(drug_interactions)
dim(gene_trials)

inf <- subset(inf1,prot!="BDNF")
inf_genes <- with(inf,gene)

gdi_inf <- subset(drug_interactions,entrez_gene_symbol%in%inf_genes)
gt_inf <- subset(gene_trials,entrez_gene_symbol%in%inf_genes)

dim(gdi_inf)
dim(gt_inf)
