#!/usr/bin/bash

module load gcc/6
R --no-save <<END
  library(XML)
  library(gap)
  library(rentrez)
  term <- "pQTLs OR (protein AND quantitative AND trait AND loci) AND human[MH] AND plasma"
  r <- entrez_search(db="pubmed",term=term,retmax=3000)
  write.table(with(r,ids),file="entrez.ids",col.names=FALSE,row.names=FALSE,quote=FALSE)
# Yao 30111768[uid]
  xlsx <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-05512-x/MediaObjects/41467_2018_5512_MOESM1_ESM.xlsx"
  st1 <- openxlsx::read.xlsx(xlsx,colNames=TRUE, startRow=3)
  xlsx <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-05512-x/MediaObjects/41467_2018_5512_MOESM7_ESM.xlsx"
  st6 <- openxlsx::read.xlsx(xlsx,colNames=TRUE, startRow=3)
  st61 <- merge(st6,st1[,1:4],by.x="Protein",by.y="Protein.Abbreviation")
  yao <- merge(gap::inf1,st61,by.x="gene",by.y="Protein-coding.Gene.Abbreviation")
# Solomon 30562114[uid]
# Sun 29875488[uid]
  xlsx <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0175-2/MediaObjects/41586_2018_175_MOESM4_ESM.xlsx"
  st4 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, startRow=5)
  st19 <- openxlsx::read.xlsx(xlsx, sheet=19, colNames=TRUE, startRow=2)
  st20 <- openxlsx::read.xlsx(xlsx, sheet=20, colNames=TRUE, startRow=3)
  knownlist <- c(with(st19,PMID), "30111768", "29875488", "30562114")
  revlist <- replace(knownlist,knownlist=="PMC4358658",25652787)
  f <- entrez_fetch(db="pubmed",id=revlist,rettype="text",retmax=3000)
  s <- entrez_summary(db="pubmed",id=revlist,rettype="text",retmax=3000)
  fields <- c("uid", "pubdate", "sortfirstauthor", "title", "source", "volume", "pages")
  e <- extract_from_esummary(s, fields)
  write.table(t(e), col.names=FALSE, row.names=FALSE, file="pubmed.list", sep="\t")
  r1 <- merge(inf1,st20,by.x="uniprot",by.y="UniProt")
  r2 <- merge(inf1,st4,by.x="uniprot",by.y="UniProt")
  intersect(with(r1,uniprot),with(r2,uniprot))
END

split --lines=250 --numeric-suffixes=1 --suffix-length=1 entrez.ids
cat x1 | tr '\n' ' ' | xsel -i

cd annotate
if [ ! -f PMC-ids.csv.gz ]; then wget ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/PMC-ids.csv.gz; fi
zgrep -e PMC4698720 -e PMC4358658 PMC-ids.csv.gz
export assoc=gwas_catalog_v1.0.2-associations_e100_r2020-08-26.tsv
export traitmap=gwas_catalog_trait-mappings_r2020-08-26.tsv
cut -f2,8,11-13,18,22,25,26,28,29,35,36 $assoc > assoc.txt
join -113 -25 -t$'\t' \
    <(sed '1d' assoc.txt | sort -t$'\t' -k13,13) \
    <(sed '1d' $traitmap | awk '/protein measure|protein levels/' | sort -t$'\t' -k5,5) | cut -f2 | uniq > pmid
cd -

# https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
# https://pubmed.ncbi.nlm.nih.gov/
# https://www.ncbi.nlm.nih.gov/pmc/pmctopmid/
# https://www.ebi.ac.uk/gwas/api/search/downloads/alternative
# https://www.ebi.ac.uk/gwas/api/search/downloads/trait_mappings
