#!/usr/bin/bash

R --no-save <<END
  library(rentrez)
  r <- entrez_search(db="pubmed",term="(pQTL OR pQTLs OR (protein AND quantitative AND trait AND locus)) AND 2008:2020[DP] AND homo_sapiens[ORGN]",
       retmax=3000)
  write.table(with(r,ids),file="rentrez.ids",col.names=FALSE,row.names=FALSE,quote=FALSE)
END
split --lines=250 --numeric-suffixes --suffix-length=1 rentrez.ids
cat x0 | tr '\n' ' ' | xsel -i
cat x1 | tr '\n' ' ' | xsel -i
cat x2 | tr '\n' ' ' | xsel -i
cat x3 | tr '\n' ' ' | xsel -i
cat x4 | tr '\n' ' ' | xsel -i

export TMPDIR=$INF/work

R --no-save <<END
  source("rsid/rentrez.ini")
  # Yao 30111768[uid]
  st1 <- openxlsx::read.xlsx("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-05512-x/MediaObjects/41467_2018_5512_MOESM1_ESM.xlsx",
         colNames=TRUE, startRow=3)
  st6 <- openxlsx::read.xlsx("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-05512-x/MediaObjects/41467_2018_5512_MOESM7_ESM.xlsx",
         colNames=TRUE, startRow=3)
  st61 <- merge(st6,st1[,1:4],by.x="Protein",by.y="Protein.Abbreviation")
  yao <- merge(gap::inf1,st61,by.x="gene",by.y="Protein-coding.Gene.Abbreviation")
  # Solomon 30562114[uid]

  # Sun 29875488[uid]
  xlsx <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0175-2/MediaObjects/41586_2018_175_MOESM4_ESM.xlsx"
  st4 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, startRow=5)
  st19 <- openxlsx::read.xlsx(xlsx, sheet=19, colNames=TRUE, startRow=2)
  st20 <- openxlsx::read.xlsx(xlsx, sheet=20, colNames=TRUE, startRow=2)
  knownlist <- c(with(st19,PMID), "30111768", "29875488", "30562114")
  revlist <- replace(knownlist,knownlist=="PMC4358658",25652787)
  fetch_knownlist(revlist)
  r1 <- merge(inf1,st20,by.x="uniprot",by.y="UniProt")
  r2 <- merge(inf1,st4,by.x="uniprot",by.y="UniProt")
  intersect(with(r1,uniprot),with(r2,uniprot))
# https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
# Add a couple of recent papers
# output list by uniprot ID
END
# LD with INF1.merge

# Additional search

cd annotate
# https://www.ncbi.nlm.nih.gov/pmc/pmctopmid/
if [ ! -f PMC-ids.csv.gz ]; then wget ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/PMC-ids.csv.gz; fi
zgrep -e PMC4698720 -e PMC4358658 PMC-ids.csv.gz
# https://www.ebi.ac.uk/gwas/api/search/downloads/alternative
export assoc=gwas_catalog_v1.0.2-associations_e99_r2020-01-16.tsv
# https://www.ebi.ac.uk/gwas/api/search/downloads/trait_mappings
export traitmap=gwas_catalog_trait-mappings_r2020-01-16.tsv
cut -f2,8,11-13,18,22,25,26,28,29,35,36 $assoc > assoc.txt
join -113 -25 -t$'\t' \
    <(sed '1d' assoc.txt | sort -t$'\t' -k13,13) \
    <(sed '1d' $traitmap | awk '/protein measure|protein levels/' | sort -t$'\t' -k5,5) | \
cut -f2 | uniq > pmid
join <(awk '!/mice|Mice|plant|Plant|rice|soybean|Soybean|tomato/' pubmed.summary | sort -k1,1) \
     <(sed '1d' assoc.txt | sort -k1,1) > pubmed.left
cut -d' ' -f1 pubmed.left | uniq | wc -l
cd -
