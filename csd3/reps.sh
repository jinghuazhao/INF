# 21-1-2020 JHZ

export TMPDIR=/rds/user/jhz22/hpc-work/work

# https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
R --no-save <<END
  library(rentrez)
  library(XML)
# set_entrez_key("")
  Sys.getenv("ENTREZ_KEY")
  entrez_dbs()
  entrez_db_links("pubmed")
  pubmed_fields <- entrez_db_searchable("pubmed")
  r <- entrez_search(db="pubmed",term="(pQTLs OR (protein quantitative trait locus)) AND 2008:2020[DP] AND homo_sapiens[ORGN]",use_history=TRUE)
  class(r)
  names(r)
  with(r,web_history)
  unlink(paste("pubmed",c("fetch","summary"),sep="."))
  for(i in seq(1,with(r,count),50))
  {
    cat(i+49, "records downloaded\r")
    f <- entrez_fetch(db="pubmed", web_history=with(r,web_history), rettype="text", retmax=50, retstart=i)
    write.table(f, col.names=FALSE, row.names=FALSE, file="pubmed.fetch", append=TRUE)
    s <- entrez_summary(db="pubmed", web_history=with(r,web_history), rettype="text", retmax=50, retstart=i)
    fields <- c("uid", "pubdate", "sortfirstauthor", "title", "source", "volume", "pages")
    e <- extract_from_esummary(s, fields)
    write.table(t(e), col.names=FALSE, row.names=FALSE, file="pubmed.summary", append=TRUE, sep="\t")
  }
  link.example <- function(id=600807)
  {
    upload <- entrez_post(db="omim", id=id)
    asthma_variants <- entrez_link(dbfrom="omim", db="clinvar", cmd="neighbor_history", web_history=upload)
    asthma_variants
    snp_links <- entrez_link(dbfrom="clinvar", db="snp",
                             web_history=asthma_variants$web_histories$omim_clinvar,
                             cmd="neighbor_history")
    all_links <- entrez_link(dbfrom='pubmed', id=id, db='all')
  }
END

export traitmap=gwas_catalog_trait-mappings_r2020-01-16.tsv
export assoc=gwas_catalog_v1.0.2-associations_e99_r2020-01-16.tsv

cut -f2,8,11-13,18,22,25,26,28,29,35,36 annotate/$assoc > annotate/assoc.txt

join -113 -25 -t$'\t' \
    <(sed '1d' annotate/assoc.txt | sort -t$'\t' -k13,13) \
    <(sed '1d' annotate/$traitmap | grep 'protein measure' | sort -t$'\t' -k5,5) > annotate/ll

join <(sort -k1,1 work/pubmed.summary | awk '!/mice|Mice|plant|Plant|rice|soybean|Soybean|tomato/') \
     <(sed '1d' annotate/assoc.txt | sort -k1,1) > work/pubmed.left
