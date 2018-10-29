# Goal: create master mapping file for Olink proteins

rm(list=ls())

#--------------------- Import ---------------------#

library(readxl)
library(stringr)

# no of tabs in the shitty olink sheet
n=12

mylist <- vector('list', n)
names(mylist) <- c('cardiometabolic', 'cell.regul', 'cvd2', 'cvd3', 'devel', 'imm.resp',
                   'imm.onc', 'inf', 'metab', 'neu', 'onc2', 'organ.damage')

for (i in 1:n){

  mylist[[i]] <- as.data.frame(
                  read_xlsx(path= "~/post-doc/o-link/all-panels/Olink-validation-data-all-panels.xlsx", 
                            sheet = i, range='A1:B93', col_names = TRUE)
                              )
  
}

#--------------------- Reformat ---------------------#

col.rename <- function(x) {
  colnames(x) <- c("target", "uniprot")
  x
}

mylist <- lapply(mylist, col.rename) 

for (i in 1:n){
  mylist[[i]] <- data.frame(
                            mylist[[i]], panel=rep(names(mylist)[i], nrow(mylist[[i]]))
                            )
}

df <- do.call('rbind', mylist)
  
#--------------------- Clean-up ---------------------#

# data input error by Olink 'o' instead of 'O'
df$target <- gsub("IL-2oRA", "IL-20RA", df$target)

# turns out TWEAK labelled with an out of date or inferior UP id
df$uniprot[grep('TWEAK', df$target)] <- "O43508"

# Clean up 1. identify bad entries: Olink have made some errors and Excel import causes some problems

# clean whitespace

df$uniprot <- gsub('\\\r\\\n', ";", df$uniprot, ignore.case = F)

df$uniprot <- gsub(', |,', ";", df$uniprot, ignore.case = F)

df$uniprot <- gsub("[[:space:]]", "", df$uniprot)

# Clean up 2. '-' represents isoform notation eg O43521-2

df$uniprot.isoform <- NA

df$uniprot.isoform[grep('-', df$uniprot)] <- grep('-', df$uniprot, v=T)

df$uniprot <- gsub("-[0-9]$", "", df$uniprot)

# Special circumstances 1: uniprot is 'NA' 

na.ind <- grep('^NA', df$uniprot)

# turns out this is NTproBNP: give it the uniprot of BNP
for (i in na.ind){
  if (grepl('brain\ natriuretic\ peptide', df$target[i])){
    df$uniprot[i] <- 'P16860'
  }
}


# Special circumstances 2:two ids for protein complex eg IL12A-IL12B
# uniprot ids sep by ';'
# df[grep(";", df$uniprot), ]

df$multiple.proteins <- FALSE
df$multiple.proteins[grep(";", df$uniprot)] <- TRUE

df$protein.1 <- df$uniprot  
df$protein.2 <- NA

df$protein.2[which(df$multiple.proteins==T)] <- str_extract(string=df$uniprot, pattern = ";[A-Z0-9]+")[which(df$multiple.proteins==T)]
df$protein.2 <- gsub("^;", "", df$protein.2)

df$protein.1[which(df$multiple.proteins==T)] <- str_extract(string=df$uniprot, pattern = "^[A-Z0-9]+;")[which(df$multiple.proteins==T)]
df$protein.1 <- gsub(";$", "", df$protein.1)

# where there are 2 uniprot ids (eg protein complex) the uniprot ids are not always in consistent order
# lets make them in consistent alphabetical order
df$uniprot.ordered <- NA
df$uniprot.ordered[!df$multiple.proteins] <- df$uniprot[!df$multiple.proteins]

alphabetize.up <- function(x) {
  if( !inherits(x, what='data.frame')){
    stop('argument must be a data.frame')
    }
  y <- paste(sort(c( x$protein.1, x$protein.2)),collapse=';')
  y
  }

inds <- which(df$multiple.proteins)

for (i in inds){
  df$uniprot.ordered[i] <- alphabetize.up(df[i,]) 
}

#annoying that p1 and p2 are arbitrary: now we've ordered things let's start over on this front
df$protein.1 <- df$protein.2 <-NULL

# now repeat the exercise for p1 and p2 using the alphabetized concatenation

df$protein.1 <- df$uniprot.ordered  
df$protein.2 <- NA

df$protein.2[which(df$multiple.proteins==T)] <- str_extract(string=df$uniprot.ordered, pattern = ";[A-Z0-9]+")[which(df$multiple.proteins==T)]
df$protein.2 <- gsub("^;", "", df$protein.2)

df$protein.1[which(df$multiple.proteins==T)] <- str_extract(string=df$uniprot.ordered, pattern = "^[A-Z0-9]+;")[which(df$multiple.proteins==T)]
df$protein.1 <- gsub(";$", "", df$protein.1)

# col to identify dup proteins and which panels

dup.prots <- union( which( duplicated(df$uniprot.ordered)), which( duplicated(df$uniprot.ordered, fromLast = T)) )
df$prot.on.multiple.panel <- FALSE
df$prot.on.multiple.panel[dup.prots] <- TRUE

df$panels.with.prot <- NA

  
tmp.list <- split( df[dup.prots,], f=df$uniprot.ordered[dup.prots] )

mylist <- lapply(tmp.list, FUN = function(x) paste( as.character(x$panel), collapse=";" ) )

for (i in dup.prots){
  uprot <- df$uniprot.ordered[i]
  df[i, "panels.with.prot"] <- mylist[[uprot]]
}

#--------------------- Gene symbol annotation ---------------------#

# matching to gene symbols: do for p1 and p2

library(biomaRt)

#ensembl <- useMart(biomart="ensembl",
#                   dataset="hsapiens_gene_ensembl",
#                   host='http://jul2018.archive.ensembl.org')

#filters <- listFilters(ensembl)

#x <- getBM(attributes = c('uniprotswissprot', 'hgnc_symbol', 'entrezgene', 'chromosome_name'), 
#             filters = 'uniprotswissprot', 
#             values = df$protein.1, 
#             mart = ensembl)

# some UP ids not found by BioMart: turns out we have outdated IDs 
#df[which(!df$protein.1 %in% x$uniprotswissprot),]

#--------------------- Try an archived version of Ensembl ---------------------#

# find urls for old ensembl versions
listEnsemblArchives()

# hg19/GRCh37
ensembl.hg19 <- useMart(biomart= "ENSEMBL_MART_ENSEMBL",
                     dataset="hsapiens_gene_ensembl",
                     host = 'http://grch37.ensembl.org')

# note attribute names differ in the older release
gene.pos <- getBM(attributes = c('uniprotswissprot', 'hgnc_symbol', # 'entrezgene',
                          'chromosome_name', 'start_position', 'end_position'), 
           filters = 'uniprotswissprot', 
           values = unique(df$protein.1), 
           mart = ensembl.hg19)

# there are some duplicated genes

dup.ind <- union( which(duplicated(gene.pos$hgnc_symbol)),
                  which(duplicated(gene.pos$hgnc_symbol, fromLast = T))
              )

# strange chr names

strange.ind <- which(!gene.pos$chromosome_name %in% c(1:22, 'X', 'Y'))

to.cut <- intersect(dup.ind, strange.ind)

gene.pos2 <- gene.pos[-to.cut,]

#-------------------------------------------------------------------------------#

df2 <- merge(x = df, y = gene.pos2, by.x = "protein.1", by.y = "uniprotswissprot", all = TRUE)

#-------------------------------------------------------------------------------#

inf <- df2[df2$panel=="inf",]

inf  <- data.frame(target.short= gsub("^.+\\(|)$", "", inf$target), inf)

# load the olink inf eset

library(Biobase)

eset <- readRDS("~/post-doc/o-link/esets/round2/post-qc/eset.inf1.flag.out.outlier.out.rds")

features <- fData(eset)

features <-  features[, c("olink.id", "uniprot.id")]

features$common.name <- gsub("^[0-9]+_", "", features$olink.id)

features$uniprot.id <- gsub("[[:space:]]", "", features$uniprot.id)

if (all(features$uniprot.id %in% inf$uniprot)==FALSE){
  warning("not all the uniprot ids in the data olink suppliedare in the annotation file")
  
  tmp <- features[!features$uniprot.id %in% inf$uniprot, ] 
  
}

inf.2 <- merge(x = inf, y = features[,c("olink.id", "uniprot.id")], by.x = "protein.1", by.y = "uniprot.id", all = TRUE)

# which proteins have not been mapped to a gene

non.mapped <- which(is.na(inf.2$hgnc_symbol))
inf.2[non.mapped,]

# turns out there are alternative up ids for these, which explains why biomart couldn't map them to gene symbols

inf.2$alternate.uniprot <- NA

inf.2$alternate.uniprot[which(inf.2$target.short=="FGF-5")] <- "P12034"
inf.2$alternate.uniprot[which(inf.2$target.short=="CD6")] <- "P30203"


posn <- getBM(attributes = c('uniprotswissprot', 'hgnc_symbol', # 'entrezgene',
                                 'chromosome_name', 'start_position', 'end_position'), 
                  filters = 'uniprotswissprot', 
                  values = inf.2$alternate.uniprot, 
                  mart = ensembl.hg19)

for (i in non.mapped){
  up <- inf.2$alternate.uniprot[i]
  inf.2$chromosome_name[i] <- posn$chromosome_name[which(posn$uniprotswissprot==up)]
  inf.2$start_position[i] <- posn$start_position[which(posn$uniprotswissprot==up)]
  inf.2$end_position[i] <- posn$end_position[which(posn$uniprotswissprot==up)]
}

# check it worked
inf.2[non.mapped,]

inf.final <- inf.2[ ,c("target", "target.short", "uniprot", "panel", "prot.on.multiple.panel",
                       "panels.with.prot",  "hgnc_symbol", "chromosome_name", "start_position", "end_position", "olink.id", "alternate.uniprot")]

write.table(inf.final, file="~/post-doc/o-link/scallop/olink.inf.panel.annot.txt",
            row.names=F, col.names =T, sep="\t")