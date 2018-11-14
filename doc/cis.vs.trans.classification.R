# "Thu Nov  8 12:13:07 2018"

# author jp549@cam.ac.uk

# identify cis vs trans hits

# rule: a cis acting variant lies within the region
# from 1MB upstream of the start position to 1MB downstream of the end position 
# of the gene that encodes the protein being tested in the GWAS

# All signals that are outside this window will be defined as trans


rm(list=ls())

options(stringsAsFactors = F)

# Take Jing's results file "/scratch/jhz22/INF/work/INTERVAL.jma.dat"
results <- read.table("/scratch/jhz22/INF/work/INTERVAL.jma.dat", head=T)

# get mapping to uniprot ids
ids <- read.table("/scratch/jhz22/INF/inf1.list", head=F)

colnames(ids) <- c('casual.name', 'UPid')

# map on the uniprot ids
res2 <- merge(x= results, y= ids, by.x='prot', by.y='casual.name', all.x=T)

# Take Jimmy's mapping file (on my laptop, ~/post-doc/o-link/scallop/olink.inf.panel.annot.txt)
# proteins <- read.table("/scratch/jp549/olink.inf.panel.annot.txt", head=T)
proteins <- read.table("/scratch/jhz22/INF/doc/olink.inf.panel.annot.tsv", head=T)
proteins[with(proteins,uniprot=="Q8NF90"),"hgnc_symbol"] <- "FGF5"
proteins[with(proteins,uniprot=="Q8WWJ7"),"hgnc_symbol"] <- "CD6"

# keep the relevant columns
col.keep <- c("target", "target.short", "uniprot", "hgnc_symbol", 
              "chromosome_name", "start_position", "end_position")

proteins <- proteins[ ,col.keep]

# add a prefix 'p.' so we know these cols refer to the protein being GWAS'd
colnames(proteins) <- paste("p", colnames(proteins), sep=".")

# map on to the results file, using uniprot as the common reference
res3 <- merge(x= res2, y= proteins, by.x='UPid', by.y='p.uniprot', all.x=T)

# classify into cis and trans

# set cis as  -1MB upstream to +1MB downstream
radius <- 1e+6

res3$cis.start <- res3$p.start_position - radius
if (any(res3$cis.start <0 )){
  res3$cis.start[which(res3$cis.start<0)] <- 0
}

res3$cis.end <- res3$p.end_position + radius

res3$cis <- rep(NA, nrow(res3))

# any variant on a different chromosome to the gene encoding the target protein is not cis
dist.inds <- which(res3$Chr != res3$p.chromosome_name)

if (length(dist.inds)>0){
  res3$cis[dist.inds] <- FALSE
}

# for ones on the same chr, we can't be sure without looking at position

same.inds <- which(res3$Chr == res3$p.chromosome_name)

if (length(same.inds)>0){
  # see if variant lies in the cis region
  res3$cis[same.inds] <- res3$bp[same.inds] > res3$cis.start[same.inds]  & res3$bp[same.inds] < res3$cis.end[same.inds]
}

res3$cis.trans <- rep(NA, nrow(res3))

res3$cis.trans[res3$cis==T] <- "cis"
res3$cis.trans[res3$cis==F] <- "trans"

# split by protein
list.by.prot <- split(res3, f= res3$p.hgnc_symbol)

# get the breakdown of cis vs trans per protein
# sapply(list.by.prot, function(x) table(x$cis.trans) )

cis.trans.per.prot <- with(res3,table(p.hgnc_symbol, cis.trans))
cis.trans.per.prot
sum(as.matrix(cis.trans.per.prot))
