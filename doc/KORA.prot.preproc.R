# Script started "Thu Oct  4 11:09:34 2018"

# To run on TRYGYVVE

# Purpose: take INF1 proteomics data from KORA cohort, and normalise for subsequent linear regression/GWAS

# Oddities: some individuals without ANY proteomic data are included in the .csv file

rm(list=ls())

setwd("/data/andmala/KORA/")

# read in the Olink INF proteomic data
x <- read.csv("K15117g_Malarstig_tra_20180502.csv", head=T)

# initial columns have id and covariate data. Protein columns have prefix "uh_o"
# Find the indices of these columns
prot.ind <- grep("^uh_o_", colnames(x) )

# subset to get just the protein data
pr <- x[ ,prot.ind]

if (ncol(pr) != 92){
  stop("Incorrect number of proteins: there should be 92")
}

total.na <- sum( apply(pr, 1, FUN= function(x) all(is.na(x)) ) )

message(paste(total.na, "samples have no protein values whatsoever"))

# indices of samples with no Olink data
total.na.ind <- which( apply(pr, 1, FUN= function(x) all(is.na(x)) ) )

# indices of samples with no Axiom id (ie not genotyped)
axiom.na <- which(is.na(x$zz_nr_axiom))

to.cut <- union(total.na.ind, axiom.na)

# cut samples with no protein data
y <- x[-to.cut, ]

# let's investigate whether protein values are truncated at the lower end of the distribution:
# Olink set values below lower limit of detection (LLOD) to LLOD,
# although the untruncated data can be provided

find.n.llod <- function(x){
  # find the minimum value
  x.min <-x[which.min(x)]
  if (any(x < x.min)){
    stop("Error: vector contains value lower than supposed minimum")
  }
  n.llod <- sum(x == x.min)
  n.llod
}

# find the number of individuals at LLOD for each protein
n.llod.by.prot <- sort( apply(y[,prot.ind], 2, FUN= find.n.llod ), decreasing = T )

pc.llod.by.prot <- 100*n.llod.by.prot/nrow(y)

# set a % threshold for proportion individuals below LLOD, above which to plot proteins
thr <- 20

pdf("~/KORA.pc.below.llod.pdf")
barplot(height= pc.llod.by.prot[which(pc.llod.by.prot>thr)],
     names.arg=  toupper( gsub("^uh_o_", "", names(pc.llod.by.prot[which(pc.llod.by.prot>thr)])) ),
     ylim=c(0,100),
     xlab="",
     ylab= "% individuals below LLOD", 
     main= "KORA cohort: proteins with high % of individuals below LLOD",
     las=2)
dev.off()

# Define the inverse normalisation function
#https://github.com/antagomir/scripts/blob/master/R/inverse.normalization.R
inverse_normal_transform <- function (x) { qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))) }

# Apply invnt transformation to proteins 
# nb apply transposes result: proteins become cols
pr.transformed <- apply(y[,prot.ind], 2, FUN=inverse_normal_transform)

# 
if (identical( row.names(pr.transformed), row.names(y[,-prot.ind]))){
  z <- cbind(y[,-prot.ind], pr.transformed)
} else {
  stop("row names don't match")
}

colnames(z)[which(colnames(z)=="utalter")] <- "age"
colnames(z)[grep("sex", colnames(z))] <- "sex"

ind.keep <- grep("zz_nr_axiom|age|sex|^uh_o_", colnames(z))

# just keep ids, age, sex, protein levels

final <- z[,ind.keep]

write.table(final, file="~/kora.normalised.prot.txt",
            sep="\t", quote=F, row.names = F)