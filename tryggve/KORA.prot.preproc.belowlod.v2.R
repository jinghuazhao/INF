# Script modified "Thu Mar  7 16:52:15 2019"

# To run on TRYGYVVE

# Purpose: take INF1 proteomics data from KORA cohort (untruncated NPX values),
# and normalise for subsequent linear regression/GWAS

# Oddities: some individuals without ANY proteomic data are included in the .csv file

# Issues:

# 1) LIF is not in the original KORA dataset provided by Anders "K15117g_Malarstig_tra_20181129.csv": this was supplied later 
# as "K15117g_Malarstig_tra_20190219.csv" and needs to be merged in. There are also different nrow for these csv files

# 2) # NB 2 protein names have typos 'ORTN' and 'OT3' ('NRTN' and 'NT3').
# The file columns for proteins in this dataset generally begin "UH_O_" followed by protein symbol.
# In the case of these 2 proteins this formatting pattern has been messed up so the names are "UH_ORTN" and "UH_OT3".
# So it looks like the "_N" from the protein name got lost here i.e. they should be "UH_O_NRTN" and "UH_O_NT3".

#-----------------------------------------------# 
# 1) Data set-up 
#-----------------------------------------------# 

options(stringsAsFactors = F)

rm(list=ls())

# this should be the file with untruncated values (ie not censored at 'LLOD')
data.dir <- "/data/andmala/KORA/KORA_updated_olink/"
outdir <- "/data/jampet/KORA/"

setwd(data.dir)

# read in the Olink INF proteomic data
x <- read.csv("K15117g_Malarstig_tra_20181129.csv", head=T)

# important: Olink have withdrawn BDNF: the assay is dodgy
# see: https://www.olink.com/bdnf-info/

x$UH_O_BDNF <- NULL

# fix colnames
colnames(x) <- gsub("UH_ORTN", "UH_O_NRTN", colnames(x))
colnames(x) <- gsub("UH_OT3", "UH_O_NT3", colnames(x))

# Protein LIF was omitted from the first data upload
lif <- read.csv("/data/jinhua/data/KORA/K15117g_Malarstig_tra_20190219.csv")

x <- merge(x, lif, by="ZZ_nr", all.x=T)

# initial columns have id and covariate data. Protein columns have prefix "UH_O"
# Find the indices of these columns
prot.ind <- grep("^uh_o", colnames(x), ignore.case = TRUE )

# subset to get just the protein data
pr <- x[ ,prot.ind]

if (ncol(pr) != 91){
  stop(paste( "Incorrect number of proteins: there should be 91. There are", ncol(pr)) )
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
# originally we found the data was truncated when we expected it not to be

# here's a simple test: for each protein find the lowest value:
# if there are many ties then likely the data is censored

find.n.llod <- function(x){
  # find the minimum value
  x.min <- x[which.min(x)]
  if (any(x < x.min, na.rm = T)){
    stop("Error: vector contains value lower than supposed minimum")
  }
  n.llod <- sum(x == x.min, na.rm=T)
  n.llod
}

# find the number of individuals at the minimum value for each protein
n.llod.by.prot <- sort( apply(y[,prot.ind], 2, FUN= find.n.llod ), decreasing = T )

# see if the minimum value is not unique for any protein (suggesting truncation ...)
if (any( n.llod.by.prot != 1)){
  message("there are multiple minimum values for at least one protein")
  message("Checking for culprit proteins...")
  n.llod.by.prot[n.llod.by.prot != 1]
}

# ucsex:  M = 1, F= 2, then colour code by sex

y$ucsex[which(y$ucsex==1)] <- 'M'
y$ucsex[which(y$ucsex==2)] <- 'F'

#-----------------------------------------------# 
# 2) Visualisation 
#-----------------------------------------------# 

pdf(paste0(outdir,"kora.below.llod.npx.by.sex.",Sys.Date(),".pdf"), height=7, width = 6)

par(mfrow=c(4,3), mar=c(4,4,2,1))
for (i in prot.ind){
  boxplot( y[,i] ~ y$ucsex, las=1, col= c("blue", "pink"))
  mtext(text= colnames(y)[i], side=3)
}

dev.off()

pdf(paste0(outdir, "kora.below.llod.npx.hist.",Sys.Date(),".pdf"), height=7, width = 6)

par(mfrow=c(4,3), mar=c(4,4,2,1))
for (i in prot.ind){
  hist( y[,i], freq=F, breaks=50, main="", xlab='', ylab=NA, las=1 )
  
  # add a density curve
  curve(
    expr= dnorm(x, mean=mean(y[,i], na.rm=TRUE), sd=sd(y[,i], na.rm=T)),
    add=TRUE, col='darkblue', lwd=2
  )
  mtext(text= 'NPX expr', side=1, line=2, cex=0.75)
  mtext(text= colnames(y)[i], side=3)
}

dev.off()

pdf(paste0(outdir,"kora.below.llod.npx.distrib.",Sys.Date(),".pdf"), height=7, width = 6)

par(mfrow=c(4,3), mar=c(4,3,2,1))
for (i in prot.ind){
  plot(density( y[,i], na.rm=T), main="", ylab=NA, las=1 )
  mtext(text= colnames(y)[i], side=3)
}

dev.off()

#-----------------------------------------------# 
# 3) Normalisation 
#-----------------------------------------------# 

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

ind.keep <- grep("zz_nr_axiom|age|sex|^uh_o_", colnames(z), ignore.case = T)

# just keep ids, age, sex, protein levels

final <- z[,ind.keep]

#-----------------------------------------------# 
# 4) Write out 
#-----------------------------------------------# 

write.table(final, file=paste0(outdir,"kora.below.llod.normalised.prot.",Sys.Date(),".txt"),
            sep="\t", quote=F, row.names = F)