rm(list=ls())
library(RCircos)
#library(hugene10stensembl.db)
#library(SNPlocs.Hsapiens.dbSNP.20120608)
# get gene data, and format
wd <- "~/post-doc/serpina1/"
setwd(wd)
snp <- 'rs28929474'
P <- 5e-8
snp.data <- read.delim("~/post-doc/serpina1/rs28929474.txt", head=T, stringsAsFactors=F)
inds <- c('updated_rsid', 'chromosome', 'position', 'alleleA_meta', 'alleleB_meta', 'beta_meta', 'log_pvalue_meta', 'Direction', 'Somamer', 'cis_tran', 'B_allele_freq_meta')
snp.data <- snp.data[,inds]
passed.qc <- read.csv("~/post-doc/interval-cardio/Somamers_all_Passed_CV20.csv", head=T)
snp.data <- snp.data[which(snp.data$log_pvalue_meta > -log10(P)),]
qc.fails = snp.data$Somamer[!snp.data$Somamer %in% passed.qc$SOMAS_ID_round2]
# keep SERPINA1 even though QC fail so we can plot it
qc.fails = qc.fails[-grep("SERPINA", qc.fails)]
snp.data <- snp.data[-which(snp.data$Somamer %in% qc.fails), ]
rownames(snp.data) <- 1:nrow(snp.data)
sort( gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]$", "", snp.data$Somamer ) )
# get protein positions
pr.pos <- read.delim("~/post-doc/data-obj/SOMALOGIC_Master_Table_160410_1129info.tsv", head=T, stringsAsFactors= F)
pr.pos2 <- pr.pos[match( snp.data$Somamer, pr.pos$SOMAS_ID_round2, nomatch = 0), ]
if( identical(pr.pos2$SOMAS_ID_round2, snp.data$Somamer) ) snp.data.new <- cbind(snp.data,
    pr.pos2[,c("SOMAS_ID_round2", "EntrezGeneSymbol", "Target", "chromosome_name", "start_position", "end_position")])
g.data <- snp.data.new[,c("chromosome_name", "start_position", "end_position", "EntrezGeneSymbol", "log_pvalue_meta", "Direction", "beta_meta")]
g.data$chromosome_name <- paste0("chr", g.data$chromosome_name)
# get genome background
data(UCSC.HG19.Human.CytoBandIdeogram)
chr.exclude <- NULL # paste0('chr', 'Y')
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
tracks.inside <- 5
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info=cyto.info, chr.exclude, tracks.inside=tracks.inside, tracks.outside=tracks.outside)
# add link data (lines to snp)
#snpPos <- rsid2loc(snp)
#names(snpPos) <- gsub("^ch", "chr", names(snpPos))
#Link.Data <- data.frame(prbs[,1:3], Chromosome.1= rep(names(snpPos), nrow(prbs)),
#                        chromStart.1= rep(snpPos,  nrow(prbs)), chromEnd.1= rep(snpPos,  nrow(prbs)))
Link.Data <- data.frame(
  g.data[,1:3],
  Chromosome.1= paste0("chr", snp.data.new$chromosome),
  chromStart.1= snp.data.new$position,
  chromEnd.1= snp.data.new$position,
  dir=g.data$Direction,
  effectsize= g.data$beta_meta
)
# Plot
out.file <- paste0(wd, "Circos_", snp,"_liberal-v3.pdf")
pdf(file=out.file, height=8, width=8, compress=TRUE);
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size <- 0.5
RCircos.Reset.Plot.Parameters(rcircos.params)
#name.col <- 4 # 4th col in prbs ie gene name
#side <- "in"
#track.num <- 1
# add genes as objects on the plot
RCircos.Gene.Connector.Plot(genomic.data= g.data,track.num=1, side='in')
# add gene names
RCircos.Gene.Name.Plot(g.data, name.col=4,track.num=2, side='in')
# add -log10 P values
#RCircos.Line.Plot(g.data, data.col=5, track.num=3, side='in')
Link.Data <- Link.Data[with(data=Link.Data, order(as.numeric(gsub("^chr", "", chromosome_name)), start_position)), ]
mycolr <- rep(NA, nrow(Link.Data))
# missense allele A here is effect allele
mycolr[which(Link.Data$dir=="--")] <- "blue"
mycolr[which(Link.Data$dir=="++")] <- "red"
# add links
Link.Data$PlotColor <- mycolr
RCircos.Link.Plot(Link.Data, track.num= 5, by.chromosome=F, is.sort=T, lineWidth= 0.35*round(abs(Link.Data$effectsize*10)))
dev.off()
