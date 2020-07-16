# 16-7-2020 JHZ

options(width=160)

olink_panel <- function(xlsx, tabs, order, nlines, verbose)
{
  for (x in tabs) 
  {
    if (verbose) cat("\n\n", x, ":\n", rep("-", nchar(x)+1), "\n\n", sep="")
    t <- openxlsx::read.xlsx(xlsx, sheet=x, colNames=TRUE, skipEmptyRows=FALSE, cols=1:16, rows=3:95)
    s <- sub(" ", "_", x)
    if (!order) assign(s, t, envir=.GlobalEnv) else
    {
      o <- order(t[,2])
      assign(s, t[o,], envir=.GlobalEnv)
    }
    s <- get(s)
    t <- "Target"
    n <- names(s)
    if (verbose)
    {
      print(head(s["Target"],nlines),right=FALSE)
      cat("\n")
      print(head(s[setdiff(n,t)],nlines))
    }
  }
}
xlsx <- "Olink validation data all panels.xlsx"
tabs <- "Inflammation"
olink_panel(xlsx,tabs,FALSE,92,TRUE)
# TWEAK O43508 <- Q4ACW9. See https://www.uniprot.org/uniprot/ for additional information
Inflammation["UniProt.No."] <- with(Inflammation, {replace(UniProt.No.,UniProt.No.=="Q4ACW9","O43508")})
Inflammation["alias"] <- NA
Inflammation["alias"] <- with(Inflammation, {replace(alias, UniProt.No.=="O43508", "Q4ACW9")})
Inflammation["alias"] <- with(Inflammation, {replace(alias, UniProt.No.=="Q8NF90", "P12034")})
Inflammation["alias"] <- with(Inflammation, {replace(alias, UniProt.No.=="Q8WWJ7", "P30203")})
inf.orig <- Inflammation
inf <- read.table("inf1.list",header=FALSE,col.names=c("prot","UniProt"),sep="\t",as.is=TRUE)
inf1 <- merge(inf,inf.orig,by.x="UniProt",by.y="UniProt.No.")
write.csv(inf1[c("UniProt","prot","Target","alias")], file="inf1.csv", quote=FALSE, row.names=FALSE)
# UCSC hgTables
hgTables <- read.delim("hgTables.tsv",as.is=TRUE)
hgTables <- within(hgTables, UniProt <- unlist(lapply(strsplit(hgTables$name,"-"),"[",1)))
inf <- within(inf,UniProt <- replace(UniProt,UniProt=="Q8NF90","P12034"))
inf <- within(inf,UniProt <- replace(UniProt,UniProt=="Q8WWJ7","P30203"))
inf <- merge(inf,hgTables,by="UniProt",all=TRUE)
inf2 <- subset(inf,UniProt%in%inf1$UniProt|UniProt%in%c("P12034","P30203"))
write.csv(subset(inf2,!grepl("hap",X.chrom)), file="inf2.csv", quote=FALSE, row.names=FALSE)
# The tables are ordered below
tabs <-c("Cardiometabolic","Cell Regulation","CVD II","CVD III","Development","Immune Response","Immuno-Oncology",
         "Inflammation","Metabolism","Neurology","Oncology II","Organ Damage")
olink_panel(xlsx,tabs,TRUE,92,FALSE)
ls()

# Venn diagram with the SomaLogic panel
xlsx <- paste("doc/Olink validation data all panels.xlsx",sep="/")
tabs <- c("Inflammation")
olink_panel(xlsx,tabs,TRUE,92,FALSE)
library(reshape)
Inflammation <- rename(Inflammation,c(UniProt.No.="UniProt.No"))
Olink <- Inflammation[c("Target","UniProt.No")]
subset(Olink,UniProt.No=="O43508"|UniProt.No=="Q4ACW9")
toreplace <- with(Olink,UniProt.No=="Q4ACW9")
Olink[toreplace,"UniProt.No"] <- "O43508"
Olink <- within(subset(Olink,UniProt.No!="NA"), {OlinkID=UniProt.No; OlinkTarget=Target})
HOME <- Sys.getenv("HOME")
somalogic <- read.delim(paste(HOME,"SomaLogic","doc","SOMALOGIC_Master_Table_160410_1129info.tsv",sep="/"),as.is=TRUE)
SomaLogic <- within(subset(somalogic,UniProt!="NA"),{SomaLogicID=UniProt;SomaLogicTarget=Target})
olink_somalogic <- merge(Olink[c("UniProt.No","OlinkID","OlinkTarget")],
                         SomaLogic[c("SomaLogicTarget","UniProt","SomaLogicID")],
                         by.x="UniProt.No",by.y="UniProt")
# TWEAK-O43508 is missing
i <- setdiff(unique(with(olink_somalogic,UniProt.No)),"P23560")
write(i,file="i")
library(VennDiagram)
plist <- list(setdiff(Olink[["UniProt.No"]],"P23560"),setdiff(SomaLogic[["UniProt"]],"P23560"))
cnames <- c("Olink", "SomaLogic")
venn.diagram(x = plist, category.names=cnames, filename='venn_diagram.png', imagetype="png", output=TRUE)

## additional validation

Olink <- "INF/doc/olink.inf.panel.annot.tsv"
o <- read.delim(Olink, as.is=TRUE)
SomaLogic <- "SomaLogic/doc/SOMALOGIC_Master_Table_160410_1129info.tsv"
s <- read.delim(SomaLogic, as.is=TRUE)

library(reshape)
s <- rename(s, c(UniProt="uniprot"))
setdiff(intersect(o[["uniprot"]],s[["uniprot"]]),"P23560")

os <- merge(o,s,by="uniprot")
u <- setdiff(unique(os[["uniprot"]]),"P23560")
write(u,file="u")
length(u)
p <- unique(subset(os[c("uniprot","target.short")],uniprot%in%u))
dim(p)

library(dplyr)
nj <- nest_join(o,s,by="uniprot")
nj[order(nj[["uniprot"]]),"uniprot"]
unlist(lapply(nj$s,"[[",7))

library(VennDiagram)
plist <- list(setdiff(o[["uniprot"]],"P23560"),setdiff(s[["uniprot"]],c(NA,"P23560")))
cnames=c("Olink", "SomaLogic")
venn.diagram(x = plist, category.names=cnames, filename='Olink-SomaLogic-Venn-diagram.png', imagetype="png", output=TRUE)

xlsx <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0175-2/MediaObjects/41586_2018_175_MOESM4_ESM.xlsx"
t5 <- openxlsx::read.xlsx(xlsx, sheet=5, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:19), rows=c(3:2746))
t6 <- openxlsx::read.xlsx(xlsx, sheet=6, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:20), rows=c(3:167))
