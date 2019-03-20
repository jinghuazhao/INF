# 3-1-2019 JHZ

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
