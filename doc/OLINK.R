# 5-9-2018 JHZ

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
# TWEAK O43508 <- Q4ACW9
Inflammation["UniProt.No."] <- with(Inflammation, {replace(UniProt.No.,UniProt.No.=="Q4ACW9","O43508")})
Inflammation["Comment"] <- NA
Inflammation["Comment"] <- with(Inflammation, {replace(Comment, UniProt.No.=="O43508", "Q4ACW9")})
inf.orig <- Inflammation
# grep inf1 olink.prot.list.txt | sed 's/inf1_//g;s/___/\t/g' > inf1.list
inf <- read.table("inf1.list",header=FALSE,col.names=c("prot","UniProt"),sep="\t",as.is=TRUE)
inf1 <- merge(inf,inf.orig,by.x="UniProt",by.y="UniProt.No.")
# See https://www.uniprot.org/uniprot/ for additional information
write.csv(inf1[c("UniProt","prot","Target","Comment")], file="inf1.csv", quote=FALSE, row.names=FALSE)
# from CVD I analysis plan
cvd1 <- read.delim("cvd1.txt", as.is=TRUE)
cvd1 <- cvd1[c("Olink_name", "gene", "Uniprot")]
inf2 <- merge(cvd1, inf.orig, by.x="Uniprot", by.y="UniProt.No.")

tabs <-c("Cardiometabolic","Cell Regulation","CVD II","CVD III","Development","Immune Response","Immuno-Oncology",
         "Inflammation","Metabolism","Neurology","Oncology II","Organ Damage")
olink_panel(xlsx,tabs,TRUE,92,FALSE)
ls()
