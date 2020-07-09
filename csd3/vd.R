R --no-save -q <<END
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
  xlsx <- paste("doc/Olink validation data all panels.xlsx",sep="/")
  tabs <- c("Inflammation")
  olink_panel(xlsx,tabs,TRUE,92,FALSE)
  library(reshape)
  Inflammation <- rename(Inflammation,c(UniProt.No.="UniProt.No"))
  Olink <- Inflammation[c("Target","UniProt.No")]
  Olink <- within(subset(Olink,UniProt.No!="NA"), {OlinkID=UniProt.No; OlinkTarget=Target})
  HOME <- Sys.getenv("HOME")
  somalogic <- read.delim(paste(HOME,"SomaLogic","doc","SOMALOGIC_Master_Table_160410_1129info.tsv",sep="/"),as.is=TRUE)
  SomaLogic <- within(subset(somalogic,UniProt!="NA"),{SomaLogicID=UniProt;SomaLogicTarget=Target})
  olink_somalogic <- merge(Olink[c("UniProt.No","OlinkID","OlinkTarget")],
                           SomaLogic[c("SomaLogicTarget","UniProt","SomaLogicID")],
                           by.x="UniProt.No",by.y="UniProt")
  i <- setdiff(unique(with(olink_somalogic,UniProt.No)),"P23560")
  write(i,file="i")
  library(VennDiagram)
  plist <- list(setdiff(Olink[["UniProt.No"]],"P23560"),setdiff(SomaLogic[["UniProt"]],"P23560"))
  cnames <- c("Olink", "SomaLogic")
  venn.diagram(x = plist, category.names=cnames, filename='venn_diagram.png', imagetype="png", output=TRUE)
END
