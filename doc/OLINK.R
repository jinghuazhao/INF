# 23-8-2018 JHZ

options(width=160)

olink <- function(xlsx, tabs, order=TRUE, nlines=5)
{
  for (s in 1:length(tabs)) 
  {
    cat("\n\n", tabs[s], ":\n", rep("-", nchar(tabs[s])+1), "\n\n", sep="")
    t <- openxlsx::read.xlsx(xlsx, sheet=s, colNames=TRUE, skipEmptyRows=FALSE, cols=c(1:16), rows=3:95)
    if (!order) assign(tabs[s], t, envir=.GlobalEnv) else
    {
      o <- order(t[,2])
      assign(tabs[s], t[o,], envir=.GlobalEnv)
    }
    s <- get(noquote(tabs[s]))
    t <- "Target"
    n <- names(s)
    print(head(s["Target"]),right=FALSE)
    cat("\n")
    print(head(s[setdiff(n,t)],nlines))
  }
}
xlsx <- "Olink validation data all panels.xlsx"
tabs <-c("Cardiometabolic","Cell_Regulation","CVDII","CVDIII","Development","Immune_Respone","Immue_Oncology",
         "Inflammation","Metabolism","Neurology","OncologyII","Organ_Damage")
olink(xlsx,tabs,FALSE)
olink(xlsx,tabs,nlines=92)
