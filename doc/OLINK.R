# 23-8-2018 JHZ

options(width=160)

olink <- function(xlsx, tabs, order=FALSE, nlines=5)
{
  for (x in tabs) 
  {
    cat("\n\n", x, ":\n", rep("-", nchar(x)+1), "\n\n", sep="")
    t <- openxlsx::read.xlsx(xlsx, sheet=x, colNames=TRUE, skipEmptyRows=FALSE, cols=c(1:16), rows=3:95)
    s <- sub(" ", "_", x)
    if (!order) assign(s, t, envir=.GlobalEnv) else
    {
      o <- order(t[,2])
      assign(s, t[o,], envir=.GlobalEnv)
    }
    s <- get(noquote(s))
    t <- "Target"
    n <- names(s)
    print(head(s["Target"],nlines),right=FALSE)
    cat("\n")
    print(head(s[setdiff(n,t)],nlines))
  }
}
xlsx <- "Olink validation data all panels.xlsx"
tabs <-c("Cardiometabolic","Cell Regulation","CVD II","CVD III","Development","Immune Response","Immuno-Oncology",
         "Inflammation","Metabolism","Neurology","Oncology II","Organ Damage")
olink(xlsx,tabs)
tabs <- "Inflammation"
olink(xlsx,tabs,nlines=92)
