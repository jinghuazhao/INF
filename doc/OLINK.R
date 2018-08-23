require(openxlsx)
xlsx <- "Olink validation data all panels.xlsx"
# Inflammation
tabs <- c("Cardiometabolic","Cell_Regulation","CVDII","CVDIII","Development","Immune_Respone","Immue_Oncology",
          "Inflammation","Metabolism","Neurology","OncologyII","Organ_Damage")
tabs
for (s in 1:length(tabs)) 
{
  t <- read.xlsx(xlsx, sheet=s, colNames=TRUE, skipEmptyRows=FALSE, cols=c(1:16), rows=3:95)
  assign(tabs[s],t)
}
cols <- 1:5
head(CVDII[cols])
head(CVDIII[cols])
head(Inflammation[cols])
