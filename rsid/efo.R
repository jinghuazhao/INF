options(width=200)

INF <- Sys.getenv("INF")
url <- "https://jhz22.user.srcf.net/INF1.latest.xlsx"
efo <- subset(openxlsx::read.xlsx(url, sheet="EFO", colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:5), rows=c(2:79)),!is.na(MRBASEID))

library(gwasrapidd)
library(dplyr)
efo <- within(efo,{id=gsub(":","_",id)})
ebicat_trait_list <- with(efo, get_traits(efo_id=id))
efo_uri <- merge(efo,ebicat_trait_list@traits[c("efo_id","uri")],by.x="id",by.y="efo_id",all.x=TRUE)

write.table(efo_uri, file=file.path(INF,"rsid","efo.txt"),quote=FALSE,row.names=FALSE,sep="\t")
