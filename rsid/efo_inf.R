url <- "https://jhz22.user.srcf.net/INF1.latest.xlsx"
efo <- subset(openxlsx::read.xlsx(url, sheet="EFO", colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:4), rows=c(1:78)),!is.na(MRBASEID))

