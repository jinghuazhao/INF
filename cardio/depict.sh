# 21-5-2019 JHZ

export INF=$HOME/INF
if [ ! -d $INF/snps/cojo/depict ]; then
   mkdir $INF/snps/cojo/depict
fi
cd $INF/snps/cojo/depict
ln -sf $INF/snps/cojo/INF1.jma.cis
ln -sf $INF/snps/cojo/INF1.jma.tbl
sort INF1.jma.cis | \
uniq | \
join - $INF/work/INTERVAL.rsid > depict.rsid
R --no-save <<END
    require(dplyr)
    t <- read.delim("INF1.jma.tbl",as.is=TRUE)
    tbl <- within(t, {
      prot <- sapply(strsplit(Chromosome,":"),"[",1)
      Chromosome <- sapply(strsplit(Chromosome,":"),"[",2)
    })
    rsid <- read.table("depict.rsid",as.is=TRUE,col.names=c("MarkerName","rsid"))
    m <- within(nest_join(tbl,rsid),{rsid <- unlist(lapply(lapply(y,"[[",1),"[",1))})
    isna <- with(m, is.na(rsid))
    t <- within(m, {rsid[isna] <- MarkerName[isna]})
    h <- c("Allele1","Allele2","Freq1","Effect","StdErr","P.value","N","Chromosome","Position")
    s <- subset(aggregate(x = t[h], by = list(with(t,rsid)), FUN = "min"), substr(Group.1,1,2)=="rs")
    names(s)[1] <- "SNP"
    write.table(s,file="INF1.sumstats",col.names=FALSE,row.names=FALSE,quote=FALSE)
END
cd -
