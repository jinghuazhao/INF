# 9-1-2020 JHZ

zcat Whole_Blood.allpairs.txt.gz | cut -f1 | awk 'NR>1{print substr($1,1,15)}' | sort | uniq > work/wb
R --no-save -q <<END
   wb <- scan("work/wb", what="")
   g <- within(subset(grex::grex(wb),!is.na(uniprot_id)),{uniprot=trimws(uniprot_id)})
   m <- merge(gap::inf1,g,by="uniprot",all.x=TRUE)
END
