# 9-4-2019 JHZ

source tryggve/analysis.ini

function fp()
{
  (
    gunzip -c METAL/4E.BP1-1.tbl.gz | \
    head -1
    awk 'NR>1' clumping/INF1.clumped | \
    cut -d' ' -f1,3 | \
    parallel -j4 -C' ' 'zgrep -w -H {2} METAL/{1}-1.tbl.gz'
  ) | \
  sed 's|METAL/||g;s/-1.tbl.gz//g' > INF1.clumped.tbl
  (
    awk 'NR>1' INF1.clumped.tbl | \
    cut -f1,3,13 | \
    awk '{split($1,a,":");print a[1],$2,$3}'
    parallel -j4 -C' ' '
      export direction=$(zgrep -w {2} METAL/{1}-1.tbl.gz | cut -f13)
      let j=1
      for i in $(grep "Input File" METAL/{1}-1.tbl.info | cut -d" " -f7)
      do
         export n=$(awk -vj=$j "BEGIN{split(ENVIRON[\"direction\"],a,\"\");print a[j]}")
         if [ "$n" != "?" ]; then zgrep -H -w {2} $i; fi
         let j=$j+1
      done
  '
  ) | \
  sed 's|/data/jinhua/INF/sumstats||g;s/.gz//g' > INF1.clumped.all
  R -q --no-save <<\ \ END
    tbl <- read.delim("INF1.clumped.tbl",as.is=TRUE)
    tbl <- within(tbl, {
      prot <- sapply(strsplit(Chromosome,":"),"[",1)
      Chromosome <- sapply(strsplit(Chromosome,":"),"[",2)
    })
    all <- read.table("INF1.clumped.all",as.is=TRUE,
           col.names=c("SNPID", "CHR", "POS", "STRAND", "N", "EFFECT_ALLELE", "REFERENCE_ALLELE",
                       "CODE_ALL_FQ", "BETA", "SE", "PVAL", "RSQ", "RSQ_IMP", "IMP"))
    all <- within(all, {
      dir.study.prot <- sapply(strsplit(SNPID,":"),"[",1)
      p1 <- sapply(strsplit(SNPID,":"),"[",2)
      p2 <- sapply(strsplit(SNPID,":"),"[",3)
      MarkerName <- paste(p1,p2,sep=":")
      study <- sapply(strsplit(dir.study.prot,"/"),"[",2)
      study.prot <- sapply(strsplit(dir.study.prot,"/"),"[",3)
      substudy <- sapply(strsplit(study.prot,"[.]"),"[",1)
      pos <- unlist(lapply(gregexpr("[.]",study.prot),"[",1))
      prot <- substring(study.prot,pos+1)
    })
    require(rmeta)
    pdf("INF1.fp.pdf")
    xlim <- c(-1.5,1.5)
    for(i in 1:nrow(tbl))
    {
       p <- tbl[i,"prot"]
       m <- tbl[i,"MarkerName"]
       d <- gsub("[?]","",tbl[i,"Direction"])
       s <- unlist(strsplit(d,""))
       f <- as.numeric(paste0(s,1))
       A1 <- toupper(tbl[i,"Allele1"])
       A2 <- toupper(tbl[i,"Allele2"])
       print(paste0(i,"-",p,":",m))
       with(subset(all,prot==p & MarkerName==m), {
         e <- toupper(EFFECT_ALLELE)
         r <- toupper(REFERENCE_ALLELE)
         a1 <- a2 <- vector('character',length(e))
         a1 <- e
         a2 <- r
         c <- rep(1,length(e))
         j <- sapply(a1,'!=',A1)
         a1[j] <- r[j]
         a2[j] <- e[j]
         c[j] <- -1
         print(cbind(A1,A2,EFFECT_ALLELE,REFERENCE_ALLELE,BETA,a1,a2,BETA*c))
         BETA <- BETA * c
         tabletext <- cbind(c("Study",study,"Summary"),
                              c("Effect",format(BETA,digits=3),format(tbl[i,"Effect"],digits=3)),
                              c("SE",format(SE,digits=3),format(tbl[i,"StdErr"],digits=3)),
                              c("N",N,tbl[i,"N"]))
         print(tabletext)
         forestplot(tabletext,
                    c(NA,BETA,tbl[i,"Effect"]),
                    c(NA,BETA-1.96*SE,tbl[i,"Effect"]-1.96*tbl[i,"StdErr"]),
                    c(NA,BETA+1.96*SE,tbl[i,"Effect"]+1.96*tbl[i,"StdErr"]),
                    zero=0,
                    is.summary=c(TRUE,rep(FALSE,length(BETA)),TRUE),
                    boxsize=0.75,
                    col=meta.colors(box="royalblue",line="darkblue", summary="royalblue"))
         title(paste0(p," [",m," (",A1,"/",A2,")","]"))
         metaplot(BETA,SE,N,
                  labels=paste0(study," (",format(BETA,digits=3),"/",format(SE,digits=3),")"),
                  xlab="Effect distribution",ylab="",xlim=xlim,
                  summn=tbl[i,"Effect"],sumse=tbl[i,"StdErr"],sumnn=tbl[i,"N"],
                  colors=meta.colors(box="red",lines="blue", zero="green", summary="red", text="black"))
         title(paste0(p," [",m," (",A1,"/",A2,")", " N=",tbl[i,"N"],"]"))
       })
    }
    dev.off()
  END
}
