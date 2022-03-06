#!/usr/bin/bash

export deCODE=~/rds/results/public/proteomics/deCODE
export v4=SomaLogicv4.tsv

if [ ! -d ${INF}/deCODE ]; then mkdir -p ${INF}/deCODE; fi

function panel()
{
R --no-save -q <<END
  options(width=200)
  f <- file.path(Sys.getenv("deCODE"),"doc","ferkingstad21.xlsx")
  SomaLogicv4 <- openxlsx::read.xlsx(f,sheet=1,startRow=3,colNames=TRUE,cols=1:12)
  out <- file.path(Sys.getenv("INF"),"deCODE",Sys.getenv("v4"))
  write.table(SomaLogicv4,file=out,quote=FALSE,row.names=FALSE,sep="\t")
END
cat olink_inf.lst | parallel -C' ' 'ln -s ${deCODE}/{}*gz ${INF}/deCODE/{}'
cat olink_inf.lst | parallel -C' ' 'ln -s ${deCODE}/{}*gz.tbi ${INF}/deCODE/{}.gz.tbi'
}

function replication()
{
  (
    awk -vOFS="\t" '{print "prot",$0}' ${deCODE}/doc/deCODE.hdr
    join -13 -22 <(cut -f1,4,7 --output-delimiter=' ' ${INF}/deCODE/SomaLogicv4.tsv | sort -k3,3) \
                 <(cut -f2,4,20,22 --output-delimiter=' ' ${INF}/work/INF1.METAL38 | awk '{print "chr"$2":"$4,$3,$1}' | sort -k2,2) | \
    parallel -C' ' -j20 '
      tabix ${INF}/deCODE/{2}.gz {4} | grep -w {5} | awk -vOFS="\t" -vprot={2} "{print prot,\$0}"
    '
  ) > ${INF}/deCODE/replication38.tsv
  Rscript -e '
    INF <- Sys.getenv("INF")
    suppressMessages(library(dplyr))
    rep38 <- read.delim(file.path(INF,"deCODE","replication38.tsv")) %>%
             mutate(seqname=Chrom,start=as.integer(Pos-1),end=Pos) %>%
                    arrange(Chrom,Pos)
    gr <- with(rep38,GenomicRanges::GRanges(seqnames=seqname,IRanges::IRanges(start,end,names=Name))) %>% unique()
    suppressMessages(library(rtracklayer))
    path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
    ch <- import.chain(path)
    seqlevelsStyle(gr) <- "UCSC"
    gr19 <- liftOver(gr,ch)
    df <- mutate(data.frame(gr),group_name=names(gr)) %>%
          left_join(data.frame(gr19),by=c("group_name")) %>%
          mutate(chr=seqnames.x,pos37=end.y,Name=group_name) %>%
          select(Name,chr,pos37,group) %>%
          full_join(rep38) %>%
          mutate(snpid=if_else(effectAllele>otherAllele,paste0(chr,":",pos37,"_",otherAllele,"_",effectAllele),
                                                        paste0(chr,":",pos37,"_",effectAllele,"_",otherAllele)),
                 pval=gap::pvalue(Beta/SE)) %>%
          select(snpid,prot,pval)
    write.table(df,file=file.path(INF,"deCODE","replication.tsv"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  '
}

function region()
{
  (
    awk -v OFS='\t' '{print $0,"Prot","MarkerName","sentinel"}' ${deCODE}/doc/deCODE.hdr
    join -13 -21 <(cut -f1,4,7 --output-delimiter=' ' ${INF}/deCODE/SomaLogicv4.tsv | sort -k3,3) \
                 <(Rscript -e '
                      INF <- Sys.getenv("INF")
                      suppressMessages(library(dplyr))
                      load(file.path(INF,"work","novel_data.rda"))
                      metal <- read.delim(file.path(INF,"work","INF1.METAL38")) %>%
                               select(prot,MarkerName,pos38)
                      novel_data <- left_join(novel_data,metal) %>%
                                    mutate(region=paste0("chr",Chromosome,":",pos38-1e6,"-",pos38+1e6)) %>%
                                    select(uniprot,region,prot,MarkerName,rsid) %>%
                                    arrange(uniprot)
                      write.table(novel_data,row.names=FALSE,col.names=FALSE,quote=FALSE)
                   ') | \
    parallel -C' ' -j20 '
      tabix ${INF}/deCODE/{2}.gz {4} | awk -v prot={5} -vsnpid={6} -vrsid={7} -vOFS="\t" "{print \$0,prot,snpid,rsid}"
    '
  ) > ${INF}/deCODE/region38.tsv
  Rscript -e '
    INF <- Sys.getenv("INF")
    suppressMessages(library(dplyr))
    rep38 <- read.delim(file.path(INF,"deCODE","region38.tsv")) %>%
             filter(Pval<=5e-8) %>%
             select(-minus_log10_pval,-N,-ImpMAF) %>%
             mutate(seqname=Chrom,start=as.integer(Pos-1),end=Pos) %>%
                    arrange(Chrom,Pos)
    gr <- with(rep38,GenomicRanges::GRanges(seqnames=seqname,IRanges::IRanges(start,end,names=Name))) %>% unique()
    suppressMessages(library(rtracklayer))
    path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
    ch <- import.chain(path)
    seqlevelsStyle(gr) <- "UCSC"
    gr19 <- liftOver(gr,ch)
    df <- transmute(data.frame(gr),group_name=names(gr)) %>%
          left_join(data.frame(gr19)) %>%
          filter(!is.na(seqnames)) %>%
          mutate(chr=seqnames,pos37=end,Name=group_name) %>%
          select(Name,chr,pos37,group) %>%
          left_join(select(rep38,-seqname,-start,-end)) %>%
          select(-Name,-group,-Chrom)
    write.table(subset(df,!is.na(rsids)),file=file.path(INF,"deCODE","region.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
  '
  Rscript -e '
    INF <- Sys.getenv("INF")
    suppressMessages(library(dplyr))
    region <- read.delim(file.path(INF,"deCODE","region.tsv")) %>%
              mutate(rsid=gsub("_[A-Z]*_[A-Z]*","",rsids))
    key <- Sys.getenv("LDLINK_TOKEN")
    sentinels <- unique(region$sentinel)
    print(sentinels)
    blocks <- r <- list()
    check <- function(snps)
    {
      sentinel_and_snps <- c(sentinels[i],snps[grepl("^rs",snps)])
      r[[i]] <- LDlinkR::LDmatrix(snps=sentinel_and_snps,pop="EUR",r2d="r2",token=key)
      r2 <- subset(r[[i]],RS_number==sentinels[i])
      sel <- !is.na(r2) & r2>=0.8
      cat(sentinels[i],"\n")
      print(names(r2)[sel])
      print(r2[sel])
    }
    for(i in 1:length(sentinels))
    {
      blocks[[i]] <- subset(region,sentinel==sentinels[i])
      if(nrow(blocks[[i]])<1000)
      {
        check(blocks[[i]]$rsid)
     } else {
        check(blocks[[i]]$rsid[1:800])
        check(blocks[[i]]$rsid[801:1320])
     }
    }
  '
}

# sentinels
# [1] "rs1950897"   "rs7406661"   "rs200489612" "rs7213460"   "rs570025519" "rs28377109"  "rs11574915"  "rs579459"
# hightlight
# [1] "RS_number"   "rs7406661"   "rs56093546"  "rs56115403"  "rs56244095"
# [6] "rs12945886"  "rs62061425"  "rs72837687"  "rs62061426"  "rs67143157"
#[11] "rs111072793" "rs55714927"
# [1] "rs7406661" "1"         "0.976"     "0.988"     "0.988"     "0.866"
# [7] "0.923"     "0.917"     "0.911"     "0.818"     "0.859"     "0.836"
# rs579459
# [1] "RS_number"   "rs977371848" "rs579459"
# [1] "rs579459" "0.823"    "1"

function deCODE()
(
  echo -e "Protein\tSentinels\tUniProt\tSNPid\tcis/trans\tProxies\tr2\tp\tTarget Full Name\tSource\tPMID\tComment"
  join -12 -21 <(awk '$3<5e-8' ${INF}/deCODE/replication.tsv | sort -k2,2) \
               <(cut -f1,7 --output-delimiter=' ' ${INF}/deCODE/${v4} | sort -k1,1) | \
  awk '{print $4"-"$2,$0}' | \
  sort -k1,1 | \
  join - <(awk 'NR>1 {print $20"-"$1,$3,$2,$21}' ${INF}/work/INF1.METAL | sort -k1,1) | \
  awk 'a[$1]++==0' | \
  cut -d' ' -f1,2 --complement | \
  sort -k4,4 | \
  tr ' ' '\t'| \
  join -24 <(Rscript -e 'write.table(pQTLtools::inf1[c("prot","target.short")],col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")' | \
             sort -t$'\t' -k1,1) - | \
  awk -vOFS='\t' '{print $2,$6,$5,$3,$7,"as sentinels",1,$4}' | \
  sort -t$'\t' -k1,1 | \
  join -t$'\t' - \
               <(Rscript -e 'write.table(data.frame(pQTLtools::inf1[c("target.short","target")],Source="Ferkingstad et al. (2021)",PMID="",Comment=""),
                                         col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")' | \
                 sort -t$'\t' -k1,1)
)

deCODE > ${INF}/deCODE/deCODE.tsv
