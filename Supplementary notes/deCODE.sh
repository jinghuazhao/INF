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

function collect()
{
  (
    awk -vOFS="\t" '{print "prot",$0}' ${deCODE}/doc/deCODE.hdr
    join -13 -22 <(cut -f1,4,7 --output-delimiter=' ' ${INF}/deCODE/SomaLogicv4.tsv | sort -k3,3) \
                 <(cut -f2,4,20,22 --output-delimiter=' ' ${INF}/work/INF1.METAL38 | awk '{print "chr"$2":"$4,$3,$1}' | sort -k2,2) | \
    parallel -C' ' -j20 '
      tabix ${INF}/deCODE/{2}.gz {4} | grep -w {5} | awk -vOFS="\t" -vprot={2} "{print prot,\$0}"
    '
  ) > ${INF}/deCODE/replication38.tsv
}

function replication()
{
cat << 'EOL' > ${INF}/deCODE/deCODE.R
    options(width=200)
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
          mutate(chr=seqnames.x,pos=end.y,Name=group_name) %>%
          select(Name,chr,pos,group) %>%
          full_join(rep38) %>%
          mutate(snpid=if_else(effectAllele>otherAllele,paste0(chr,":",pos,"_",otherAllele,"_",effectAllele),
                                                        paste0(chr,":",pos,"_",effectAllele,"_",otherAllele)),
                 chr=gsub("chr","",chr),
                 mlog10p=-gap::log10p(Beta/SE)) %>%
          select(-c(Name,Chrom,Pos,seqname,start,end,minus_log10_pval))
    write.table(select(df,snpid,prot,Pval),file=file.path(INF,"deCODE","replication.tsv"),
                row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    df_ids <- read.table(file.path(INF,"deCODE","olink_inf.full"),
                         col.names=c("uniprot","prot","gene","chrpos","rsid","snpid")) %>%
              left_join(df) %>%
              select(-c(group,rsids,chrpos)) %>%
              select(gene,rsid,chr,pos,effectAllele,otherAllele,ImpMAF,Beta,SE,mlog10p,Pval) %>%
              left_join(select(pQTLdata::inf1,gene,target.short)) %>%
              arrange(gene,rsid,Pval) %>%
              group_by(gene,rsid) %>%
              slice_head(n=1)
    deCODE <- select(df_ids,gene,rsid,chr,pos,effectAllele,otherAllele,Beta,SE,Pval,ImpMAF,mlog10p)
    INF1_METAL <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
                  left_join(select(pQTLdata::inf1,prot,gene)) %>%
                  mutate(Allele1=toupper(Allele1), Allele2=toupper(Allele2))
    tbl <- select(INF1_METAL,prot,rsid,Chromosome,Position,cis.trans) %>%
           left_join(select(pQTLdata::inf1,gene,prot,target.short)) %>%
           left_join(select(deCODE,gene,rsid,effectAllele,otherAllele,ImpMAF,Beta,SE,Pval)) %>%
           mutate(Protein=target.short,Pval=gap::pvalue(Beta/SE),Pval=if_else(Pval=="NAeNA","NA",Pval)) %>%
           select(Protein,Chromosome,Position,rsid,effectAllele,otherAllele,ImpMAF,Beta,SE,Pval)
    write.table(tbl,file=file.path(INF,"deCODE","deCODE.csv"),row.names=FALSE,quote=FALSE,sep=",")
    INF1_aristotle <- read.delim(file.path(INF,"aristotle","INF1.merge.replication.txt-rsid")) %>%
                      mutate(b_ARISTOTLE=BETA,se_ARISTOTLE=SE,P_ARISTOTLE=PVAL) %>%
                      mutate(prot=unlist(lapply(strsplit(Protein,"-"),"[",1)),
                             rsid=unlist(lapply(strsplit(Protein,"-"),"[",2)),) %>%
                      left_join(select(pQTLdata::inf1,prot,gene,target.short)) %>%
                      left_join(select(INF1_METAL,prot,rsid,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,cis.trans)) %>%
                      select(target.short,gene,rsid,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,cis.trans,
                             EFFECT_ALLELE,REFERENCE_ALLELE,b_ARISTOTLE,se_ARISTOTLE,P_ARISTOTLE)
    INF1_deCODE <- left_join(INF1_METAL,deCODE) %>%
                   mutate(b_deCODE=Beta,se_deCODE=SE,P_deCODE=Pval,mlog10p=-gap::log10p(b_deCODE/se_deCODE)) %>%
                   select(gene,rsid,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,effectAllele,otherAllele,-Beta,-SE,-Pval,
                          mlog10p,cis.trans,b_deCODE,se_deCODE,P_deCODE) %>%
                   mutate(sw=if_else(Allele1==effectAllele,1,-1),b_deCODE=sw*b_deCODE,sw2=(sign(Effect)==sign(b_deCODE))+0,
                   col=case_when(
                            mlog10p >= -log10(5e-8) ~ "red",
                            mlog10p >= -log10(0.05/180) & mlog10p <= -log10(5e-8) ~ "orange",
                            mlog10p  < -log10(0.05/180) ~ "grey",
                            TRUE ~ "white" )) %>%
                   filter(!is.na(b_deCODE))
    subset(INF1_deCODE[c("Effect","b_deCODE","log.P.","P_deCODE","cis.trans")],cis.trans=="cis") %>% arrange(Effect)
    nrow(INF1_deCODE)
    filter(INF1_deCODE[c("gene","Effect","b_deCODE","Allele1","Allele2","effectAllele","otherAllele","sw","log.P.","mlog10p")],
           mlog10p>=-log10(5e-8)) %>%
    filter(sign(Effect)!=sign(b_deCODE))
    filter(INF1_deCODE[c("gene","Effect","b_deCODE","Allele1","Allele2","effectAllele","otherAllele","sw","log.P.","mlog10p")],
           mlog10p>=-log10(5e-2/180)) %>%
    filter(sign(Effect)!=sign(b_deCODE))
    INF1_aristotle_deCODE <- left_join(select(INF1_aristotle,target.short,gene,rsid,Allele1,Allele2,Freq1,Effect,StdErr,cis.trans,
                                              EFFECT_ALLELE,REFERENCE_ALLELE,b_ARISTOTLE,se_ARISTOTLE,P_ARISTOTLE),
                                       select(INF1_deCODE,gene,rsid,effectAllele,otherAllele,b_deCODE,se_deCODE,P_deCODE,mlog10p))
    write.table(INF1_aristotle_deCODE,file=file.path(INF,"deCODE","INF1_aristotle_deCODE.csv"),row.names=FALSE,quote=FALSE,sep=",")
    filter(INF1_aristotle, is.na(P_ARISTOTLE)) %>% nrow()
    filter(INF1_aristotle, P_ARISTOTLE <= 5e-10) %>% nrow()
    filter(INF1_aristotle, P_ARISTOTLE <= 5e-8) %>% nrow()
    filter(INF1_aristotle, P_ARISTOTLE <= 0.05/180) %>% nrow()
    filter(INF1_aristotle, P_ARISTOTLE <= 5e-2) %>% nrow()
    filter(INF1_deCODE, is.na(mlog10p)) %>% nrow()
    filter(INF1_deCODE, mlog10p>=-log10(5e-10)) %>% nrow()
    filter(INF1_deCODE, mlog10p>=-log10(5e-8)) %>% nrow()
    filter(INF1_deCODE, mlog10p>=-log10(0.05/180)) %>% nrow()
    filter(INF1_deCODE, mlog10p>=-log10(5e-2)) %>% nrow()
    filter(INF1_aristotle_deCODE,is.na(P_ARISTOTLE) & is.na(mlog10p)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE <=5e-10 | mlog10p >= -log10(5e-10)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE <=5e-8 | mlog10p >= -log10(5e-8)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE <=0.05/180 | mlog10p >= -log10(0.05/180)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE <=5e-2 | mlog10p >= -log10(5e-2)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE > 5e-8 & mlog10p >= -log10(5e-8)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE <=5e-8 & mlog10p <  -log10(5e-8)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE > 5e-8 & mlog10p <  -log10(5e-8)) %>% nrow()
    filter(select(INF1_aristotle_deCODE,target.short,gene,rsid,Allele1,Allele2,effectAllele,otherAllele,
                  b_deCODE,se_deCODE,P_deCODE,P_ARISTOTLE,cis.trans),
           P_ARISTOTLE <=5e-8 & P_deCODE > 5e-8)
    SF_INF_deCODE <- function(d,f)
    {
      png(file=file.path(INF,"deCODE",f),res=300,width=15,height=15,units="in")
      attach(d)
      par(mar=c(5,5,1,1))
      plot(Effect,b_deCODE,pch=19,cex=2,col=col,xaxt="n",ann=FALSE,cex.axis=2)
      axis(1,cex.axis=2,lwd.tick=0.5)
      legend(x=1,y=-1,c("P<=5e-8","P >5e-8 & P <= 0.05/180","P>0.05/180"),
             box.lwd=0,cex=2,col=c("red","orange","grey"),pch=19)
      mtext("deCODE",side=2,line=3,cex=2)
      mtext("Meta-analysis",side=1,line=3,cex=2)
      abline(h=0,v=0)
      detach(d)
      dev.off()
    }
    all <- INF1_deCODE %>%
           filter(!is.na(b_deCODE)) %>%
           select(Effect,b_deCODE)
    cor(all,use="everything")
    cis <- INF1_deCODE %>%
           filter(!is.na(b_deCODE) & cis.trans=="cis") %>%
           select(Effect,b_deCODE)
    cor(cis,use="everything")
    trans <- INF1_deCODE %>%
             filter(!is.na(b_deCODE) & cis.trans=="trans") %>%
             select(Effect,b_deCODE)
    cor(trans,use="everything")
    SF_INF_deCODE(INF1_deCODE,"SF-INF-deCODE.png")
    SF_INF_deCODE(filter(INF1_deCODE,cis.trans=="cis"),"SF-INF-deCODE-cis.png")
    SF_INF_deCODE(filter(INF1_deCODE,cis.trans=="trans"),"SF-INF-deCODE-trans.png")
EOL
R --no-save < ${INF}/deCODE/deCODE.R
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
    with(region,table(sentinel,Prot))
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
    r1 <- subset(region,sentinel=="rs7406661" & rsid=="rs56115403")
    r2 <- subset(region,sentinel=="rs579459" & rsid=="rs977371848")
    r <- rbind(r1,r2)
    library(pQTLtools)
    si <- merge(SomaScanV4.1,inf1,by.x="UniProt.ID",by.y="uniprot") %>%
          mutate(SeqID=gsub("-","_",SeqID))
    annex <- left_join(r,select(si,SeqID,prot),by=c('Prot'='prot')) %>%
             select(MarkerName,SeqID,Pval)
    write.table(annex,file=file.path(INF,"deCODE","annex.tsv"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  '
}

# >     with(region,table(sentinel,Prot))
#              Prot
# sentinel      CCL19 IL.10 IL.12B MCP.3 SLAMF1 TWEAK  uPA
#   rs1950897       0     0     38     0      0     0    0
#   rs7406661       0     0      0     0      0     0   60
#   rs200489612     0     0      0     0      5     0    0
#   rs7213460       0     0      0  1320      0     0    0
#   rs570025519     0     0      0     0      1     0    0
#   rs28377109      0     1      0     0      0     0    0
#   rs11574915     19     0      0     0      0     0    0
#   rs579459        0     0      0     0      0    21    0
#
# hightlight
# [1] "RS_number"   "rs7406661"   "rs56093546"  "rs56115403"  "rs56244095"
# [6] "rs12945886"  "rs62061425"  "rs72837687"  "rs62061426"  "rs67143157"
#[11] "rs111072793" "rs55714927"
# [1] "rs7406661" "1"         "0.976"     "0.988"     "0.988"     "0.866"
# [7] "0.923"     "0.917"     "0.911"     "0.818"     "0.859"     "0.836"
# subset(region,sentinel=="rs7406661" & rsid%in%c("rs56115403","rs56244095"))
#      chr   pos37     Pos      rsids effectAllele otherAllele   Beta      Pval       SE Prot        MarkerName  sentinel       rsid
# 67 chr17 7063898 7160579 rs56115403            A           G 0.0881 2.552e-18 0.010092  uPA chr17:7063667_C_T rs7406661 rs56115403
# 68 chr17 7063899 7160580 rs56244095            A           C 0.0881 2.552e-18 0.010092  uPA chr17:7063667_C_T rs7406661 rs56244095
# rs579459
# [1] "RS_number"   "rs977371848" "rs579459"
# [1] "rs579459" "0.823"    "1"
#  subset(region,sentinel=="rs579459" & rsid=="rs977371848")
#       chr     pos37       Pos       rsids effectAllele otherAllele   Beta      Pval       SE  Prot         MarkerName sentinel        rsid
# 1464 chr9 136141870 133266456 rs977371848            C           T 0.0661 4.796e-08 0.012109 TWEAK chr9:136154168_C_T rs579459 rs977371848

function deCODE()
(
  echo -e "Protein\tSentinels\tUniProt\tSNPid\tcis/trans\tProxies\tr2\tp\tTarget Full Name\tSource\tPMID\tComment"
  join -12 -21 <(awk '$3<=5e-8' ${INF}/deCODE/replication.tsv ${INF}/deCODE/annex.tsv | sort -k2,2) \
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

function run_deCODE()
{
  deCODE | \
  awk -vFS="\t" -v OFS="\t" '
  {
    if($1=="TWEAK" && $2=="rs579459") {$6="rs977371848";$7=0.823}
    if($1=="uPA" && $2=="rs7406661") {$6="rs56115403";$7=0.988}
    print
  }' > ${INF}/deCODE/deCODE.tsv
}

replication
