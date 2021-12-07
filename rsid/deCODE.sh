#!/usr/bin/bash

export proteomics_results=~/rds/results/public/proteomics/deCODE
export v4=SomaLogicv4.tsv

if [ ! -d ${INF}/deCODE ]; then mkdir -p ${INF}/deCODE; fi

function panel()
{
R --no-save -q <<END
  options(width=200)
  f <- file.path(Sys.getenv("proteomics_results"),"ferkingstad21.xlsx")
  SomaLogicv4 <- openxlsx::read.xlsx(f,sheet=1,startRow=3,colNames=TRUE,cols=1:12)
  out <- file.path(Sys.getenv("INF"),"deCODE",Sys.getenv("v4"))
  write.table(SomaLogicv4,file=out,quote=FALSE,row.names=FALSE,sep="\t")
END
}

function replication()
(
  gunzip -c ${proteomics_results}/2155/10000_28_CRYBB2_CRBB2.txt.gz | head -1 | awk -vOFS="\t" '{print "prot",$0}'
  join -13 -22 <(cut -f1,4,7 --output-delimiter=' ' ${INF}/deCODE/SomaLogicv4.tsv | sort -k3,3) \
               <(cut -f2,4,5,20 --output-delimiter=' ' ${INF}/work/INF1.METAL | awk '{print $2":"$3,$4,$1}' | sort -k2,2) | \
  parallel -C' ' -j20 '
    tabix ${INF}/deCODE/{2}*.gz {4} | grep -w {5} | awk -vOFS="\t" -vprot={2} "{print prot,\$0}"
  '
) > ${INF}/deCODE/replication.tsv

function deCODE()
(
echo -e "Protein\tSentinels\tUniProt\tSNPid\tcis/trans\tProxies\tr2\tp\tTarget Full Name\tSource\tPMID\tComment"
join -12 -21 <(awk 'NR>1 && $9<=5e-8 {
                    split($4,a,":")
                    if(a[3]<a[4]) snpid=a[1]":"a[2]"_"a[3]"_"a[4];
                    else snpid==a[1]":"a[2]"_"a[4]"_"a[3]
                    print snpid, $1, $9
                 }' ${INF}/deCODE/replication.tsv | sort -k2,2) \
             <(cut -f1,7 --output-delimiter=' ' ${INF}/deCODE/${v4} | sed 's/-/_/' | sort -k1,1) | \
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
