# 14-7-2021 JHZ

if [ ! -d ${INF}/aristotle ]; then mkdir ${INF}/aristotle; fi
(
  gunzip -c ${INF}/sumstats/ARISTOTLE/*txt.gz | \
  head -1 | \
  awk -vOFS="\t" '{print "Protein",$0}'
  cut -f5,6 ${INF}/work/INF1.merge | \
  sed '1d' | \
  tr '\t' ' ' | \
  parallel -j8 -C' ' '
    echo {1}-{2}
    zgrep -w {2} ${INF}/sumstats/ARISTOTLE/ARISTOTLE.{1}.txt.gz;
    awk -vOFS="\t" -vprot={1} "{print prot,\$0}"
  '
) > ${INF}/aristotle/INF1.merge.replication.txt

# snpid --> rsid
cd ${INF}/aristotle
for f in INF1.merge.replication.txt
do
  awk -vOFS="\t" '{if(NF==1) printf $1 OFS; else print}' ${f} > ${f}-rsid
  (
  cat ${INF}/work/INF1.merge.rsid | \
  parallel --dry-run -C' ' "
    export s={1};
    export r={2};
    sed -i 's/'\"\${s}\"'/'\"\${r}\"'/g' ${f}-rsid
  "
  ) | bash
done
awk -vOFS="\t" '{if(NR>1) {split($1,a,"-");$1=a[1]};print}' INF1.merge.replication.txt-rsid | xsel -i
for p in 5e-10 1e-5 5e-2
do
  echo ${p}
  cut -f12 ${INF}/aristotle/INF1.merge.replication.txt-rsid | awk -v p=${p} '$1<p{print $1}' | wc -l
done
cd -

R --no-save -q <<END
  library(dplyr)
  INF <- Sys.getenv("INF")
  INF1_METAL <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
                mutate(Allele1=toupper(Allele1), Allele2=toupper(Allele2))
  aristotle <- read.delim(file.path(INF,"aristotle","INF1.merge.replication.txt-rsid")) %>%
               rename(n=N)
  INF1_aristotle <- mutate(INF1_METAL,Protein=paste0(prot,"-",rsid)) %>%
                    left_join(aristotle) %>%
                    select(Protein,Allele1,Allele2,EFFECT_ALLELE,REFERENCE_ALLELE,Freq1,Effect,StdErr,CODE_ALL_FQ,BETA,SE,cis.trans) %>%
                    mutate(sw=if_else(Allele1==REFERENCE_ALLELE,-1,1)) %>%
                    mutate(BETA=sw*BETA)
  all <- INF1_aristotle %>%
         filter(!is.na(BETA)) %>%
         select(Effect,BETA)
  cor(all,use="everything")
  cis <- INF1_aristotle %>%
         filter(!is.na(BETA) & cis.trans=="cis") %>%
         select(Effect,BETA)
  cor(cis,use="everything")
  trans <- INF1_aristotle %>%
           filter(!is.na(BETA) & cis.trans=="trans") %>%
           select(Effect,BETA)
  cor(trans,use="everything")
END
