#!/usr/bin/bash

function setup()
{
  if [ ! -f ${INF}/work/INF1.merge.genes ]; then
     grep -v BDNF ${INF}/doc/olink.inf.panel.annot.tsv | \
     cut -f3,8,9,10 | \
     sed 's/"//g' | \
     sort -k1,1 | \
     join -12 ${INF}/work/inf1.tmp - > ${INF}/work/INF1.merge.genes
  fi

  for type in cis trans pan
  do
    if [ ! -d ${INF}/mr/${type} ]; then mkdir -p ${INF}/mr/${type}; fi; 
    export suffix=${type};
    export lp=-7.30103
    cut -f3 ${INF}/work/INF1.METAL | sed '1d' | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    parallel -j5 --env suffix -C' ' '
      echo --- {2} ---
      (
        case ${suffix} in
        cis)
          gunzip -c ${INF}/METAL/{2}-1.tbl.gz | cut -f1-6,10-12,18 | \
          awk -vchr={3} -vstart={4} -vend={5} -vM=1e6 -vsuffix=${suffix} "(\$1==chr && \$2>=start-M && \$2 <= end+M)"
          ;;
        trans)
          gunzip -c ${INF}/METAL/{2}-1.tbl.gz | cut -f1-6,10-12,18 | \
          awk -vchr={3} -vstart={4} -vend={5} -vM=1e6 -vlogp=${lp} -vsuffix=${suffix} "!(\$1==chr && \$2>=start-M && \$2 <= end+M) && \$9<=logp"
          ;;
        pan)
          gunzip -c ${INF}/METAL/{2}-1.tbl.gz | cut -f1-6,10-12,18 | \
          awk -vchr={3} -vstart={4} -vend={5} -vM=1e6 -vlogp=${lp} -vsuffix=${suffix} "\$9<=logp"
          ;;
        esac
      ) >  ${INF}/mr/${suffix}/{2}-${suffix}.mri
      (
        echo -e "prot\trsid\tChromosome\tPosition\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tlogP\tN"
        awk "{\$4=toupper(\$4);\$5=toupper(\$5);print}" ${INF}/mr/${suffix}/{2}-${suffix}.mri | \
        sort -k3,3 | \
        join -23 ${INF}/work/INTERVAL.rsid - | \
        awk -v prot={2} "{\$1=prot;print}" | \
        tr " " "\t"
      ) | gzip -f > ${INF}/mr/${suffix}/{2}-${suffix}.mrx
    '
  done
}

function mr()
{
  parallel --env INF -C' ' '
    export MRBASEID={1}; 
    export prot={2}; 
    export type={3}; 
    export prefix={1}-{2}-{3};
    export suffix=${type};
    echo ${prefix}
    R --no-save <${INF}/rsid/mr.R 2>&1 | \
    tee ${INF}/mr/${type}/${prefix}.log
    for f in result loo single;
    do
      export z=${INF}/mr/${suffix}/${prefix}-${f};
      if [ -f ${z}.txt ]; then
        awk -vFS="\t" "NR==1||(\$9<=0.05 && \$9!=\"NA\")" ${z}.txt > ${z}.sig;
        export l=$(wc -l ${z}.sig | cut -d" " -f1);
        if [ ${l} -le 1 ]; then rm ${z}.sig; fi;
      fi
    done
  ' ::: $(awk -vFS="\t" 'NR>1 {print $4}' ${INF}/rsid/efo.txt) \
    ::: $(sed '1d' ${INF}/work/INF1.merge | cut -f5 | sort -k1,1 | uniq) \
    ::: cis trans pan
}

mr

function collect()
{
  for type in cis trans pan
  do
    export nrows=$(sed '1d' ${INF}/rsid/efo.txt | wc -l | cut -d' ' -f1)
    for i in $(seq ${nrows})
    do
      export trait=$(sed '1d' ${INF}/rsid/efo.txt | awk -vFS="\t" -vnr=${i} 'NR==nr{print $2}')
      export id=$(sed '1d' ${INF}/rsid/efo.txt | awk -vFS="\t" -vnr=${i} 'NR==nr{print $4}')
      echo ${id} -- ${trait}
      (
        cat ${INF}/mr/${type}/*result.txt | head -1 | \
        cut -f1,2 --complement | awk -v OFS="\t"  '{print $0, "cistrans"}'
        grep -h -w ${id} ${INF}/mr/${type}/*result.txt | grep -e Inverse -e ratio | \
        cut -f1,2 --complement | awk -v OFS="\t" -v type=${type} '{print $0, type}'
      ) > ${INF}/mr/${id}-${type}.result
      (
        cat ${INF}/mr/${type}/*single.txt | head -1 | \
        cut -f3,4 --complement | awk -v OFS="\t"  '{print $0, "cistrans"}'
        grep -h -w ${id} ${INF}/mr/${type}/*single.txt | awk '$NF!="NA"' | \
        cut -f3,4 --complement | awk -v OFS="\t" -v type=${type} '{print $0, type}'
      ) > ${INF}/mr/${id}-${type}.single
    done
  done

  (
    cat ${INF}/mr/*result | head -1
    grep -h -v cistrans ${INF}/mr/*result
  ) > ${INF}/mr/efo-result.txt

  export all=$(ls ${INF}/mr/cis/*result.txt ${INF}/mr/trans/*result.txt ${INF}/mr/pan/*result.txt | wc -l)
  export p=$(bc -l <<< 0.05/${all})
  echo ${all} ${p}
  awk -vp=${p} -vFS="\t" -vOFS="\t" '
  {
    if (index(FILENAME,"cis")) tag="cis";
    else if (index(FILENAME,"trans")) tag="trans";
    else if (index(FILENAME,"pan")) tag="pan"
  # id.exposure id.outcome outcome exposure method nsnp b se pval
    if ($(NF-1)<p) {print $1,$2,$3,$4,tag,$6,$7,$8,$9}
  }' ${INF}/mr/*result | \sed 's/|| id:/\t/' | xsel -i
}

R --no-save -q <<END
   library(dplyr)
   library(ggplot2)
   options(width=120)
   INF <- Sys.getenv("INF")
   efo <- read.delim(file.path(INF,"rsid","efo.txt")) %>%
          mutate(x=1:n()) %>%
          select(MRBASEID,trait,x)
   d3 <- read.delim(file.path(INF,"mr","efo-result.txt")) %>%
         filter(exposure=="IL.12B") %>%
         mutate(MRBASEID=unlist(lapply(strsplit(outcome,"id:"),"[",2)),y=b,
                col=case_when(cistrans=="cis" ~ "red",
                              cistrans=="trans" ~ "blue",
                              cistrans=="pan" ~ "black")) %>%
         mutate(cistrans=recode(cistrans,cis="(a) cis",trans="(b) trans",pan="(c) pan")) %>%
         select(-outcome,-method) %>%
         left_join(efo) %>%
         arrange(trait)
   p <- ggplot(d3,aes(y = trait, x = y))+
   theme_bw()+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y = element_text(size=14),text = element_text(size=15))+
   facet_wrap(~cistrans,ncol=3,scales="free_x")+
   geom_segment(aes(x = b-1.96*se, xend = b+1.96*se, yend = trait, colour=cistrans), show.legend=FALSE)+
   geom_vline(lty=2, aes(xintercept=0), colour = "red")+
   geom_point()+
   xlab("Effect size")+
   ylab("")
   ggsave(p,filename=file.path(INF,"mr","mr-IL.12B.png"),device="png",dpi=300,height=15,width=12)
   library(ggforestplot)
   select <- function(cistrans)
   {
      d <- filter(d3,cistrans==cistrans)
      forestplot(d, name=trait, estimate=b, se=se, pvalue=pval) +
      xlab("Effect size")
   }
   select("cis")
   select("trans")
   select("pan")
END
# d3[d3$cistrans=="pan","trait"] <- d3[d3$cistrans=="pan","x"]
# d3[d3$cistrans=="trans","trait"] <- d3[d3$cistrans=="trans","x"]
# unlist(gregexpr("[|]{1}","abc||"))
#  
# uncomment if clumping outside TwoSampleMR:
# cut -f3 mr/{2}-${suffix}.mri > mr/{2}-${suffix}.mrs
# plink --bfile INTERVAL/cardio/INTERVAL --extract mr/{2}-${suffix}.mrs \
#       --geno 0.1 --mind 0.1 --maf 0.005 --indep-pairwise 1000kb 1 0.01 --out mr/{2}-${suffix}
# plink-1.9 --bfile $f --clump $rt.tab --chr ${1} --from-bp ${2} --to-bp ${3} --clump-field P --clump-kb 500 --clump-p1 5e-8 --clump-r2 0 \
#           --clump-snp-field snpid --out $f
#   grep -w -f mr/{2}-${suffix}.prune.in mr/{2}-${suffix}.mri | \
#   join -23 <(zgrep "chr{3}" ${SUMSTATS}/snp150.snpid_rsid.gz) - | \
#   awk "{\$3=\"chr\"\$1\":\"\$2;print}" | \
#   join -23 -12 snp_pos - | \

# ukb-b-19657 (FEV1), N=421,986
# https://epigraphdb.org/pqtl/IL12B
# https://www.targetvalidation.org/evidence/ENSG00000113302/EFO_0000540?view=sec:known_drug