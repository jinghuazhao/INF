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
    if [ ! -d ${INF}/mr/${type} ]; then mkdir -p ${INF}/mr/${type}; fi 
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
        if [ ${l} -le 1 ]; then rm ${z}.sig; fi
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

  for type in cis trans pan
  do
  (
    cat ${INF}/mr/*${type}.result | head -1
    grep -h -v cistrans ${INF}/mr/*${type}.result
  ) > ${INF}/mr/${type}-efo-result.txt
  done

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

# pheatmap
R --no-save -q <<END
   INF <- Sys.getenv("INF")
   library(dplyr)
   library(stringr)
   efo <- read.delim(file.path(INF,"rsid","efo.txt"))
   mr <- read.delim(file.path(INF,"mr","efo-result.txt")) %>%
         filter(cistrans=="cis") %>%
         mutate(fdr=p.adjust(pval,method="fdr"),
                delimiter=str_locate(outcome,"\\|\\|"),
                mrbaseid=substring(outcome,delimiter[,2]+2),
                mrbaseid=gsub("id:","",mrbaseid)
               ) %>%
         left_join(efo,by=c("mrbaseid"="MRBASEID")) %>%
         left_join(gap::inf1[c("prot","target.short")],by=c("exposure"="prot")) %>%
         mutate(outcome=paste0(id," (",trait,")"),
                exposure=target.short,
                group=as.numeric(cut(b,breaks=quantile(b,seq(0,1,0.125)))),
                log10p=sign(b)*(-log10(pval))
               ) %>%
         select(exposure,outcome,method,b,se,pval,nsnp,fdr,group,log10p)
   exposure <- unique(with(mr,exposure))
   outcome <- unique(with(mr,outcome))
   n <- length(exposure)
   m <- length(outcome)
   mr_mat <- matrix(NA,m,n)
   colnames(mr_mat) <- exposure
   rownames(mr_mat) <- outcome
   for(k in 1:nrow(mr))
   {
      t <- mr[k,c('exposure','outcome','b','group','log10p')]
      i <- t[['outcome']]
      j <- t[['exposure']]
      v <- t[['log10p']]
      mr_mat[i,j] <- v
   }
   options(width=200)
   subset(mr,fdr<=0.05)
   library(pheatmap)
   png(file.path(INF,"mr","efo-cis.png"),res=300,width=30,height=15,units="in")
   pheatmap(mr_mat,cluster_rows=FALSE,cluster_cols=FALSE,angle_col="315",fontsize_row=18,fontsize_col=18)
   dev.off()
END

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

# All GSMR

export suffix=cis
if [ ! -d ${INF}/mr/gsmr ]; then mkdir -p ${INF}/mr/gsmr; fi

function opengwas()
{
  cd ${INF}/OpenGWAS
  for id in $(sed '1d' ${INF}/rsid/efo.txt | cut -f4)
  do
    export f=https://gwas.mrcieu.ac.uk/files/${id}/${id}.vcf.gz
    if [ ! -f ${f} ]; then wget ${f}; fi
    if [ ! -f ${f}.tbi ]; then wget ${f}.tbi; fi
  done
  cd -
}

function mrx()
{
  if [ ! -d ${INF}/mr/gsmr/mrx ]; then mkdir -p ${INF}/mr/gsmr/mrx; fi
  awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
  parallel -j5 --env suffix -C' ' '
      echo --- {2} ---
      (
        gunzip -c ${INF}/METAL/{2}-1.tbl.gz | cut -f1-6,10-12,18 | \
        awk -vchr={3} -vstart={4} -vend={5} -vM=1e6 -vsuffix=${suffix} "(\$1==chr && \$2>=start-M && \$2 <= end+M)"
      ) >  ${INF}/mr/gsmr/{2}-${suffix}.mri
      (
        echo -e "prot\trsid\tChromosome\tPosition\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tlogP\tN"
        awk "{\$4=toupper(\$4);\$5=toupper(\$5);print}" ${INF}/mr/${suffix}/{2}-${suffix}.mri | \
        sort -k3,3 | \
        join -23 ${INF}/work/INTERVAL.rsid - | \
        awk -v prot={2} "{\$1=prot;print}" | \
        tr " " "\t"
      ) | gzip -f > ${INF}/mr/gsmr/mrx/{2}-${suffix}.mrx
    '
}

function ref_prot_outcome_gsmr()
{
  export suffix=cis
  awk '$21==ENVIRON["suffix"] {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
  awk -vM=1e6 '{print $2, $3, $4-M, $5+M}' | \
  while read prot chr start end
  do
    export prot=${prot}
    export chr=${chr}
    export start=${start}
    export end=${end}
    echo ${chr} ${start} ${end} > ${INF}/mr/gsmr/ref/${prot}.bed1
    plink2 --bfile ${INF}/work/INTERVAL --extract bed1 ${INF}/mr/gsmr/ref/${prot}.bed1 \
           --make-bed --rm-dup force-first list --out ${INF}/mr/gsmr/ref/${prot}
  done
  awk '$21==ENVIRON["suffix"] {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
  parallel -j15 --env INF --env suffix -C' ' '
    (
      echo -e "SNP A1 A2 freq b se p N"
      cut -f1,2 --complement ${INF}/mr/gsmr/mrx/{2}-${suffix}.mri | \
      awk "{\$2=toupper(\$2);\$3=toupper(\$3);\$7=10^\$7;print}"
    ) | gzip -f > ${INF}/mr/gsmr/prot/{2}.gz
  '
  if [ ! -d ${INF}/mr/gsmr/trait ]; then mkdir -p ${INF}/mr/gsmr/trait; fi
  awk -vFS="\t" 'NR>1 {print $4,$5+$6}' ${INF}/rsid/efo.txt | \
  while read efo N
  do
    export efo=${efo}
    export N=${N}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 "{print \$1, \$2, \$3\":\"\$4-M\"-\"\$5+M}" | \
    parallel -C' ' -j15 --env INF --env efo --env N --env suffix '
      (
        echo -e "SNP A1 A2 freq b se p N"
        bcftools query -f "%CHROM %POS %ID %ALT %REF [%AF] [%ES] [%SE] [%LP] [$N] \n" -r {3} ${INF}/OpenGWAS/${efo}.vcf.gz | \
        awk "{
          if(\$4<\$5) snpid=\"chr\"\$1\":\"\$2\"_\"\$4\"_\"\$5; else snpid=\"chr\"\$1\":\"\$2"_"\$4"_"\$5
          \$7=10^-\$7
          print snpid, \$4, \$5, \$6, \$7, \$8, \$9, \$10
        }"
      ) | \
      gzip -f> ${INF}/mr/gsmr/trait/${efo}-{2}.gz
    '
  done
# sbatch --wait ${INF}/rsid/mr.sb
  (
    cat *.gsmr | head -1
    grep -e Exposure -e nan -v *gsmr | tr ':' '\t' | cut -f1 --complement
  ) > gsmr.txt
  R --no-save -q <<\ \ END
    INF <- Sys.getenv("INF")
    efo <- read.delim(file.path(INF,"rsid","efo.txt"))
    gsmr <- read.delim("gsmr.txt")
    library(dplyr)
    gsmr_efo <- gsmr %>%
                left_join(pQTLtools::inf1[c("prot","target.short")], by=c("Exposure"="prot")) %>%
                left_join(efo,by=c("Outcome"="MRBASEID")) %>%
                rename(protein=Exposure,MRBASEID=Outcome) %>%
                mutate(protein=target.short,fdr=p.adjust(p,method="fdr")) %>%
                select(protein,MRBASEID,trait,bxy,se,p,nsnp,fdr,Ncases,Ncontrols,id,uri,Zhengetal) %>%
                arrange(fdr)
    subset(gsmr_efo[setdiff(names(gsmr_efo),c("Zhengetal","uri"))],fdr<=0.05)
    write.table(gsmr_efo,"gsmr-efo.txt",row.names=FALSE,quote=FALSE,sep="\t")
  END
}
# zgrep SAMPLE ${INF}/OpenGWAS/*.vcf.gz | cut -f1,7,8 > ${INF}/OpenGWAS/ieu.sample
# sed 's/.vcf.gz:/\t/;s/TotalControls=/\t/;s/,TotalCases=/\t/;s/,StudyType/\t/' ieu.sample | cut -f1,3,4 > ${INF}/OpenGWAS/ieu.N
#   zgrep -h SAMPLE ${INF}/OpenGWAS/${efo}.vcf.gz | cut -f1,7,8 > ${INF}/mr/gsmr/trait/${efo}.sample

R --no-save -q <<END
   INF <- Sys.getenv("INF")
   library(dplyr)
   library(stringr)
   gsmr <- read.delim(file.path(INF,"mr","gsmr","out","5e-8","gsmr-efo.txt")) %>%
           mutate(outcome=paste0(id," (",trait,")"),
                  exposure=protein,
                  group=as.numeric(cut(bxy,breaks=quantile(bxy,seq(0,1,0.125))))) %>%
           select(exposure,outcome,bxy,se,p,nsnp,fdr,group)
   exposure <- unique(with(gsmr,exposure))
   outcome <- unique(with(gsmr,outcome))
   n <- length(exposure)
   m <- length(outcome)
   gsmr_mat <- matrix(NA,m,n)
   colnames(gsmr_mat) <- exposure
   rownames(gsmr_mat) <- outcome
   for(k in 1:nrow(gsmr))
   {
      t <- gsmr[k,c('exposure','outcome','bxy','group','fdr')]
      i <- t[['outcome']]
      j <- t[['exposure']]
      v <- t[['bxy']]
      gsmr_mat[i,j] <- v
   }
   options(width=200)
   subset(gsmr,fdr<=0.05)
   library(pheatmap)
   png(file.path(INF,"mr","gsmr","out","gsmr-efo.png"),res=300,width=30,height=15,units="in")
   pheatmap(gsmr_mat,cluster_rows=TRUE,cluster_cols=TRUE,angle_col="315",fontsize_row=24,fontsize_col=24)
   dev.off()
END

# --- HGI

function hgi()
{
  if [ ! -d ${INF}/mr/gsmr/hgi ]; then mkdir -p ${INF}/mr/gsmr/hgi; fi
  export suffix=cis
  awk '$21==ENVIRON["suffix"] {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
  awk -vM=1e6 '{print $2, $3, $4-M, $5+M}' | \
  while read prot chr start end
  do
    export prot=${prot}
    export chr=${chr}
    export start=${start}
    export end=${end}
    echo ${chr} ${start} ${end} > ${INF}/mr/gsmr/hgi/${prot}.bed1
    plink2 --bfile ${INF}/work/INTERVAL --extract bed1 ${INF}/mr/gsmr/hgi/${prot}.bed1 \
           --make-bed --rm-dup force-first list --out ${INF}/mr/gsmr/hgi/${prot}
  done
  awk '$21==ENVIRON["suffix"] {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
  parallel -j15 --env INF --env suffix -C' ' '
    (
      echo -e "SNP A1 A2 freq b se p N"
      cut -f1,2 --complement ${INF}/mr/gsmr/mrx/{2}-${suffix}.mri | \
      awk "{\$2=toupper(\$2);\$3=toupper(\$3);\$7=10^\$7;print}"
    ) | gzip -f > ${INF}/mr/gsmr/hgi/{2}.gz
  '
  export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20210415/results/20210607
  for trait in A2 B1 B2 C2
  do
    export trait=${trait}
    export src=${HGI}/COVID19_HGI_${trait}_ALL_leave_23andme_20210607.b37.txt.gz
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 "{print \$1, \$2, \$3\":\"\$4-M\"-\"\$5+M}" | \
    parallel -C' ' -j15 --env INF --env trait --env suffix '
      (
        echo -e "SNP A1 A2 freq b se p N"
        tabix -h ${src} {3} | awk "
        NR>1{
               if (\$4<\$3) snpid=\"chr\"\$1\":\"\$2\"_\"\$4\"_\"\$3;
               else snpid=\"chr\"\$1\":\"\$2\"_\"\$3\"_\"\$4
               print snpid,\$4,\$3,\$14,\$7,\$8,\$9,\$10+\$11
            }"
      ) | \
      gzip -f> ${INF}/mr/gsmr/trait/${trait}-{2}.gz
    '
  done
}

function gsmr_collect()
# Accepts single argument as GWAS p value used.
{
  export f=$1
  cat *.gsmr | grep -v -e nsnp -e nan | sort -k5,5gr > ${f}.txt

  R --no-save -q <<\ \ END
    f <- Sys.getenv("f")
    gsmr <- read.table(paste0(f,".txt"),col.names=c("prot","trait","b","se","p","nsnp"))
    write.table(within(gsmr,{fdr <- p.adjust(p,method="fdr")}),file=paste0(f,".tsv"),quote=FALSE,row.names=FALSE,sep="\t")
  END
}

## rsid version -- at the expense of removing duplicates

function exposure()
{
  if [ ! -d ${INF}/mr/gsmr/prot ]; then mkdir -p ${INF}/mr/gsmr/prot; fi
  awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
  parallel -j15 --env INF --env suffix -C' ' '
    (
      echo -e "SNP A1 A2 freq b se p N"
      gunzip -c ${INF}/mr/gsmr/mrx/{2}-${suffix}.mrx | \
      sed "1d" | \
      cut -f1,3,4 --complement | \
      awk "a[\$1]++==0{\$7=10^\$7;print}"
    ) | gzip -f > ${INF}/mr/gsmr/prot/{2}-${suffix}.gz
  '
}
