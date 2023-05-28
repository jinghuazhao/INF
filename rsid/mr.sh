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
         left_join(gap.datasets::inf1[c("prot","target.short")],by=c("exposure"="prot")) %>%
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
   pheatmap(mr_mat,cluster_rows=FALSE,cluster_cols=FALSE,angle_col="270",fontsize_row=24,fontsize_col=24)
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
  for id in $(cat <(sed '1d' ${INF}/rsid/efo.txt | grep -v ukb | cut -f4) <(cut -f4 ${INF}/OpenGWAS/ukb-replacement.txt | grep -v finn))
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
        awk -v prot={2} "{\$1=prot;print}" | \        tr " " "\t"
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
  cat <(awk -vFS="\t" 'NR>1 && $4 !~ /ukb/ {print $4,2/(1/$5+1/$6)}' ${INF}/rsid/efo.txt) \
      <(awk '!/finn/' ${INF}/OpenGWAS/ukb-replacement.txt | awk '$5!=""' | awk -vFS="\t" '{print $4,2/(1/$5+1/$6)}') | \
  while read efo N
  do
    export efo=${efo}
    export N=${N}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 "{print \$1, \$2, \$3\":\"\$4-M\"-\"\$5+M}" | \
    parallel -C' ' -j15 --env INF --env efo --env N --env suffix '
      (
        echo -e "SNP A1 A2 freq b se p N"
        bcftools query -f "%CHROM %POS %ID %ALT %REF [%AF] [%ES] [%SE] [%LP] [%SS] \n" -r {3} ${INF}/OpenGWAS/${efo}.vcf.gz | \
        awk -vN=${N} "{
          if(\$4<\$5) snpid=\"chr\"\$1\":\"\$2\"_\"\$4\"_\"\$5;
                 else snpid=\"chr\"\$1\":\"\$2\"_\"\$5\"_\"\$4
          \$9=10^-\$9
          if (\$10==\".\") \$10=N
          print snpid, \$1, \$2, \$4, \$5, \$6, \$7, \$8, \$9, \$10
        }" | sort -k2,2n -k3,3n -k8,8g| cut -d" " -f2,3 --complement | awk "a[\$1]++==0"
      ) | \
      gzip -f> ${INF}/mr/gsmr/trait/${efo}-{2}.gz
    '
  done
  awk -vFS="\t" '/finn/{print $4,2/(1/$5+1/$6)}' ${INF}/OpenGWAS/ukb-replacement.txt | \
  while read efo N
  do
    export efo=${efo}
    export N=${N}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 "{print \$1, \$2, \$3\":\"\$4-M\"-\"\$5+M}" | \
    parallel -C' ' -j15 --env INF --env efo --env N --env suffix '
      export uniprot={1}
      export prot={2}
      export id={3}
      echo {1} {2} {3}
      Rscript -e "
        require(TwoSampleMR)
        require(dplyr)
        efo <- Sys.getenv(\"efo\")
        id <- Sys.getenv(\"id\")
        n <- Sys.getenv(\"N\")
        od <- extract_outcome_data(id,efo) %>%
              distinct() %>%
              mutate(snpid=gap::chr_pos_a1_a2(chr,pos,effect_allele.outcome,other_allele.outcome),
                     effect_allele.outcome=toupper(effect_allele.outcome),
                     other_allele.outcome=toupper(other_allele.outcome)) %>%
              select(snpid,effect_allele.outcome,other_allele.outcome,eaf.outcome,beta.outcome,se.outcome,pval.outcome,samplesize.outcome) %>%
              setNames(c(\"SNP\",\"A1\",\"A2\",\"freq\",\"b\",\"se\",\"p\",\"N\")) %>%
              group_by(SNP) %>%
              slice(which.min(p)) %>%
              data.frame()
        od[is.na(od\$N),\"N\"] <- n
        write.table(od,quote=FALSE,row.names=FALSE)
      " | \
      gzip -f> ${INF}/mr/gsmr/trait/${efo}-{2}.gz
    '
  done
# sbatch --wait ${INF}/rsid/mr.sb
 (
    cat *.gsmr | head -1
    grep -e Exposure -e nan -e TNFB -v *gsmr | tr ':' '\t' | cut -f1 --complement
  ) > gsmr.txt
  (
    grep -e Psoriasis -e HIV -e Anky -e Celiac -e Sarcoidosis -e Sicca -e biliary -e Viral ${INF}/OpenGWAS/efo-update.txt | \
    cut -f1 | \
    grep -f - -v  gsmr.txt
  ) > gsmr-reduce.txt
  awk 'NR>1{print $1,$2}' gsmr.txt | \
  parallel  -C' ' 'sed "1d" ${INF}/mr/gsmr/trait/{1}-{2}-rsid.txt | sort -k7,7g | awk -vprot={1} -vid={2} "NR==1{print prot,id,\$0}"' > gsmr-top.txt
  awk 'NR>1{print $1,$2}' gsmr-reduce.txt | \
  parallel  -C' ' 'sed "1d" ${INF}/mr/gsmr/trait/{1}-{2}-rsid.txt | sort -k7,7g | awk -vprot={1} -vid={2} "NR==1{print prot,id,\$0}"' > gsmr-reduce-top.txt
  R --no-save -q <<\ \ END
    library(dplyr)
    INF <- Sys.getenv("INF")
    metal <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
             filter(cis.trans=="cis") %>%
             select(prot,rsid,Chromosome) %>%
             rename(pqtl=rsid,chr=Chromosome)
    gsmr <- function (input="gsmr.txt",output="gsmr-efo.txt",top="gsmr-top.txt")
    {
      gsmr <- read.delim(input)
      efo <- read.delim(file.path(INF,"OpenGWAS","efo-update.txt")) %>%
             select(id,trait,ncase,ncontrol) %>%
             mutate(Ntotal=ncase+ncontrol) %>%
             filter(id %in% unique(gsmr$Outcome))
      gwas <- read.table(top,col.names=c("prot","id","qtl","a1_qtl","a2_qtl","freq","b_qtl","se_qtl","p_qtl","n_qtl")) %>%
              select(prot,id,qtl,p_qtl)
      gsmr_efo <- left_join(gsmr,pQTLdata::inf1[c("prot","target.short")], by=c("Exposure"="prot")) %>%
                  left_join(filter(metal,prot %in% (gsmr$Exposure)), by=c("Exposure"="prot")) %>%
                  left_join(efo,by=c("Outcome"="id")) %>%
                  left_join(gwas,by=c("Exposure"="prot","Outcome"="id")) %>%
                  rename(protein=Exposure,id=Outcome,Disease=trait,Ncase=ncase,Ncontrol=ncontrol) %>%
                  mutate(protein=target.short,fdr=p.adjust(p,method="fdr")) %>%
                  select(protein,Disease,id,nsnp,fdr,Ncase,Ncontrol,Ntotal,bxy,se,pqtl,p,qtl,p_qtl,chr,pqtl) %>%
                  arrange(fdr)
      write.table(gsmr_efo,output,row.names=FALSE,quote=FALSE,sep="\t")
    }
    gsmr()
    gsmr(input="gsmr-reduce.txt",output="gsmr-efo-reduce.txt",top="gsmr-reduce-top.txt")
    old <- function()
    {
      non_ukb <- read.table(file.path(INF,"OpenGWAS","ukb-replacement.txt"),col.names=c("MRBASEID","x1","x2","New","y1","y2"),sep="\t")
      efo <- read.delim(file.path(INF,"rsid","efo.txt")) %>%
            left_join(non_ukb) %>%
            mutate(MRBASEID=if_else(is.na(New),MRBASEID,New),
                   Ncases=if_else(is.na(y1),Ncases,y1),
                   Ncontrols=if_else(is.na(y2),Ncontrols,y2)) %>%
                   select(-New,-x1,-x2,-y1,-y2)
      gsmr_efo <- gsmr %>%
                  left_join(pQTLdata::inf1[c("prot","target.short")], by=c("Exposure"="prot")) %>%
                  left_join(efo,by=c("Outcome"="MRBASEID")) %>%
                  rename(protein=Exposure,MRBASEID=Outcome) %>%
                  mutate(protein=target.short,fdr=p.adjust(p,method="fdr")) %>%
                  select(protein,MRBASEID,trait,bxy,se,p,nsnp,fdr,Ncases,Ncontrols,id,uri,Zhengetal) %>%
                  arrange(fdr)
      subset(gsmr_efo[setdiff(names(gsmr_efo),c("Zhengetal","uri"))],fdr<=0.05)
      write.table(gsmr_efo,"gsmr-efo.txt",row.names=FALSE,quote=FALSE,sep="\t")
    }
  END
  # LZ
  cut -f1,2,3,8 ${INF}/mr/gsmr/gsmr-efo.txt | \
  awk -vFS="\t" -vOFS='\t' '$4<0.05{ gsub(/-/,".",$1);print $1"-"$2,$3}' | \
  parallel -C '\t' 'echo {2}; ls ${INF}/mr/gsmr/trait/{1}; cp ${INF}/mr/gsmr/trait/{1}/{1}*.pdf ${INF}/mr/gsmr/{1}-{2}.pdf'
  # r(Top SNPs)
  cut -f1,2,3,8 ${INF}/mr/gsmr/gsmr-efo.txt | \
  awk -vFS="\t" -vOFS='\t' '$4<0.05{ gsub(/-/,".",$1);print $1"-"$2,$3}' | \
  parallel -C '\t' 'echo {2}; cp ${INF}/mr/gsmr/out/{1}*.top ${INF}/mr/gsmr/{1}-{2}.top'
  # r(pQTL,GSMR SNPs)
  cut -f1,2,3,8 ${INF}/mr/gsmr/gsmr-efo.txt | \
  awk -vFS="\t" -vOFS='\t' '$4<0.05{ gsub(/-/,".",$1);print $1"-"$2,$3}' | \
  parallel -C '\t' 'echo {2}; cp ${INF}/mr/gsmr/out/{1}*.r ${INF}/mr/gsmr/{1}-{2}.r'
}
# zgrep SAMPLE ${INF}/OpenGWAS/*.vcf.gz | cut -f1,7,8 > ${INF}/OpenGWAS/ieu.sample
# sed 's/.vcf.gz:/\t/;s/TotalControls=/\t/;s/,TotalCases=/\t/;s/,StudyType/\t/' ieu.sample | cut -f1,3,4 > ${INF}/OpenGWAS/ieu.N
#   zgrep -h SAMPLE ${INF}/OpenGWAS/${efo}.vcf.gz | cut -f1,7,8 > ${INF}/mr/gsmr/trait/${efo}.sample

# https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics
export TMPDIR=${HPC_WORK}/work
export b38tob37=~/hpc-work/bin/hg38ToHg19.over.chain.gz

function harmonise()
{
  export hdata=${1}
  export GCST=$(echo ${1} | sed 's|/|\t|g' | cut -f10)
  export fields=${2}
  echo ${hdata} $GCST ${fields}
  case ${GCST} in
  GCST90014023)
    if [ ! -f ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}-b37.bed ]; then
       cat <(echo \#chrom Start End SNPid | tr ' ' '\t') \
           <(gunzip -c ${hdata} | awk -vOFS="\t" 'NR>1&&$3!="NA"{print "chr"$3,$4,$4+1,$1}') > ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.bed
       liftOver -bedPlus=4 \
                ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.bed \
                ${b38tob37} \
                ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}-b37.bed \
                ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}-b37.unlifted.bed
    fi
    cat <(echo chr pos $(gunzip -c ${hdata} | head -1 | cut -f${fields} | cut -f1 --complement)) \
        <(join -13 -21 <(awk '!/Un/ && !/_random/' ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}-b37.bed | sed '1d;s/chr//' | cut -f1,2,4 | sort -k3,3) \
                       <(gunzip -c ${hdata} | sed '1d' | cut -f${fields} | sort -k1,1) | \
          cut -d' ' -f1 --complement | awk '!/X|Y/' | sort -k1,1n -k2,2n) | tr ' ' '\t' | \
    bgzip -f > ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.tsv.gz
    tabix -S1 -s1 -b2 -e2 -f ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.tsv.gz
    ;;
  GCST90061440 | GCST90019016 | GCST90014325)
    gunzip -c ${hdata} | head -3 | cut -f${fields}
    cat <(gunzip -c ${hdata} | head -1 | cut -f${fields}) \
        <(gunzip -c ${hdata} | sed '1d' | cut -f${fields} | awk '$1!="NA" && !/X|Y/' | sort -k3,3 -k4,4n) | \
    bgzip -f > ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.tsv.gz
    tabix -S1 -s3 -b4 -e4 -f ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.tsv.gz
    ;;
  *)
    ;;
  esac
}
for h in ${INF}/OpenGWAS/GCST90061440/harmonised/34033851-GCST90061440-EFO_1001486.h.tsv.gz \
         ${INF}/OpenGWAS/GCST90019016/harmonised/34927100-GCST90019016-EFO_0000676.h.tsv.gz \
         ${INF}/OpenGWAS/GCST90014023/harmonised/34012112-GCST90014023-EFO_0001359.h.tsv.gz \
         ${INF}/OpenGWAS/GCST90014325/harmonised/34103634-GCST90014325-EFO_0000270.h.tsv.gz
do
 export f=$( gunzip -c ${h} | head -1 | tr '\t' '\n' | \
 awk '/hm_variant_id|hm_chrom|hm_pos|hm_rsid|hm_other_allele|hm_effect_allele|hm_beta|hm_effect_allele_frequency|p_value|standard_error/{printf ","NR}'|\
 sed 's/^[,]//')
 harmonise ${h} ${f}
done

function GCST90010715()
{
  gunzip -c ${INF}//OpenGWAS/jia/jia_gwas_snptest_control_freq.txt.gz | \
  sed '1d' | \
  awk '
  {
    if ($4<$5) snpid = "chr" $2 ":" $3 "_" $4 "_" $5; else snpid = "chr" $2 ":" $3 "_" $5 "_" $4
    g=$6+$7+$8
    EAF=(0.5*$6+$7)/g
    print snpid, $1, EAF
  }' | \
  sort -k1,1 > ${INF}/OpenGWAS/jia/jia.eaf
  export efo=GCST90010715
  cat <(echo snpid rsid chr pos a1 a2 freq b se p N) \
      <(sed '1d' ${INF}/OpenGWAS/${efo}/${efo}_buildGRCh37.tsv | \
        awk '
        {
          chr=$3+0; pos=$4; a1=$6; a2=$5
          if (a1<a2) snpid="chr" chr ":" pos "_" a1 "_" a2; else snpid="chr" chr ":" pos "_" a2 "_" a1
          print snpid, $2, chr, pos, a1, a2, "NA", $12, $13, $11, 2/(1/3305+1/9196)
        }' | \
        sort -k1,1 -k9,9g | \
        awk 'a[$1]++==0' | \
        sort -k1,1 | \
        join - <(cut -d' ' -f1,3 ${INF}/OpenGWAS/jia/jia.eaf) | \
        sort -k3,3n -k4,4n | \
        awk '{$7=$12;print}' | \
        cut -d' ' -f12 --complement) | \
  tr ' ' '\t' | \
  bgzip -f > ${INF}/OpenGWAS/${efo}/harmonised/${efo}.tsv.gz
  tabix -S1 -s3 -b4 -e4 -f ${INF}/OpenGWAS/${efo}/harmonised/${efo}.tsv.gz
}
# gunzip -c jia/jia_gwas_snptest_control_freq.txt.gz | head -1 | tr ' ' '\n' | awk '{print "#"NR,$1}'
#1 rsid
#2 chromosome
#3 position
#4 alleleA
#5 alleleB
#6 controls_AA
#7 controls_AB
#8 controls_BB
#9 controls_NULL
#10 controls_total
#11 controls_maf

# alternate_ids	variant_id	chromosome	position	alleleA	alleleB	all_maf	all_OR	all_OR_lower	all_OR_upper	p_value	frequentist_add_beta_1	frequentist_add_se_1
# rs2977608	rs2977608	01	768253	A	C	0.242297	1.03793	0.971733	1.10864	0.104282	0.071745	0.044166

function GCST90134602()
{
  export GCST=GCST90134602
  export HGI=${INF}/OpenGWAS/${GCST}/${GCST}_buildGRCh38.tsv
  cat <(echo \#chrom Start End SNPid | tr ' ' '\t') \
      <(awk -vOFS="\t" 'NR>1{gsub(/23/,"X",$1);print "chr"$1,$2,$2+1,$5}' ${HGI}) > ${INF}/OpenGWAS/${GCST}/${GCST}.bed
  liftOver -bedPlus=4 \
           ${INF}/OpenGWAS/${GCST}/${GCST}.bed \
           ${b38tob37} \
           ${INF}/OpenGWAS/${GCST}/${GCST}-b37.bed \
           ${INF}/OpenGWAS/${GCST}/${GCST}-b37.unlifted.bed

  cat <(echo chr pos A2 A1 b se p cases controls N freq rsid) \
      <(join -j3 <(sed '1d;s/chr//' ${INF}/OpenGWAS/${GCST}/${GCST}-b37.bed | cut -f1,2,4 | sort -k3,3) \
                 <(sed '1d' ${HGI} | cut -f3-5,7,8,9,10-12,14,15 | sort -k3,3) | \
        sort -k1,1 -k8,8g | \
        awk 'a[$1]++==0' | \
        cut -d' ' -f1 --complement | \
        awk '{gsub(/X/,"23",$1);print}' | \
        sort -k1,1n -k2,2n) | \
  tr ' ' '\t' | \
  bgzip -f > ${INF}/OpenGWAS/${GCST}/${GCST}.tsv.gz
  tabix -S1 -s1 -b2 -e2 -f ${GCST}.tsv.gz

  gzip -f ${INF}/OpenGWAS/${GCST}/${GCST}.bed \
          ${INF}/OpenGWAS/${GCST}/${GCST}-b37.unlifted.bed
}
# chromosome      base_pair_location      other_allele    effect_allele   SNP     all_meta_N      beta    standard_error  p_value all_inv_var_meta_cases  meta_controls       all_inv_var_meta_effective      all_inv_var_het_p       all_meta_AF     variant_id
# 1       758351  A       G       1:758351:A:G    12      -2.5351e-02     3.3686e-02      4.5172e-01      7077    51575   5827    4.1912e-01      1.249e-01  rs12238997
#
# gunzip -c OpenGWAS/GCST90134602/GCST90134602.tsv.gz | head
# chr     pos     A2      A1      b       se      p       cases   controls        N       freq    rsid
# 1       714596  T       C       -3.3903e-02     7.0748e-02      6.3179e-01      5743    41335   4676    3.347e-02       rs149887893

function HanY()
{
  gunzip -c ${INF}/OpenGWAS/GCST010043/HanY_prePMID_asthma_Meta-analysis_UKBB_TAGC.txt.gz | \
  awk -vOFS='\t' '
    {
      if (NR==1) print "snpid chr pos a1 a2 z p N";
      else
         {
            chr=$2;pos=$3;a1=toupper($4);a2=toupper($5)
            if (a1<a2) snpid="chr"chr":"pos"_"a1"_"a2;
                  else snpid="chr"chr":"pos"_"a2"_"a1
            print snpid, chr, pos, a1, a2, $6, $7, $10
         }
    }' | \
  sort -k1,1 -k7,7g | \
  awk 'a[$1]++==0' | \
  sort -k2,2n -k3,3n | \
  bgzip -f > ${INF}/OpenGWAS/GCST010043/GCST010043.tsv.gz
  tabix -S1 -s2 -b3 -e3 ${INF}/OpenGWAS/GCST010043/GCST010043.tsv.gz
}
# GCST010043/HanY_prePMID_asthma_Meta-analysis_UKBB_TAGC.txt.gz
# SNP     CHR     BP      EA      NEA     Z       P       Direction_UKBB_TAGC     P_het   N
# 1:100004463_TA_T        1       100004463       T       TA      0.374   0.7081  -+      0.2873  536345
# GCST010042/HanY_prePMID_asthma_UKBB.txt.gz
# SNP     CHR     BP      EA      NEA     EAF     INFO    OR      OR_95L  OR_95U  P       N
# 1:692794_CA_C   1       692794  CA      C       0.881938        0.824483        1.00032752742223        0.980323146514236       1.0207401158248 0.97    393859

function mr_rsid()
# data with RSid for TwoSampleMR
{
  for r in {1..59}
  do
    export region=$(awk -vr=${r} 'NR==r{print $4":"$5"-"$6}' ${INF}/TNFB/cis.dat)
    export prot=$(awk -vr=${r} 'NR==r{print $3}' ${INF}/TNFB/cis.dat)
    echo ${prot} ${region}
    if [ ${prot} == "TNFB" ]; then continue; fi
    for OpenGWAS in $(sed '1d' ${INF}/OpenGWAS/efo-update.txt | cut -f1 | awk '/ieu|ebi|bbj/')
    do
      export N=$(grep -w ${OpenGWAS} ${INF}/OpenGWAS/efo-update.txt | awk -vFS='\t' '{print 2/(1/$3+1/$4)}')
      (
        echo -e "SNP A1 A2 freq b se p N"
        bcftools query -f "%ID %ALT %REF [%AF] [%ES] [%SE] [%LP] [%SS]\n" -r ${region} ${INF}/OpenGWAS/${OpenGWAS}.vcf.gz | \
        awk -vN=${N} '{$7=10^-$7;if ($8==".") $8=N;print}'
      ) > ${INF}/mr/gsmr/trait/${prot}-${OpenGWAS}-rsid.txt
    done
    for GCST in $(sed '1d' ${INF}/OpenGWAS/efo-update.txt | cut -f1 | awk '/^GCST/')
    do
      export N=$(grep -w ${GCST} ${INF}/OpenGWAS/efo-update.txt | awk -vFS='\t' '{print 2/(1/$3+1/$4)}')
      (
        echo -e "SNP A1 A2 freq b se p N"
        case ${GCST} in
        GCST90010715)
          tabix ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.tsv.gz ${region} | \
          awk '{print $2,$5,$6,$7,$8,$9,$10,$11}'
          ;;
        GCST90014325 | GCST90061440)
          tabix ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.tsv.gz ${region} | \
          awk -vN=${N} '{print $2,$6,$5,$8,$7,$10,$9,N}'
          ;;
        GCST90019016)
          tabix ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.tsv.gz ${region} | \
          awk -vN=${N} '{print $2,$6,$5,$8,$7,$9,$10,N}'
          ;;
        GCST90014023)
          tabix ${INF}/OpenGWAS/${GCST}/harmonised/${GCST}.tsv.gz ${region} | \
          awk -vN=${N} '{print $3,$7,$6,$9,$8,$11,$10,N}'
          ;;
        GCST90134602)
          tabix ${INF}/OpenGWAS/${GCST}/${GCST}.tsv.gz ${region} | \
          awk '{print $12,$4,$3,$11,$5,$6,$7,$10}'
          ;;
        *)
          ;;
        esac
      ) > ${INF}/mr/gsmr/trait/${prot}-${GCST}-rsid.txt
    done
    for finngen in $(sed '1d' ${INF}/OpenGWAS/efo-update.txt | cut -f1 | awk -vFS='\t' '/^finngen/')
    do
      export N=$(grep -w ${finngen} ${INF}/OpenGWAS/efo-update.txt | awk -vFS='\t' '{print 2/(1/$3+1/$4)}')
      (
         echo -e "SNP A1 A2 freq b se p N"
         tabix ${INF}/OpenGWAS/${finngen}.gz ${region} | \
         awk -vN=${N} '{print $5,$4,$3,$11,$9,$10,$7,N}'
      ) > ${INF}/mr/gsmr/trait/${prot}-${finngen}-rsid.txt
    done
    for efo in $(sed '1d' ${INF}/OpenGWAS/efo-update.txt | cut -f1)
    do
       sed '1d' ${INF}/mr/gsmr/trait/${prot}-${efo}-rsid.txt | \
       sort -k7,7g | \
       head -1 > ${INF}/mr/gsmr/trait/${prot}-${efo}-top.txt
       if [ $(wc -l ${INF}/mr/gsmr/trait/${prot}-${efo}-top.txt | cut -d' ' -f1) -eq 0 ]; then rm ${INF}/mr/gsmr/trait/${prot}-${efo}-top.txt; fi
    done
  done
}
# gunzip -c GCST90134602.tsv.gz | head -1 | tr '\t' '\n' | awk '{print NR,$1}'

#!/usr/bin/bash

#SBATCH --account CARDIO-SL0-CPU
#SBATCH --partition cardio
#SBATCH --qos=cardio
#SBATCH --mem=28800

##SBATCH --account=PETERS-SL3-CPU
##SBATCH --partition=cclake-himem
#SBATCH --job-name=_trait
##SBATCH --mem=6840
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/mr/gsmr/slurm/_trait_%A_%a.o
#SBATCH --error=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/mr/gsmr/slurm/_trait_%A_%a.e

#SBATCH --time=12:00:00

function efo_update()
# switch to gsmr_trait in mr.sb apeared problematic
{
  export EFO_UPDATE=${INF}/OpenGWAS/efo-update.txt
  sed '1d' ${EFO_UPDATE} | grep -e ebi -e ieu | awk -vFS="\t" '{print $1,2/(1/$3+1/$4)}' | \
  while read efo N
  do
    export efo=${efo}
    export N=${N}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | awk '!/TNFB/' | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 "{print \$1, \$2, \$3\":\"\$4-M\"-\"\$5+M}" | \
    parallel -C' ' -j15 --env INF --env efo --env N --env suffix '
      (
        echo -e "SNP A1 A2 freq b se p N"
        bcftools query -f "%CHROM %POS %ID %ALT %REF [%AF] [%ES] [%SE] [%LP] [%SS] \n" -r {3} ${INF}/OpenGWAS/${efo}.vcf.gz | \
        awk -vN=${N} "{
               if(\$4<\$5) snpid=\"chr\"\$1\":\"\$2\"_\"\$4\"_\"\$5;
               else snpid=\"chr\"\$1\":\"\$2\"_\"\$5\"_\"\$4
               \$9=10^-\$9
               if (\$10==\".\") \$10=N
               print snpid, \$1, \$2, \$4, \$5, \$6, \$7, \$8, \$9, \$10
             }" | sort -k2,2n -k3,3n -k9,9gr | cut -d" " -f2,3 --complement | awk "a[\$1]++==0"
      ) | \
      gzip -f > ${INF}/mr/gsmr/trait/{2}-${efo}.gz
    '
  done
  sed '1d' ${EFO_UPDATE} | grep -v -e ebi -e ieu -e finn | awk -vFS="\t" '{print $1,2/(1/$3+1/$4)}' | \
  while read efo N
  do
    export efo=${efo}
    export N=${N}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | awk '!/TNFB/' | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 "{print \$1, \$2, \$3\":\"\$4-M\"-\"\$5+M}" | \
    parallel -C' ' -j15 --env INF --env efo --env cases --env controls '
      echo ${efo} {2}
      (
        echo SNP A1 A2 freq b se p N
        case ${efo} in
        GCST90019016)
           tabix ${INF}/OpenGWAS/${efo}/harmonised/${efo}.tsv.gz {3} | \
           awk -vefo=${efo} -vN=${N} "
           {
               chr=\$3; pos=\$4; a1=\$6; a2=\$5
               if (a1<a2) snpid=\"chr\" chr \":\" pos \"_\" a1 \"_\" a2; else snpid=\"chr\" chr \":\" pos \"_\" a2 \"_\" a1
               print snpid, a1, a2, \$8, \$7, \$9, \$10, N
           }" | \
           sort -k1,1 -k7,7g | \
           awk "a[\$1]++==0"
           ;;
        GCST90010715)
           tabix ${INF}/OpenGWAS/${efo}/harmonised/${efo}.tsv.gz {3} | \
           awk -vefo=${efo} -vN=${N} "
           {
               chr=\$3; pos=\$4; a1=\$5; a2=\$6
               if (a1<a2) snpid=\"chr\" chr \":\" pos \"_\" a1 \"_\" a2; else snpid=\"chr\" chr \":\" pos \"_\" a2 \"_\" a1
               print snpid, a1, a2, \$7, \$8, \$9, \$10, \$11
           }" | \
           sort -k1,1 -k7,7g | \
           awk "a[\$1]++==0"
           ;;
        GCST90061440 | GCST90014325)
           tabix ${INF}/OpenGWAS/${efo}/harmonised/${efo}.tsv.gz {3} | \
           awk -vefo=${efo} -vN=${N} "
           {
               chr=\$3; pos=\$4; a1=\$6; a2=\$5
               if (a1<a2) snpid=\"chr\" chr \":\" pos \"_\" a1 \"_\" a2; else snpid=\"chr\" chr \":\" pos \"_\" a2 \"_\" a1
               print snpid, a1, a2, \$8, \$7, \$10, \$9, N
           }" | \
           sort -k1,1 -k7,7g | \
           awk "a[\$1]++==0"
           ;;
        GCST90014023)
           tabix ${INF}/OpenGWAS/${efo}/harmonised/${efo}.tsv.gz {3} | \
           awk -vefo=${efo} -vN=${N} "
           {
               chr=\$1; pos=\$2; a1=\$7; a2=\$6
               if (a1<a2) snpid=\"chr\" chr \":\" pos \"_\" a1 \"_\" a2; else snpid=\"chr\" chr \":\" pos \"_\" a2 \"_\" a1
               print snpid, \$7, \$6, \$9, \$8, \$11, \$10, N
           }" | \
           sort -k1,1 -k7,7g | \
           awk "a[\$1]++==0"
           ;;
        GCST90134602)
           tabix ${INF}/OpenGWAS/${efo}/${efo}.tsv.gz {3} | \
           awk "{
               chr=\$1;pos=\$2;a1=\$4;a2=\$3;
               if (a1<a2) snpid=\"chr\" chr \":\" pos \"_\" a1 \"_\" a2; else snpid=\"chr\" chr \":\" pos \"_\" a2 \"_\" a1
               print snpid, a1, a2, \$11, \$5, \$6, \$7, \$10
           }"
        ;;
      # chr pos A2 A1 b se p cases controls N freq rsid
        *)
        ;;
        esac
      ) | \
      gzip -f > ${INF}/mr/gsmr/trait/{2}-${efo}.gz
    '
  done
  sed '1d' ${EFO_UPDATE} | grep finn | awk -vFS="\t" '{print $1,2/(1/$3+1/$4)}' | \
  while read efo N
  do
    export efo=${efo}
    export N=${N}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 "{print \$1, \$2, \$3\":\"\$4-M\"-\"\$5+M}" | \
    parallel -C' ' -j15 --env INF --env efo --env N '
      echo ${efo} {2} {3}
      cat <(echo SNP A1 A2 freq b se p N) \
          <(tabix ${INF}/OpenGWAS/${efo}.gz {3} | \
            awk -vN=${N} "{
               chr=\$1;pos=\$2;a1=toupper(\$4);a2=toupper(\$3)
               if (a1<a2) snpid=\"chr\" chr \":\" pos \"_\" a1 \"_\" a2;
                     else snpid=\"chr\" chr \":\" pos \"_\" a2 \"_\" a1
               print snpid,a1,a2,\$11,\$9,\$10,\$7,N
               }" | \
            sort -k1,1 -k7,7g | \
            awk "a[\$1]++==0"
           ) | \
      gzip -f > ${INF}/mr/gsmr/trait/{2}-${efo}.gz
    '
  done
# chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls
  sed '1d' ${EFO_UPDATE} | awk -vFS="\t" '{print $1,2/(1/$3+1/$4)}' | \
  while read efo N
  do
    export efo=${efo}
    export N=${N}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    parallel -C' ' -j15 --env INF --env efo '
      echo ${efo} ${INF}/mr/gsmr/trait/{2}-${efo}.gz > ${INF}/mr/gsmr/trait/gsmr_{2}-${efo}
    '
  done
}

function gsmr_mr_heatmap()
{
#!/usr/bin/bash

#SBATCH --account=PETERS-SL3-CPU
#SBATCH --partition=cclake-himem
#SBATCH --job-name=_mr
#SBATCH --mem=6840
#SBATCH --time=12:00:00
#SBATCH --output=_mr_%A_%a.o
#SBATCH --error=_mr_%A_%a.e

Rscript -e '
   INF <- Sys.getenv("INF")
   suppressMessages(library(dplyr))
   library(stringr)
   gsmr <- read.delim(file.path(INF,"mr","gsmr","gsmr-efo.txt")) %>%
           mutate(outcome=paste0(Disease),
                  exposure=protein,
                  z=bxy/se,
                  or=exp(bxy)
           ) %>%
           left_join(gap.datasets::inf1[c("target.short","gene")],by=c("exposure"="target.short")) %>%
           select(gene,outcome,z,or,bxy,se,p,nsnp,fdr)
   gene <- unique(with(gsmr,gene))
   outcome <- unique(with(gsmr,outcome))
   n <- length(gene)
   m <- length(outcome)
   gsmr_mat <- matrix(NA,m,n)
   colnames(gsmr_mat) <- gene
   rownames(gsmr_mat) <- outcome
   gsmr_mat_fdr <- gsmr_mat
   for(k in 1:nrow(gsmr))
   {
      t <- gsmr[k,c("gene","outcome","or","fdr")]
      i <- t[["outcome"]]
      j <- t[["gene"]]
      gsmr_mat[i,j] <- t[["or"]]
      gsmr_mat_fdr[i,j] <- t[["fdr"]]
   }
   rownames(gsmr_mat) <- gsub("\\b(^[a-z])","\\U\\1",rownames(gsmr_mat),perl=TRUE)
   z <- rowSums(is.na(gsmr_mat))
   gsmr_mat_reduce <- gsmr_mat[z<40,]
   gsmr_mat_reduce_fdr <- gsmr_mat_fdr[z<40,]
   rm(gene,outcome)
   options(width=200)
   subset(gsmr,fdr<=0.05)
   library(grid)
   library(pheatmap)
   png(file.path(INF,"mr","gsmr","gsmr-efo.png"),res=300,width=30,height=18,units="in")
   setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))),
           action="prepend")
   pheatmap(gsmr_mat,cluster_rows=FALSE,cluster_cols=FALSE,angle_col="315",fontsize_row=30,fontsize_col=30,
            display_numbers = matrix(ifelse(!is.na(gsmr_mat) & abs(gsmr_mat_fdr) <= 0.01, "*", ""), nrow(gsmr_mat)), fontsize_number=20)
   setHook("grid.newpage", NULL, "replace")
   grid.text("Proteins", y=-0.07, gp=gpar(fontsize=48))
   grid.text("Immune-mediated outcomes", x=-0.07, rot=90, gp=gpar(fontsize=48))
   dev.off()
   png(file.path(INF,"mr","gsmr","gsmr-efo-reduce.png"),res=300,width=30,height=18,units="in")
   setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))),
           action="prepend")
   pheatmap(gsmr_mat_reduce,cluster_rows=FALSE,cluster_cols=FALSE,angle_col="315",fontsize_row=30,fontsize_col=30,
            display_numbers = matrix(ifelse(!is.na(gsmr_mat_reduce) & abs(gsmr_mat_reduce_fdr) <= 0.02, "*", ""), nrow(gsmr_mat_reduce)),
            fontsize_number=20)
   setHook("grid.newpage", NULL, "replace")
   grid.text("Proteins", y=-0.07, gp=gpar(fontsize=48))
   grid.text("Immune-mediated outcomes", x=-0.07, rot=90, gp=gpar(fontsize=48))
   dev.off()
   write.table(colnames(gsmr_mat),quote=FALSE,row.names=FALSE)
   others <- function()
   {
   # moved in the latest decision
     tnfb <- filter(gsmr,gene=="LTA" & fdr<=0.05) %>% rename(Effect=bxy,StdErr=se)
     attach(tnfb)
     png(file.path(INF,"mr","gsmr","out","TNFB.png"),height=10,width=18,units="in",res=300)
     requireNamespace("meta")
     mg <- meta::metagen(Effect,StdErr,sprintf("%s",gsub("IGA","IgA",gsub("\\b(^[a-z])","\\U\\1",outcome,perl=TRUE))),sm="OR",title="TNFB")
     meta::forest(mg,colgap.forest.left = "0.5cm",fontsize=24,
                  leftcols=c("studlab"),leftlabs=c("Outcome"),
                  rightcols=c("effect","ci","pval"),rightlabs=c("OR","95% CI","GSMR P"),digits=2,digits.pval=2,scientific.pval=TRUE,
                  plotwidth="5inch",sortvar=Effect,
                  common=FALSE, random=FALSE, print.I2=FALSE, print.pval.Q=FALSE, print.tau2=FALSE,addrow=TRUE,backtransf=TRUE,spacing=1.6)
     with(mg,cat("prot =", p, "MarkerName =", m, "Q =", Q, "df =", df.Q, "p =", pval.Q,
                 "I2 =", I2, "lower.I2 =", lower.I2, "upper.I2 =", upper.I2, "\n"))
     dev.off()
     requireNamespace("rmeta")
     gap::ESplot(data.frame(id=outcome,b=Effect,se=StdErr),fontsize=22)
     rmeta::metaplot(Effect,StdErr,
                     labels=outcome,
                     xlab="Effect size",ylab="",xlim=c(-1.5,1),cex=2,
                     summlabel="",summn=NA,sumse=NA,sumn=NA,lwd=2,boxsize=0.6,
                     colors=rmeta::meta.colors(box="red",lines="blue", zero="black", summary="red", text="black"))
     detach(tnfb)
   }
'

Rscripe -e '
   INF <- Sys.getenv("INF")
   suppressMessages(library(dplyr))
   library(stringr)
   mr <- read.delim(file.path(INF,"mr","gsmr","mr-efo-mr.tsv")) %>%
         mutate(disease=gsub("\\s\\(oligoarticular or rheumatoid factor-negative polyarticular\\)","",disease)) %>%
         mutate(disease=gsub("\\s\\(non-Lofgren's syndrome\\)","",disease)) %>%
         mutate(disease=gsub("_PSORIASIS","",disease)) %>%
         mutate(outcome=disease,
                exposure=gene,
                or=exp(b),
                group=cut(or,breaks=c(0,0.49,0.99,1.49,2))) %>%
         select(gene,outcome,or,b,se,pval,nsnp,fdr,group) %>%
         mutate(or=ifelse(!is.na(or) & or<=2,or,NA))
   options(width=200)
   subset(mr,fdr<=0.05)
   gene <- unique(with(mr,gene))
   outcome <- unique(with(mr,outcome))
   n <- length(gene)
   m <- length(outcome)
   mr_mat <- matrix(NA,m,n)
   colnames(mr_mat) <- gene
   rownames(mr_mat) <- outcome
   mr_mat_fdr <- mr_mat
   for(k in 1:nrow(mr))
   {
      t <- mr[k,c("gene","outcome","or","group","fdr")]
      i <- t[["outcome"]]
      j <- t[["gene"]]
      mr_mat[i,j] <- t[["or"]]
      mr_mat_fdr[i,j] <- t[["fdr"]]
   }
   rownames(mr_mat) <- gsub("\\b(^[a-z])","\\U\\1",rownames(mr_mat),perl=TRUE)
   rm(gene,outcome)
   library(grid)
   library(pheatmap)
   png(file.path(INF,"mr","gsmr","mr-efo.png"),res=300,width=30,height=18,units="in")
   setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))),
           action="prepend")
   pheatmap(mr_mat,cluster_rows=FALSE,cluster_cols=FALSE,angle_col="315",fontsize_row=30,fontsize_col=30,
            display_numbers = matrix(ifelse(!is.na(mr_mat) & abs(mr_mat_fdr) <= 0.05, "*", ""), nrow(mr_mat)), fontsize_number=20)
   setHook("grid.newpage", NULL, "replace")
   grid.text("Proteins", y=-0.07, gp=gpar(fontsize=48))
   grid.text("Immune-mediated outcomes", x=-0.07, rot=90, gp=gpar(fontsize=48))
   dev.off()
   write.table(colnames(mr_mat),quote=FALSE,row.names=FALSE)
'
}

function mr_recollect()
{
  export rt=${INF}/mr/gsmr
  (
    cat ${rt}/mr/*mr | head -1
    ls ${rt}/mr/*mr | grep -v TNFB | xargs -l -I {} sed '1d' {}
  ) > $rt/mr-efo.mr
  (
    cat ${rt}/mr/*het | head -1
    ls ${rt}/mr/*het | grep -v TNFB | xargs -l -I {} sed '1d' {}
  ) > $rt/mr-efo.het
  (
    cat ${rt}/mr/*single | head -1
    ls ${rt}/mr/*single | grep -v TNFB | xargs -l -I {} sed '1d' {}
  ) > $rt/mr-efo.single
  Rscript -e '
    options(width=200)
    library(dplyr)
    library(openxlsx)
    rt <- Sys.getenv("rt")
    genes <- c("IL12B","CD40","IL18R1","CD6","CD5","CXCL5")
    IVW <- read.delim(file.path(rt,"mr-efo.mr")) %>%
           mutate(fdr=p.adjust(pval,method="fdr")) %>%
           left_join(select(pQTLdata::inf1,prot,gene)) %>%
           select(gene,id.outcome,b,se,pval,fdr,nsnp,disease)
    cat(nrow(filter(IVW,fdr<=0.05)),"\n")
    write.table(IVW,file.path(rt,"mr-efo-mr.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
    Heterogeneity <- read.delim(file.path(rt,"mr-efo.het")) %>%
                     mutate(fdr=p.adjust(Q_pval,method="fdr")) %>%
                     left_join(select(pQTLdata::inf1,prot,gene)) %>%
                     select(gene,id.outcome,method,Q,Q_df,Q_pval,fdr,disease)
    cat(nrow(filter(Heterogeneity,fdr<=0.05)),"\n")
    write.table(Heterogeneity,file.path(rt,"mr-efo-het.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
    Single <- read.delim(file.path(rt,"mr-efo.single")) %>%
             left_join(select(pQTLdata::inf1,prot,gene)) %>%
             select(gene,id.outcome,SNP,b,se,p)
#   GSMR <- read.delim(file.path(rt,"out","gsmr-efo.txt")) %>%
#           rename(target.short=protein) %>%
#           left_join(select(pQTLdata::inf1,target.short,gene)) %>%
#           select(gene,MRBASEID,trait,bxy,se,p,nsnp,fdr,Ncases,Ncontrols,id,uri,Zhengetal)
    GSMR <- read.delim(file.path(rt,"gsmr-efo.txt")) %>%
            rename(target.short=protein) %>%
            left_join(select(pQTLdata::inf1,target.short,gene))
    GSMR_genes <- filter(GSMR,fdr<=0.05) %>% pull(gene) %>% unique
    IVW_genes <- filter(IVW,gene%in%genes)
    Heterogeneity_genes <- filter(Heterogeneity,gene%in%genes)
    Single_genes <- filter(Single,gene%in%genes)
    xlsx <- file.path(rt,"gsmr-mr.xlsx")
    wb <- createWorkbook(xlsx)
    hs <- createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=12, fontName="Arial Narrow", fgFill="#4F80BD")
    for (sheet in c("GSMR","IVW","Heterogeneity","Single","IVW_genes","Heterogeneity_genes","Single_genes"))
    {
      addWorksheet(wb,sheet,zoom=150)
      writeData(wb,sheet,sheet,xy=c(1,1),headerStyle=createStyle(textDecoration="BOLD",
                fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
      body <- get(sheet)
      writeDataTable(wb, sheet, body, xy=c(1,2), headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
      freezePane(wb, sheet, firstCol=TRUE, firstActiveRow=3)
      width_vec <- apply(body, 2, function(x) max(nchar(as.character(x))+2, na.rm=TRUE))
    # width_vec_header <- nchar(colnames(body))+2
      setColWidths(wb, sheet, cols = 1:ncol(body), widths = width_vec)
      writeData(wb, sheet, tail(body,1), xy=c(1, nrow(body)+2), colNames=FALSE, borders="rows", borderStyle="thick")
    }
    saveWorkbook(wb, file=xlsx, overwrite=TRUE)
  '
}

# --- HGI

function hgi()
{
  if [ ! -d ${INF}/mr/gsmr/hgi ]; then mkdir -p ${INF}/mr/gsmr/hgi; fi
  export suffix=cis
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
               print snpid,\$4,\$3,\$14,\$7,\$8,\$9,2/(1/\$10+1/\$11)
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

--- specific runs ---

for dir in out trait; do if [ ! -d ${INF}/mr/GWAS/${dir} ]; then mkdir -p ${INF}/mr/GWAS/${dir}; fi; done
export M=1e6
function RA()
# No gain from this data so remove the sumstats
{
  cut -f3 ${INF}/work/INF1.METAL | sed '1d' | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
  parallel -C' ' --env INF --env M '
  (
    echo SNP A1 A2 freq b se p N
    gunzip -c ${INF}/OpenGWAS/Ha_et_al_2020_RA_Trans.txt.gz | \
    awk -v chr={3} -v start={4} -v end={5} -v M=${M} "
    {
      if (NR>1)
      {
        split(toupper(\$1),a,\":\");
        a1=toupper(\$2)
        a2=toupper(\$3)
        if(a1<a2) snpid=\"chr\" chr \":\" a[2] \"_\" a1 \"_\" a2;
        else snpid=\"chr\" chr \":\" a[2] \"_\" a2 \"_\" a1;
        if(a[1]==chr && a[2] >= start-M && a[2] < end+M) print a[1], a[2], snpid, a1, a2, \$4, \$5, \$6, \$7, \$8
      }
    }" | sort -k1,1n -k2,2n | cut -d" " -f1-2 --complement
  ) | gzip -f > ${INF}/mr/GWAS/trait/${trait}-{2}.gz
'
}
export trait=RA
RA
function ukb_b_9125()
{
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 "{print \$1, \$2, \$3\":\"\$4-M\"-\"\$5+M}" | \
    parallel -C' ' -j15 --env INF --env efo --env N --env suffix '
      (
        echo -e "SNP A1 A2 freq b se p N"
        bcftools query -f "%CHROM %POS %ID %ALT %REF [%AF] [%ES] [%SE] [%LP] [$N] \n" -r {3} ${INF}/OpenGWAS/${trait}.vcf.gz | \
        awk "{
          if(\$4<\$5) snpid=\"chr\"\$1\":\"\$2\"_\"\$4\"_\"\$5;
                 else snpid=\"chr\"\$1\":\"\$2\"_\"\$5\"_\"\$4
          \$9=10^-\$9
          print snpid, \$4, \$5, \$6, \$7, \$8, \$9, \$10
        }"
      ) | \
      gzip -f> ${INF}/mr/GWAS/trait/${trait}-{2}.gz
    '
}
# bcftools query -r
# 6 42738749 rs9381218 C T 0.455975 9.15774e-05 0.000226441 0.161151 10285.1
export trait=ukb-b-9125
export N=$(awk 'BEGIN{print 2/(1/457732+1/5201)}')
ukb_b_9125
# Q9P0M4 IL.17C 16 88704999 88706881
# MarkerName      Allele1 Allele2 Freq1   Effect  StdErr  P-value TotalSampleSize
# 5:29439275      t       c       0.4862  -0.0085 0.0145  0.5596  291180
for i in {1..70}; do
    export i=${i}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk 'NR==ENVIRON["i"]{print $2}' | \
    while read prot; do
      export prot=${prot}
      echo ${trait} ${INF}/mr/GWAS/trait/${trait}-${prot}.gz > ${INF}/mr/GWAS/trait/gsmr_${trait}_${prot}
      gcta-1.9 --mbfile ${INF}/mr/gsmr/ref/gsmr_${prot} \
               --gsmr-file ${INF}/mr/gsmr/prot/gsmr_${prot} ${INF}/mr/GWAS/trait/gsmr_${trait}_${prot} \
               --gsmr-direction 0 \
               --clump-r2 0.1 --gwas-thresh 5e-8 --diff-freq 0.4 --heidi-thresh 0.05 --gsmr-snp-min 10 --effect-plot \
               --out ${INF}/mr/GWAS/out/${trait}-${prot}

      R --no-save -q <<\ \ \ \ \ \ END
        INF <- Sys.getenv("INF")
        trait <- Sys.getenv("trait")
        p <- Sys.getenv("prot")
        source(file.path(INF,"rsid","gsmr_plot.r"))
        gsmr_data <- read_gsmr_data(paste0(INF,"/mr/GWAS/out/",trait,"-",p,".eff_plot.gz"))
        gsmr_summary(gsmr_data)
        pdf(paste0(INF,"/mr/GWAS/out/",trait,"-",p,".eff_plot.pdf"))
        par(mar=c(6,6,5,1),mgp=c(4,1,0),xpd=TRUE)
        plot_gsmr_effect(gsmr_data, p, trait, colors()[75])
        dev.off()
      END
    done
done
# done previously
# echo ${INF}/mr/gsmr/ref/${prot} > ${INF}/mr/gsmr/ref/gsmr_${prot}
# echo ${prot} ${INF}/mr/gsmr/prot/${prot}.gz > ${INF}/mr/gsmr/prot/gsmr_${prot}

# The following code checks for validity of GSMR results (FDR<=0.05).

Rscript -e '
  options(width=200)
  INF <- Sys.getenv("INF")
  suppressMessages(library(dplyr))
  d <- file.path(INF,"mr","gsmr","trait")
  inf1 <- select(gap.datasets::inf1,prot,target.short)
  gsmr_efo <- read.delim(file.path(INF,"mr","gsmr","gsmr-efo.txt")) %>%
              filter(fdr<=0.01) %>%
              left_join(inf1,by=c("protein"="target.short")) %>%
              mutate(file_gwas=paste(prot,id,"rsid.txt",sep="-"),
                     bfile=file.path(INF,"INTERVAL","per_chr",paste0("interval.imputed.olink.chr_",chr)),
                     proxy=NA,p_proxy=NA,rsq=NA)
  gsmr_efo_reduce <- read.delim(file.path(INF,"mr","gsmr","gsmr-efo-reduce.txt")) %>%
                     filter(fdr<=0.01) %>%
                     left_join(inf1,by=c("protein"="target.short")) %>%
                     mutate(file_gwas=paste(prot,id,"rsid.txt",sep="-"),
                            bfile=file.path(INF,"INTERVAL","per_chr",paste0("interval.imputed.olink.chr_",chr)),
                            proxy=NA,p_proxy=NA,rsq=NA)
  qtl_lookup <- function(dat,panel="INTERVAL",pthreshold=1e-5,pop="EUR",plink_bin=NULL,xlsx=NULL)
  {
    for(i in 1:nrow(dat))
    {
       z <- slice(dat,i)
       pqtl <- z[["pqtl"]]
       cat(z[["file_gwas"]],"\n")
       gwas <- read.table(file.path(d,basename(z[["file_gwas"]])),header=TRUE) %>%
               arrange(p)
       h <-  filter(gwas,p<=pthreshold)
       panel_snps <- c(pqtl,pull(h,SNP))
       if (panel=="1000Genomes")
       {
          if (length(panel_snps)>500) stop("too many SNPs -- put a more stringent ptreshold")
          r <- ieugwasr::ld_matrix(panel_snps,pop="EUR")
       } else
       r <- ieugwasr::ld_matrix(panel_snps,with_alleles=TRUE,pop=pop,bfile=z[["bfile"]],plink_bin=plink_bin)
       cn <- colnames(r)
       inside <- pqtl==gsub("_[A-Z]*","",cn)
       nn <- c(cn[inside],cn[!inside])
       if (length(r)==1) r_mat <- r else r_mat <- r[nn,nn]
       r2_mat <- r_mat^2
       colnames(r2_mat) <- gsub("_[A-Z]*","",colnames(r2_mat))
       rownames(r2_mat) <- gsub("_[A-Z]*","",rownames(r2_mat))
       success <- 1
       snps <- intersect(pull(h,SNP),colnames(r2_mat))
       if (length(snps)<1) next else if (length(snps)==1) dat[i,c("proxy","p_proxy","rsq")] <- c(z[c("qtl","p_qtl")],1)
       else while(length(snps)>1)
       {
         success <- -1
         proxy <- snps[1]
         snps <- setdiff(snps,proxy)
         r2_i <- r2_mat[z[["pqtl"]],proxy]
         p_proxy <- filter(gwas,SNP==proxy) %>%
                    slice(which.min(p)) %>%
                    pull(p)
         cat("Same locus", i,z[["protein"]],z[["id"]],z[["Disease"]],z[["pqtl"]],z[["qtl"]],proxy,r2_i,z[["p_qtl"]],"\n",sep="\t")
         if(!is.null(r2_i)&!is.na(r2_i)) if(r2_i>0.8) {success <- 1; break}
       }
       if (success==-1) while (length(snps)>1)
       {
         success <- -1
         proxy <- snps[1]
         snps <- setdiff(snps,proxy)
         r2_i <- r2_mat[z[["pqtl"]],proxy]
         p_proxy <- filter(gwas,SNP==proxy) %>%
                    slice(which.min(p)) %>%
                    pull(p)
         cat("Independent locus",i,z[["protein"]],z[["id"]],z[["Disease"]],z[["pqtl"]],z[["qtl"]],proxy,r2_i,z[["p_qtl"]],"\n",sep="\t")
         if(!is.null(r2_i)&!is.na(r2_i)) if(r2_i<0.8) {success <- 1; break}
       }
       if (success==1)
       {
         dat[i,"proxy"] <- proxy
         dat[i,"p_proxy"] <- p_proxy
         dat[i,"rsq"] <- r2_i
       }
    }
    if (!is.null(xlsx))
    {
      suppressMessages(library(openxlsx))
      wb <- createWorkbook(xlsx)
      hs <- createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=12, fontName="Arial Narrow", fgFill="#4F80BD")
      proxies <- select(dat,protein,id,Disease,fdr,pqtl,p,qtl,p_qtl,proxy,p_proxy,rsq)
      for (sheet in "proxies")
      {
        addWorksheet(wb,sheet,zoom=150)
        writeData(wb,sheet,sheet,xy=c(1,1),headerStyle=createStyle(textDecoration="BOLD",
                 fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
        body <- get(sheet)
        writeDataTable(wb, sheet, body, xy=c(1,2), headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
        freezePane(wb, sheet, firstActiveCol=3, firstActiveRow=3)
        width_vec <- apply(body, 2, function(x) max(nchar(as.character(x))+2, na.rm=TRUE))
      # width_vec_header <- nchar(colnames(body))+2
        setColWidths(wb, sheet, cols = 1:ncol(body), widths = width_vec)
        writeData(wb, sheet, tail(body,1), xy=c(1, nrow(body)+2), colNames=FALSE, borders="rows", borderStyle="thick")
      }
      saveWorkbook(wb, file=xlsx, overwrite=TRUE)
    }
    dat
  }
  proxies <- qtl_lookup(gsmr_efo,plink_bin="/rds/user/jhz22/hpc-work/bin/plink",
                        xlsx=file.path(INF,"mr","gsmr","r2_INTERVAL.xlsx")) %>%
             select(protein,id,Disease,fdr,pqtl,p,qtl,p_qtl,proxy,p_proxy,rsq)
  write.table(proxies,file=file.path(INF,"mr","gsmr","r2_INTERVAL.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
  proxies <- qtl_lookup(gsmr_efo_reduce,plink_bin="/rds/user/jhz22/hpc-work/bin/plink",
                        xlsx=file.path(INF,"mr","gsmr","r2_INTERVAL_reduce.xlsx")) %>%
             select(protein,id,Disease,fdr,pqtl,p,qtl,p_qtl,proxy,p_proxy,rsq)
  write.table(proxies,file=file.path(INF,"mr","gsmr","r2_INTERVAL_reduce.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
# file_r2=paste(prot,id,Disease,pqtl,"r2.txt",sep="-")
# write.table(data.frame(d,id,file_r2,file_gwas),row.names=FALSE,col.names=FALSE,sep=",")
' | \
parallel --dry-run --csv "
  export d={1}
  export id={2}
  export r2=\${d}{3}
  export gwas=\${d}{4}
  export ids=~/INF/work/INTERVAL.rsid
  echo \${id} {3}
  awk '\$2>=0.8 {gsub(/_[A-Z]*/,\"\",\$1);print \$1,\$2}' \"\${r2}\" | sort -k2,2g
  awk '\$2>=0.8 {gsub(/_[A-Z]*/,\"\",\$1);print \$1,\$2}' \"\${r2}\" | \
  cut -d' ' -f1 | sed 's/\"//g' | grep -f - -w \"\${gwas}\" | sort -k7,7g
  awk '\$2>=0.8 {gsub(/_[A-Z]*/,\"\",\$1);print \$1,\$2}' \"\${r2}\" | \
  cut -d' ' -f1 | sed 's/\"//g' | grep -f - -w \"\${gwas}\" | cut -d' ' -f1 | grep -f - -w \"\${ids}\"
" | sed 's/^[ ]//' > script

# awk -vFS="\t" 'NR>1 && $8<0.05 {gsub(/-/,".",$1);print "trait/"$1"-"$2"-rsid.txt"}' gsmr-efo.txt | xargs -l -I{} awk '$7<1e-5' {}
# No LD information from OpenGWAS for the following protein-id:
# incl <- c("CXCL5-ebi-a-GCST004131","CD6-ebi-a-GCST004133","CD6-ebi-a-GCST003156","CD5-ieu-a-1112",
#           "CCL4-ieu-b-69","CCL4-ieu-a-996","CCL4-ebi-a-GCST004133","CCL4-ebi-a-GCST004131")
# x <- sapply(incl, function(x) grepl(x,z$f))
# tokeep <- as.logical(apply(x,1,sum))
# if (all) gsmr_efo <- z else gsmr_efo <- filter(z,tokeep)
#      proxy <- z[["qtl"]]
#      r2_i <- r2_mat[z[["pqtl"]],z[["qtl"]]]
