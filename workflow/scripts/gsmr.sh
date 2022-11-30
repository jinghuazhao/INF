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
      gsmr_efo <- left_join(gsmr,pQTLtools::inf1[c("prot","target.short")], by=c("Exposure"="prot")) %>%
                  left_join(filter(metal,prot %in% (gsmr$Exposure)), by=c("Exposure"="prot")) %>%
                  left_join(efo,by=c("Outcome"="id")) %>%
                  left_join(gwas,by=c("Exposure"="prot","Outcome"="id")) %>%
                  rename(protein=Exposure,id=Outcome,Disease=trait,Ncase=ncase,Ncontrol=ncontrol) %>%
                  mutate(protein=target.short,fdr=p.adjust(p,method="fdr")) %>%
                  select(protein,Disease,id,nsnp,fdr,Ncase,Ncontrol,Ntotal,bxy,se,pqtl,p,qtl,p_qtl,chr) %>%
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
                  left_join(pQTLtools::inf1[c("prot","target.short")], by=c("Exposure"="prot")) %>%
                  left_join(efo,by=c("Outcome"="MRBASEID")) %>%
                  rename(protein=Exposure,MRBASEID=Outcome) %>%
                  mutate(protein=target.short,fdr=p.adjust(p,method="fdr")) %>%
                  select(protein,MRBASEID,trait,bxy,se,p,nsnp,fdr,Ncases,Ncontrols,id,uri,Zhengetal) %>%
                  arrange(fdr)
      subset(gsmr_efo[setdiff(names(gsmr_efo),c("Zhengetal","uri"))],fdr<=0.05)
    }
  END
