#!/bin/bash

(
  seq 91 | \
  parallel -C ' ' -j8 '
    export protein___uniprot=$(cut -d " " -f1-28 --complement ${INF}/h2/s.sample | head -1 | tr " " "\n" | grep -v BDNF | awk -vi={} "NR==i")
    echo ${protein___uniprot}
    Rscript -e "
      options(width=200)
      suppressMessages(library(coda))
      INF <- Sys.getenv(\"INF\")
      protein___uniprot <- Sys.getenv(\"protein___uniprot\")
      eps <- 0.02
      U <- as.mcmc(scan(file.path(INF,\"h2\",paste0(protein___uniprot,\".BLR_varU.dat\")))[-(1:10000)])
      E <- as.mcmc(scan(file.path(INF,\"h2\",paste0(protein___uniprot,\".BLR_varE.dat\")))[-(1:10000)])
      e <- as.mcmc(cbind(U, E, h2 = U / ((1 + eps) * U + E)))
      summary(e)$statistics
      HPDinterval(e)
      sumstats <- data.frame(id=c(\"U\",\"E\",\"h2\"),as.matrix(summary(e)[[\"statistics\"]]))
      write.table(sumstats,file.path(INF,\"h2\",paste0(protein___uniprot,\".BLR.txt\")),row.names=FALSE,quote=FALSE)
      png(file.path(INF,\"h2\",paste0(protein___uniprot,\".BLR.png\")),res=300,height=10,width=8,units=\"in\")
      plot(e)
      dev.off()
    "
  '
) > ${INF}/h2/BLR.log

(
  echo prot uniprot h2_BLR se_BLR h2_GCTA se_GCTA
  join <(grep h2 ${INF}/h2/*___*txt | sed 's|'"${INF}/h2/"'||;s/___/ /;s/.BLR.txt:h2//' | cut -d' ' -f1-4 | sort -k1,1) \
       <(sed '1d' ${INF}/h2/h2.dat | sort -k1,1) \
) > ${INF}/h2/BLR_GCTA.dat

Rscript -e '
  library(dplyr)
  library(data.table)
  INF <- Sys.getenv("INF")
  BLR_GCTA <- read.table(file.path(INF,"h2","BLR_GCTA.dat"),as.is=TRUE,header=TRUE) %>% arrange(desc(h2_GCTA)) %>% mutate(x=1:n())
  h2 <- data.table::melt(BLR_GCTA[c("prot","uniprot","h2_GCTA","h2_BLR")],id.vars=c("prot","uniprot"),variable.name=c("source")) %>%
        rename(h2=value) %>%
        mutate(source=gsub("h2_","",source))
  se <- data.table::melt(BLR_GCTA[c("prot","uniprot","se_GCTA","se_BLR")],id.vars=c("prot","uniprot"),variable.name=c("source")) %>%
        rename(se=value) %>%
        mutate(source=gsub("se_","",source))
  dat <- left_join(h2,se) %>%
         left_join(select(BLR_GCTA,prot,x))
  library(ggplot2)
  p <- ggplot(dat,aes(y = x, x = h2))+
       theme_bw()+
       theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             axis.line = element_line(colour = "black"),
             axis.ticks.length=unit(.1, "cm"),
             axis.title = element_text(size=16),
             axis.text.x = element_text(size=14),
             axis.text.y = element_text(size=12),
             strip.text = element_text(face="bold",size=18)
            )+
       facet_wrap(~source,ncol=2,scales="free_x")+
       geom_errorbar(aes(xmin = h2-1.96*se, xmax = h2+1.96*se, colour=source), show.legend=FALSE)+
       geom_vline(lty=2, aes(xintercept=0), colour = "blue")+
       geom_vline(lty=2, aes(xintercept=1), colour = "red")+
       scale_y_continuous(breaks=dat$x,labels=dat$prot,position="left")+
       geom_point(size=2)+
       xlab("Heritability estimates based on GCTA and BLR")+
       ylab("Ordered protein")
  ggplot2::ggsave(file.path(INF,"h2","h2_GCTA_BLR.png"),height=20,width=10)
'
