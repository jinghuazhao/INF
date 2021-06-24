options(width=200)

library(pQTLtools)
library(Biobase)
library(dplyr)
library(ggplot2)

HOME <- Sys.getenv("HOME")
dir <- paste0(HOME,"/rds/post_qc_data/interval/phenotype/olink_proteomics/post-qc/")
rds <- readRDS(paste0(dir,"eset.inf1.flag.out.outlier.in.rds"))
x <- get.prop.below.LLOD(rds)
annot <- fData(x)
annot$MissDataProp <- as.numeric( gsub("\\%$", "", annot$MissDataProp) )

INF <- Sys.getenv("INF")
np <- read.table(paste(INF, "work", "INF1.merge.nosig", sep="/"), head=F,
                 col.names = c("prot", "uniprot.id"))
annot$pQTL <- rep(NA, nrow(annot))
no.pQTL.ind <- which(annot$uniprot.id %in% np$uniprot)
annot$pQTL[no.pQTL.ind] <- "red"
annot$pQTL[-no.pQTL.ind] <- "blue"
annot <- annot[order(annot$pc.belowLOD.new, decreasing = T),]
annot <- annot[-grep("^BDNF$", annot$ID),]
annot <- annot %>% select(ID,MissDataProp) %>% rename(prot=ID)

hlp <- read.csv(file.path(INF,"ldak","h2-ldak-pve.csv")) %>%
       left_join(annot) %>%
       left_join(pQTLtools::inf1[c("prot","target.short")]) %>%
       arrange(h2_interval) %>%
       mutate(x=1:n())

interval <- hlp[c("target.short", "h2_interval", "SE_h2_interval", "MissDataProp", "x")]; names(interval)[2:3] <- c("h2","h2se")
interval <- data.frame(interval,source="(a) INTERVAL")
scallop <- hlp[c("target.short", "h2_scallop", "SE_h2_scallop", "MissDataProp", "x")]; names(scallop)[2:3] <- c("h2","h2se")
scallop <- data.frame(scallop,source="(b) SCALLOP")
pve <- hlp[c("target.short", "pve", "SE_pve", "MissDataProp", "x")]; names(pve)[2:3] <- c("h2","h2se")
pve <- data.frame(pve,source="(c) PVE")
isp <- rbind(interval,scallop,pve)

p <- ggplot(isp,aes(y = x, x = h2))+
     theme_bw()+
     theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           axis.line = element_line(colour = "black"), axis.ticks.length=unit(.1, "cm"),
           axis.text.y = element_text(color=if_else(isp$MissDataProp>=80, "red", "black"), size=7))+
     facet_wrap(~source,ncol=3,scales="free_x")+
     geom_segment(aes(x = h2-1.96*h2se, xend = h2+1.96*h2se, yend = x, colour=source), show.legend=FALSE)+
     geom_vline(lty=2, aes(xintercept=0), colour = "red")+
     scale_y_continuous(breaks=isp$x,labels=isp$target.short,position="left")+
     geom_point()+
     xlab("Heritability/PVE")+
     ylab("")
ggsave(p,filename=file.path(INF,"h2","h2-pve-ggplot2.png"),device="png",dpi = 300, units="in", width=12, height=20)

function test()
{
p <- ggplot(isp,aes(x=x, y=h2))+
     theme_bw()+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           axis.line = element_line(colour = "black"), axis.ticks.length=unit(.2, "cm"), axis.text.x = element_text(angle = 90))+
     facet_wrap(~source,ncol=1,scales="free_y")+
     geom_segment(aes(y = h2-1.96*h2se, yend = h2+1.96*h2se, xend=x, colour=source), show.legend=FALSE)+
     geom_point()+
     scale_x_continuous(breaks=isp$x,labels=isp$target.short,position="bottom")+
     xlab("")+
     ylab("Heritability/PVE")

ggplot(hlp,aes(x=x,y=h2_interval)) + geom_point()
# ggplot(hlp,aes(x=x,y=h2_interval)) + geom_point() + geom_point(y=hlp$h2_scallop)+geom_point(y=hlp$pve)
# geom_hline(lty=2, aes(yintercept=0), colour = "red")+
}
