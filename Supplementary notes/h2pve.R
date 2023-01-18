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
np <- read.table(paste(INF, "work", "INF1.merge.nosig", sep="/"), header=FALSE,
                 col.names = c("prot", "uniprot.id"))
annot$pQTL <- rep(NA, nrow(annot))
no.pQTL.ind <- which(annot$uniprot.id %in% np$uniprot)
annot$pQTL[no.pQTL.ind] <- "red"
annot$pQTL[-no.pQTL.ind] <- "blue"
annot <- annot[order(annot$pc.belowLOD.new, decreasing = T),]
annot <- annot[-grep("^BDNF$", annot$ID),]
annot <- annot %>% select(ID,MissDataProp) %>% rename(prot=ID)

INF1 <- pQTLdata::inf1[c("prot","target.short","alt_name")] %>%
        mutate(target.short=if_else(is.na(alt_name),target.short,alt_name))
hlp <- read.csv(file.path(INF,"ldak","h2-ldak-pve.csv")) %>%
       left_join(annot) %>%
       left_join(INF1) %>%
       arrange(h2_interval) %>%
       mutate(x=1:n())

interval <- hlp[c("target.short", "h2_interval", "SE_h2_interval", "MissDataProp", "x")]; names(interval)[2:3] <- c("h2","h2se")
interval <- data.frame(interval,source="(a) INTERVAL h2")
scallop <- hlp[c("target.short", "h2_scallop", "SE_h2_scallop", "MissDataProp", "x")]; names(scallop)[2:3] <- c("h2","h2se")
scallop <- data.frame(scallop,source="(b) SCALLOP")
pve <- hlp[c("target.short", "pve", "SE_pve", "MissDataProp", "x")]; names(pve)[2:3] <- c("h2","h2se")
pve <- data.frame(pve,source="(c) PVE")  %>% mutate(source=gsub("(c) PVE","PVE",source,fixed=TRUE))
isp <- rbind(interval,scallop,pve)

pve_order <- filter(pve,!is.na(h2)) %>%
             arrange(desc(h2)) %>%
             mutate(xtick=1:n())
attach(pve_order)
png(file.path(INF,"h2","SF-PVE.png"),width=15,height=8,units="in",pointsize=8,res=300)
    par(mar=c(10,5,1,1))
    plot(h2,cex=2,pch=19,xaxt="n",xlab="",ylab="",cex.axis=1.2)
    segments(xtick,h2-1.96*h2se,xtick,h2+1.96*h2se)
    axis(1, at=xtick, labels=target.short, lwd.tick=0.5, lwd=0, las=2, hadj=1, cex.axis=1.2)
    mtext("PVE",side=2,line=2.5,cex=1.5)
    mtext("Ordered protein",side=1,line=8.5,cex=1.5,font=1)
dev.off()
detach(pve_order)

p <- ggplot(pve_order,aes(xtick, h2))+
     theme_bw()+
     theme(
           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           axis.title.x = element_text(size=20),
           axis.title.y = element_text(size=20),
           axis.line = element_line(colour = "black"),
           axis.ticks.length=unit(.3, "cm"),
           axis.text.x = element_text(angle=90,color=if_else(pve_order$MissDataProp>=80, "red", "black"), size=16, hjust=1),
           axis.text.y = element_text(size=16)
          )+
     geom_point(size=4)+
     geom_errorbar(aes(ymin=h2-1.96*h2se,ymax=h2+1.96*h2se),position=position_dodge(0.05),width=0.3)+
     scale_y_continuous(breaks=seq(0,0.5,by=0.05),position="left")+
     scale_x_continuous(breaks=pve_order$xtick,labels=pve_order$target.short,position="bottom")+
     ylab("Proportion of variance explained (PVE)")+
     xlab("Ordered protein")
ggsave(p,filename=file.path(INF,"h2","SF-PVE-ggplot2.png"),device="png",dpi = 300, units="in", width=20, height=12)

p <- ggplot(isp,aes(y = x, x = h2))+
     theme_bw()+
     theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           axis.line = element_line(colour = "black"),
           axis.ticks.length=unit(.1, "cm"),
           axis.title = element_text(size=18),
           axis.text.x = element_text(size=16),
           axis.text.y = element_text(color=if_else(isp$MissDataProp>=80, "red", "black"), size=12),
           strip.text = element_text(face="bold",size=18)
          )+
     facet_wrap(~source,ncol=3,scales="free_x")+
     geom_errorbar(aes(xmin = h2-1.96*h2se, xmax = h2+1.96*h2se, colour=source), show.legend=FALSE)+
     geom_vline(lty=2, aes(xintercept=0), colour = "red")+
     scale_y_continuous(breaks=isp$x,labels=isp$target.short,position="left")+
     geom_point(size=2)+
     xlab("Heritability estimates using INTERVAL data and PVE according to pQTLs")+
     ylab("Protein")
ggsave(p,filename=file.path(INF,"h2","h2-pve-ggplot2.png"),device="png",dpi = 300, units="in", width=12, height=20)

pve <- hlp[c("target.short", "pve", "SE_pve", "MissDataProp", "x")]; names(pve)[2:3] <- c("h2","h2se")
pve <- data.frame(pve,source="(b) PVE")
isp2 <- rbind(interval,pve)

p <- ggplot(isp2,aes(y = x, x = h2))+
     theme_bw()+
     theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           axis.line = element_line(colour = "black"),
           axis.ticks.length=unit(.1, "cm"),
           axis.title = element_text(size=18),
           axis.text.x = element_text(size=16),
           axis.text.y = element_text(color=if_else(isp2$MissDataProp>=80, "red", "black"), size=12),
           strip.text = element_text(face="bold",size=18)
          )+
     facet_wrap(~source,ncol=2,scales="free_x")+
     geom_segment(aes(x = h2-1.96*h2se, xend = h2+1.96*h2se, yend = x, colour=source), show.legend=FALSE)+
     geom_vline(lty=2, aes(xintercept=0), colour = "red")+
     scale_y_continuous(breaks=isp2$x,labels=isp2$target.short,position="left")+
     geom_point(size=2)+
     xlab("Heritability estimates using INTERVAL data and PVE according to pQTLs")+
     ylab("Protein")
ggsave(p,filename=file.path(INF,"h2","h2-pve.png"),device="png",dpi = 300, units="in", width=12, height=20)

test <- function()
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
