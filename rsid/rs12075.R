fetch_region <- function(src="chen")
{
  chr <- 1
  start <- 159175353
  end <- 159525679
  M <- 1e6
  r <- paste0(chr,":",start-M,"-",end+M)
  for (i in 1:6)
  {
    f <- file.path(INF,"METAL","gwas2vcf","LP",paste0(p[i],".tsv.gz"))
    t <- seqminer::tabix.read.table(tabixFile = f, tabixRange = r, stringsAsFactors = FALSE)
    names(t) <- c("chr","pos","snpid","a1","a2","af","b","se","log10p","n")
    t <- rename_with(t, ~ paste0(paste0(g[i],"."),.x),.cols=4:10)
    if(i==1) d <- t
    else d <- left_join(d,t)
    assign(g[i],t)
    cat(p[i],g[i],"\n")
    print(subset(t,snpid=="chr1:159175354_A_G"))
    print(dim(t))
  }
  if (src!="chen")
  {
    f <- file.path(INF,"work","mono.tsv.gz")
    mono <- seqminer::tabix.read.table(tabixFile = f, tabixRange = r, stringsAsFactors = FALSE)
    names(mono) <- c("VARIANT","ID_dbSNP49","CHR","BP","REF","ALT","EFFECT_INT","SE_INT","MLOG10P_INT","ALT_FREQ_INT","INFO_INT")
    mono <- rename(mono, snpid=VARIANT,rsid=ID_dbSNP49,chr=CHR,pos=BP,a1=ALT,a2=REF,b=EFFECT_INT,se=SE_INT,log10p=MLOG10P_INT,af=ALT_FREQ_INT) %>%
            select(snpid, rsid, chr, pos, a2, a1, b, se, af, log10p) %>%
            mutate(snpid=if_else(a1<a2,paste0("chr",chr,":",pos,"_",a1,"_",a2),paste0("chr",chr,":",pos,"_",a2,"_",a1)),log10p=-log10p)
  } else {
    f <- file.path(INF,"work","mono-chen.tsv.gz")
    mono <- seqminer::tabix.read.table(tabixFile = f, tabixRange = r, stringsAsFactors = FALSE) %>%
    setNames(c("snpid","rsid","chr","pos","a1","a2","af","b","se","p")) %>%
    mutate(log10p=log10(p))
  }
  wide <- right_join(d,mono)[grepl("chr|pos|snpid|rsid|a1|a2|af|b|se|log10p",names(d))] %>%
          arrange(chr,pos)
  for(i in 1:6)
  {
    t <- wide[grepl(paste("chr","pos",g[i],sep="|"),names(wide))]
    names(t) <- c("chr","pos","a1","a2","af","b","se","log10p")
    t <- mutate(t,track=g[i])
    assign(g[i],t)
    if(i==1) d <- t
    else d <- bind_rows(d,t)
  }
  mono <- mutate(mono,track="Monocytes") %>%
          select(-snpid,-rsid)
  long <- bind_rows(d,mono[c("chr","pos","a1","a2","af","b","se","log10p","track")])
  save(wide,long,file=file.path(INF,"hotspots","rs12075.rda"),compress="xz")
  list(wide=wide,long=long)
}

plot_region <- function(dat)
{
  suppressMessages(library(ggplot2))
  p <- ggplot(dat, aes(pos,-log10p)) +
       geom_blank() +
       geom_point(aes(colour = factor(pos==159.175354), alpha = 0.7)) +
       geom_point(data = dat, aes(colour = factor(pos==159.175354))) +
       theme_light() +
       scale_x_continuous(limits = region_coords, expand = c(0, 0)) +
       facet_grid(track ~ ., scales = "free_y", space = "free") +
       theme(
               plot.margin = unit(c(0.1, 1, 0.1, 1), "line"),
               axis.title = element_text(size=18),
               axis.text.x = element_text(size=16),
               axis.text.y = element_text(size=12),
               legend.position = "none",
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.text.y = element_text(colour = "black"),
               strip.background = element_rect(fill = "white")
            ) +
       scale_color_manual(values=c("#377eb8","#e41a1c","#984ea3","#4daf4a","#ff7f00","#a65628","#fed976")) +
       xlab("Position (Mb)") +
       ylab(expression(paste("-", log[10], " p-value")))
}

suppressMessages(library(dplyr))

INF <- Sys.getenv("INF")
g <- c("CCL2", "CCL7", "CCL8", "CCL11", "CCL13", "CXCL6")
p <- c("MCP.1", "MCP.3", "MCP.2", "CCL11", "MCP.4", "CXCL6")
dat <- fetch_region()
wide <- with(dat,wide) %>% mutate(pos=pos/1e6)
long <- with(dat,long) %>%
        mutate(lcl=b-1.96*se,ucl=b+1.96*se,track=factor(track,levels=c(g,"Monocytes")),pos=pos/1e6)
region_coords <- c(first(wide$pos),last(wide$pos))
l <- plot_region(long)
ggsave(file.path(INF,"hotspots","rs12075.png"),l,height=20,width=10)

rs12075 <- filter(long,pos==159.175354)
f <- ggplot(data=rs12075, aes(y=1:7, x=b))+
     geom_point(size=2)+
     geom_errorbarh(aes(xmax = ucl, xmin = lcl, height=0.001))+
     scale_x_continuous(name="Effect size")+
     scale_y_continuous(breaks=1:7,label=with(rs12075,track),name="",trans="reverse")+
     geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)+
     theme_minimal()+
     theme(text=element_text(size=16, color="black"),
           axis.title = element_text(size=18),
           axis.text.x = element_text(size=16),
           panel.grid=element_blank(),
           panel.spacing = unit(1, "lines"))
ggsave(file.path(INF,"hotspots","rs12075-forest.png"),f,height=10,width=10)
  
# arrangement of plots
require(cowplot)
fl <- plot_grid(f, l, nrow = 1, labels = "AUTO", label_size = 16)
ggsave(file.path(INF,"hotspots","rs12075-forest-assoc.png"),height=20,width=10)

circos_plot <- function()
# circos plot
{
   t.genes <- filter(pQTLdata::hg19, pQTLdata::hg19$SYMBOL %in% g) %>%
              distinct() %>%
              select(chr,start,end,SYMBOL) %>%
              rename(t.chr=chr,t.start=start,t.end=end,t.gene=SYMBOL) %>%
              mutate(t.gene=as.character(unlist(t.gene)))
   suggestive <- transmute(long, chr=paste0("chr",chr),start=pos-1,end=pos,p=10^log10p,from="ACKR1",to=as.character(track)) %>%
                 filter(p<=1e-5 & to!="Monocytes") %>%
                 left_join(t.genes,by=c('to'='t.gene'))
   labels <- head(suggestive,1) %>%
             transmute(t.chr=chr,t.start=start,t.end=end,t.gene=from) %>%
             bind_rows(t.genes)
   suppressMessages(library(circlize))
   png(file.path(INF,"hotspots","hotspot-rs12075.png"),res=300,width=10,height=10,units="in")
   circos.clear()
   circos.par(start.degree = 90, track.height = 0.1, cell.padding = c(0, 0, 0, 0))
   circos.initializeWithIdeogram(species="hg19", track.height = 0.05, ideogram.height = 0.06, labels.cex=2, chromosome.index=paste0("chr",1:22))
   circos.genomicLabels(labels,labels.column=4, cex=1.2, font=4, side="inside")
   circos.genomicLink(suggestive[c("chr","start","end")], suggestive[c("t.chr","t.start","t.end")], col="blue", directional=1, border = 10, lwd = 2)
   dev.off()
}

circos_plot()

run_gassoc <- TRUE

if (run_gassoc)
{
  wide_z <- with(dat,wide) %>%
            mutate(
            CCL2=CCL2.b/CCL2.se,
            CCL7=CCL7.b/CCL7.se,
            CCL8=CCL8.b/CCL8.se,
            CCL11=CCL11.b/CCL11.se,
            CCL13=CCL13.b/CCL13.se,
            CXCL6=CXCL6.b/CXCL6.se,
            Monocytes=b/se,
            s=CCL2+CCL7+CCL8+CCL11+CCL13+CXCL6+Monocytes) %>%
            filter(!is.na(s)) %>%
            select(rsid,chr,pos,CCL2,CCL7,CCL8,CCL11,CCL13,CXCL6,Monocytes) %>%
            rename(marker=rsid)
  plink_bin <- "/rds/user/jhz22/hpc-work/bin/plink"
  chr <- 1
  bfile <- file.path(INF,"INTERVAL","per_chr",paste0("interval.imputed.olink.chr_",chr))
  r <- ieugwasr::ld_matrix(select(wide_z,marker),with_alleles=TRUE,pop="EUR",bfile=bfile,plink_bin=plink_bin)
  rnames <- gsub("_[A-Z]*","",colnames(r))
  wide_z <- subset(wide_z,marker %in% rnames)
  rsids <- intersect(rnames,with(wide_z,marker))
  ld <- r
  colnames(ld) <- rownames(ld) <- rnames
  ld <- ld[rsids,rsids]
  d <- subset(wide_z,marker %in% rsids)
  z <- d[c("CCL2","CCL7","CCL8","CCL11","CCL13","CXCL6")]
  rownames(z) <- with(d,marker)
  library(gassocplot2)
  pdf(file.path(INF,"hotspots","SF-rs12075-gassoc.pdf"),height=20,width=8)
  sap <- stack_assoc_plot(d[c("marker","chr","pos")], z, ld, traits=names(z), ylab="-log10(P)", legend=TRUE)
  grid::grid.draw(sap)
  dev.off()
}

chen_vars <- c("snpid","rsid","chr","pos","a1","a2","af","b","se","p")
wbc <- read.delim(file.path(INF,"work","wbc-chen.tsv.gz")) %>% setNames(chen_vars)
mono <- read.delim(file.path(INF,"work","mono-chen.tsv.gz")) %>% setNames(chen_vars)
baso <- read.delim(file.path(INF,"work","baso-chen.tsv.gz")) %>% setNames(chen_vars)

if (run_gassoc)
{
  wbc_z <- wbc %>% mutate(wbc=b/se) %>%
           select(snpid,rsid,wbc,chr,pos) %>% rename(marker=rsid)
  mono_z <- mono %>% mutate(mono=b/se) %>% select(snpid,mono)
  baso_z <- baso %>% mutate(baso=b/se) %>% select(snpid,baso)

  blood_traits <- wbc_z %>% left_join(mono_z) %>% left_join(baso_z) %>% filter(marker!="rs6413465") %>% select(chr,pos,marker,mono,baso,wbc)
  r <- ieugwasr::ld_matrix(select(blood_traits,marker),with_alleles=TRUE,pop="EUR",bfile=bfile,plink_bin=plink_bin)
  rnames <- gsub("_[A-Z]*","",colnames(r))
  blood_traits <- subset(blood_traits,marker %in% rnames)
  rsids <- intersect(rnames,with(blood_traits,marker))
  ld <- r
  colnames(ld) <- rownames(ld) <- rnames
  ld <- ld[rsids,rsids]
  d <- subset(blood_traits,marker %in% rsids)
  z <- d[c("mono","baso","wbc")] %>% setNames(c("Monocyte count","Basophil count","WBC"))
  rownames(z) <- with(d,marker)
  library(gassocplot2)
  pdf(file.path(INF,"hotspots","SF-rs12075-traits-gassoc.pdf"),height=20,width=8)
  sap <- stack_assoc_plot(d[c("marker","chr","pos")], z, ld, traits=names(z), ylab="-log10(P)", top.marker="rs12075",legend=TRUE)
  grid::grid.draw(sap)
  dev.off()
}

# stack_assoc_plot_save(sap, paste0("rs12075-gassoc.png"), 7, width=8, dpi=300)
# ggsave("rs12075-gassoc.png", plot=grid::grid.draw(sap), dpi=300, height=20, width=8, units="in", limitsize=FALSE)

wbc_rs12075 <- filter(wbc,pos==159175354) %>% mutate(track="WBC") %>% select(track, b,se)
mono_rs12075 <- filter(mono,pos==159175354) %>% mutate(track="Monocyte count") %>% select(track, b,se)
baso_rs12075 <- filter(baso,pos==159175354) %>% mutate(track="Basophil count") %>% select(track, b,se)
d <- bind_rows(filter(rs12075[c("track","b","se")],track!="Monocytes"), wbc_rs12075, mono_rs12075, baso_rs12075)
library(gap)
png(file.path(INF,"hotspots","SF-rs12075-forest.png"),height=7,width=7.5,units="in",res=300)
mr_forestplot(d, colgap.forest.left="0.05cm", fontsize=14, digits=3,
              leftlabs=c("Outcome","b","SE"),
              rightcols=c("ci","pval"), rightlabs=c("95%CI","P"),digits.pval=2,scientific.pval=TRUE,
              common=FALSE, random=FALSE, print.I2=FALSE, print.pval.Q=FALSE, print.tau2=FALSE,
              spacing=1.6,digits.TE=3,digits.se=3,xlab="Effect size")
dev.off()
