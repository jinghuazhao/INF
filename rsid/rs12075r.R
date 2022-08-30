blood_cell_trait <- function(d,trait,r,track="")
{
# f <- file.path(INF,"work",paste0(trait,".tsv.gz"))
# trait <- seqminer::tabix.read.table(tabixFile = f, tabixRange = r, stringsAsFactors = FALSE)
# names(trait) <- c("VARIANT","ID_dbSNP49","CHR","BP","REF","ALT","EFFECT_INT","SE_INT","MLOG10P_INT","ALT_FREQ_INT","INFO_INT")
# trait <- rename(trait, snpid=VARIANT,rsid=ID_dbSNP49,chr=CHR,pos=BP,
#                 a1=ALT,a2=REF,b=EFFECT_INT,se=SE_INT,log10p=MLOG10P_INT,af=ALT_FREQ_INT) %>%
#         select(snpid, rsid, chr, pos, a2, a1, b, se, af, log10p) %>%
#         mutate(snpid=if_else(a1<a2,paste0("chr",chr,":",pos,"_",a1,"_",a2),paste0("chr",chr,":",pos,"_",a2,"_",a1)),log10p=-log10p)
  f <- file.path("~/rds/results/public/gwas/blood_cell_traits/astle_2016/raw_results/blood_cell_traits/bgzipped_narrow_form",paste0(trait,".tsv.gz"))
  sumstats <- seqminer::tabix.read.table(tabixFile = f, tabixRange = r, stringsAsFactors = FALSE) %>%
              setNames(c("snpid","rsid","chr","pos","a2","a1","ALT_MINOR","DIRECTION","b","se","P","log10p","af","MA_FREQ")) %>%
              mutate(snpid=if_else(a1<a2,paste0("chr",chr,":",pos,"_",a1,"_",a2),paste0("chr",chr,":",pos,"_",a2,"_",a1)),log10p=-log10p) %>%
              select(snpid,rsid,chr,pos,a1,a2,b,se,af,log10p)
  wide <- sumstats %>%
          setNames(c("snpid","rsid","chr","pos",paste0(trait,".a1"),paste0(trait,".a2"),paste0(trait,".b"),paste0(trait,".se"),
                   paste0(trait,".af"),paste0(trait,".log10p"))) %>%
          arrange(chr,pos)
  sumstats <- mutate(sumstats,track=track) %>%
              select(-snpid,-rsid)
  list(long=sumstats,wide=wide)
}

fetch_region <- function()
{
  suppressMessages(library(dplyr))
  chr <- 1
  start <- 159175353
  end <- 159525679
  M <- 1e6
  r <- paste0(chr,":",start-M,"-",end+M)
  for (i in 1:length(g))
  {
    f <- file.path(INF,"METAL","gwas2vcf","LP",paste0(p[i],".tsv.gz"))
    t <- seqminer::tabix.read.table(tabixFile = f, tabixRange = r, stringsAsFactors = FALSE) %>%
         setNames(c("chr","pos","snpid","a1","a2","af","b","se","log10p","n")) %>%
         mutate(track=g[i])
#   t <- rename_with(t, ~ paste0(paste0(g[i],"."),.x),.cols=4:10)
    if(i==1) d <- t
    else d <- left_join(d,t)
    assign(g[i],t)
    cat(p[i],g[i],"\n")
    print(subset(t,snpid=="chr1:159175354_A_G"))
    print(dim(t))
  }
  baso <- blood_cell_trait(d,"baso",r,track="Basophils")
  mono <- blood_cell_trait(d,"mono",r,track="Monocytes")
  wbc <- blood_cell_trait(d,"wbc",r,track="WBC")
  long <- bind_rows(d,baso$long,mono$long,wbc$long)
  wide <- left_join(d,baso$wide) %>% left_join(mono$wide) %>% left_join(wbc$wide)
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

INF <- Sys.getenv("INF")
g <- "CCL13"
p <- "MCP.4"
dat <- fetch_region()
wide <- with(dat,wide) %>% mutate(pos=pos/1e6)
long <- with(dat,long) %>%
        mutate(lcl=b-1.96*se,ucl=b+1.96*se,track=factor(track,levels=c(g,"Basophils","Monocytes", "WBC")),pos=pos/1e6)
region_coords <- c(first(wide$pos),last(wide$pos))
l <- plot_region(long)
ggsave("rs12075.png",l,height=20,width=10)

rs12075 <- filter(long,pos==159.175354)
f <- ggplot(data=rs12075, aes(y=1:4, x=b))+
     geom_point(size=2)+
     geom_errorbarh(aes(xmax = ucl, xmin = lcl, height=0.001))+
     scale_x_continuous(name="Effect size")+
     scale_y_continuous(breaks=1:4,label=with(rs12075,track),name="",trans="reverse")+
     geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)+
     theme_minimal()+
     theme(text=element_text(size=16, color="black"),
           axis.title = element_text(size=18),
           axis.text.x = element_text(size=16),
           panel.grid=element_blank(),
           panel.spacing = unit(1, "lines"))
  
require(cowplot)
fl <- plot_grid(f, l, nrow = 1, labels = "AUTO", label_size = 16)
ggsave("rs12075-forest-assoc.png",height=20,width=10)
