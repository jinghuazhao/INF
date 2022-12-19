fetch_region <- function()
{
  suppressMessages(library(dplyr))
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
  wide <- d[grepl("chr|pos|snpid|rsid|a1|a2|af|b|se|log10p",names(d))] %>%
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
  list(wide=wide,long=d)
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
g <- c("CCL2", "CCL7", "CCL8", "CCL11", "CCL13", "CXCL6")
p <- c("MCP.1", "MCP.3", "MCP.2", "CCL11", "MCP.4", "CXCL6")
dat <- fetch_region()
wide <- with(dat,wide) %>% mutate(pos=pos/1e6)
long <- with(dat,long) %>%
        mutate(lcl=b-1.96*se,ucl=b+1.96*se,track=factor(track,levels=g),pos=pos/1e6)
region_coords <- c(first(wide$pos),last(wide$pos))
l <- plot_region(long)

circos_plot <- function()
# circos plot
{
   t.genes <- filter(pQTLtools::hg19, pQTLtools::hg19$SYMBOL %in% g) %>%
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

# gunzip -c ~/rds/results/public/gwas/blood_cell_traits/astle_2016/raw_results/blood_cell_traits/gzipped_interval/mono.tsv.gz | \
# bgzip -f > ${INF}/work/mono.tsv.gz
# tabix -S1 -s3 -b4 -e4 ${INF}/work/mono.tsv.gz
