fetch_region <- function()
{
  suppressMessages(library(dplyr))
  INF <- Sys.getenv("INF")
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
  f <- file.path(INF,"work","mono.tsv.gz")
  mono <- seqminer::tabix.read.table(tabixFile = f, tabixRange = r, stringsAsFactors = FALSE)
  names(mono) <- c("VARIANT","ID_dbSNP49","CHR","BP","REF","ALT","EFFECT_INT","SE_INT","MLOG10P_INT","ALT_FREQ_INT","INFO_INT")
  mono <- rename(mono, snpid=VARIANT,rsid=ID_dbSNP49,chr=CHR,pos=BP,a1=ALT,a2=REF,b=EFFECT_INT,se=SE_INT,log10p=MLOG10P_INT,af=ALT_FREQ_INT) %>%
          select(snpid, rsid, chr, pos, a2, a1, b, se, af, log10p) %>%
          mutate(snpid=if_else(a1<a2,paste0("chr",chr,":",pos,"_",a1,"_",a2),paste0("chr",chr,":",pos,"_",a2,"_",a1)),log10p=-log10p)
  wide <- right_join(d,mono)[grepl("chr|pos|snpid|rsid|b|se|log10p",names(d))] %>%
          mutate(is_lead_var=if_else(snpid=="chr1:159175354_A_G",2,1)) %>%
          arrange(chr,pos)
  for(i in 1:6)
  {
    t <- wide[grepl(paste("chr","pos",g[i],"is_lead_var",sep="|"),names(wide))]
    names(t) <- c("chr","pos","b","se","log10p","is_lead_var")
    t <- mutate(t,track=g[i])
    assign(g[i],t)
    if(i==1) d <- t
    else d <- bind_rows(d,t)
  }
  mono <- mutate(mono,track="Monocytes",is_lead_var=if_else(snpid=="chr1:159175354_A_G",2,1)) %>%
          select(-snpid,-rsid)
  long <- bind_rows(d,mono[c("chr","pos","a1","a2","af","b","se","log10p","is_lead_var","track")])
  list(wide=wide,long=long)
}

plot_region <- function(dat,g)
{
  suppressMessages(library(ggplot2))
  p <- ggplot(dat, aes(pos,-log10p)) +
       geom_blank() +
       geom_point(aes(colour = factor(is_lead_var), alpha = 0.7)) +
       geom_point(data = dat %>% filter(is_lead_var %in% c(1,2)), aes(colour = factor(is_lead_var))) +
       theme_light() +
       ylab(expression(paste("-", log[10], " p-value"))) +
       scale_x_continuous(limits = region_coords, expand = c(0, 0)) +
       facet_grid(track ~ ., scales = "free_y", space = "free") +
       theme(
               plot.margin = unit(c(0.1, 1, 0.1, 1), "line"),
               axis.text.x = element_blank(),
               legend.position = "none",
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.text.y = element_text(colour = "black"),
               strip.background = element_rect(fill = "white")
            ) +
       scale_color_manual(values=c("#377eb8","#e41a1c","#984ea3","#4daf4a","#ff7f00","#a65628","#fed976")) +
       xlab("Position")
}

g <- c("CCL2", "CCL7", "CCL8", "CCL11", "CCL13", "CXCL6")
p <- c("MCP.1", "MCP.3", "MCP.2", "CCL11", "MCP.4", "CXCL6")
dat <- fetch_region()
wide <- with(dat,wide)
long <- with(dat,long) %>%
        mutate(lcl=b-1.96*se,ucl=b+1.96*se,track=factor(track,levels=c(g,"Monocytes")))
region_coords <- c(first(wide$pos),last(wide$pos))
dat_text <- data.frame(label = c(paste0("\nred: ", "rs12075", "\nblue: ", "others")))
l <- plot_region(long) +
     geom_text(data    = dat_text,
               mapping = aes(x = -Inf, y = -Inf, label = label),
               hjust   = -0.1,
               vjust = -3.4
              )
ggsave("rs12075.png",l,height=20,width=10)

rs12075 <- filter(long,pos==159175354)
m <- with(rs12075,weighted.mean(b,1/se^2))
format_p <- function(p) paste("p =", substring(prettyNum(p, digits=2, scientific=TRUE), 2))
f <- ggplot(data = rs12075,
     aes(x = 7:1, y = b, label=paste0(track,": ",format_p(10^log10p)))) +
     geom_pointrange(aes(ymin=lcl, ymax=ucl), size=0.7) +
     geom_text(y=-15, hjust=0) +
     geom_abline(intercept = 0, slope = 0) +
     geom_abline(intercept = m, slope = 0, lty = 2) +
     ylim(c(-15,5)) +
     coord_flip() +
     theme_bw() + theme(
                         axis.text.y = element_blank(),
                         panel.grid = element_blank(),
                         panel.border = element_blank(),
                         text = element_text(size=17)
                       ) +
     xlab("Phenotypes") +
     ylab("Effect size")

# arrangement of plots
require(cowplot)
fl <- plot_grid(f, l, nrow = 1, labels = "AUTO", label_size = 12, align = "h")
ggsave("rs12075-forest-assoc.png",height=20,width=10)

# gunzip -c ~/rds/results/public/gwas/blood_cell_traits/astle_2016/raw_results/blood_cell_traits/gzipped_interval/mono.tsv.gz | \
# bgzip -f > ${INF}/work/mono.tsv.gz
# tabix -S1 -s3 -b4 -e4 ${INF}/work/mono.tsv.gz
