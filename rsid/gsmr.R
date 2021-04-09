library(cowplot)
library(ggplot2)
library(gap)

INF <- Sys.getenv("INF")
# https://wilkelab.org/cowplot/articles/shared_legends.html

hgi <- read.table(file.path(INF,"HGI-ukb","INF1_A2-B2-C2.gsmr"), header=TRUE)
pdf(file.path(INF,"HGI","INF1_A2-B2-C2.pdf"))
mr <- gsmr(hgi, "LIF.R", c("A2","B2","C2"))
write.table(with(mr,r), file=file.path(INF,"HGI","INF1_A2-B2-C2.csv"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
ES_plots <- with(mr,plots)
legend <- get_legend(ES_plots[[1]])
plots <- align_plots(ES_plots[[1]],ES_plots[[2]],ES_plots[[3]], align = 'v', axis = 'l')
prow <- plot_grid(plots[[2]], plots[[3]], legend, labels = c("B2", "C2"), rel_widths = c(1, 1, .3), nrow = 1)
plot_grid(plots[[1]], prow, labels="A2", ncol=1)
dev.off()

dat <- read.table(file.path(INF,"gsmr","INF1_CAD-FEV1.gsmr"), header=TRUE)
pdf(file.path(INF,"gsmr","INF1_CAD-FEV1.pdf"))
mr <- gsmr(dat, "LIF.R", c("CAD","FEV1"))
write.table(with(mr,r), file=file.path(INF,"gsmr","INF1_CAD-FEV1.csv"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
ES_plots <- with(mr,plots)
library(cowplot)
prow <- plot_grid(
    ES_plots[[1]] + theme(legend.position="none"),
    ES_plots[[2]] + theme(legend.position="none"),
    align = 'vh',
    labels = c("CAD", "FEV1"),
    hjust = -1,
    nrow = 1
  )
legend <- get_legend(ES_plots[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(prow, legend, rel_widths = c(3, .4))
dev.off()
