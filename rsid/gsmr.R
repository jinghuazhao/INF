# Adapted from code by Felix Day 16/9/2015

weighted.median <- function(x, w)
{
   x.order <- x[order(x)]
   w.order <- w[order(x)]
   ws <- cumsum(w.order)-0.5*w.order
   ws <- ws / sum(w.order)
   below <- max(which(ws < 0.5))
   x.order[below] + (x.order[below+1]-x.order[below]) * (0.5-ws[below]) / (ws[below+1]-ws[below])
}

mr.boot = function(bXG, sebXG, bYG, sebYG, w, n.boot=1000, method="median")
{
   boot <- vector('numeric')
   for(i in 1:n.boot)
   {
       bXG.boot <- rnorm(length(bXG), mean = bXG, sd = sebXG)
       bYG.boot <- rnorm(length(bYG), mean = bYG, sd = sebYG)
       if(method=="Egger")
       {
         bYG.boot = bYG.boot*sign(bXG.boot)
         bXG.boot = abs(bXG.boot)
         boot[i] = coef(summary(lm(bYG.boot~bXG.boot,weights=sebYG^-2)))[2,1]
       } else {
         bIV.boot <- bYG.boot / bXG.boot
         boot[i] <- weighted.median(bIV.boot, w)
      }
   }
   sd(boot)
}

gsmr <- function(data, X, Y, alpha=0.05, other_plots=FALSE)
{
   c <- abs(qnorm(alpha/2))
   nxy <- vector("list", length(X)*length(Y))
   k <- 1
   ES_plot <- vector('list')
   for(i in Y)
   {
       dat <- subset(data.frame(SNP = data[["SNP"]],
                                bzx = eval(parse(text = paste0("data$b.", X))), SEzx = eval(parse(text = paste0("data$SE.", X))),
                                bzy = eval(parse(text = paste0("data$b.", i))), SEzy = eval(parse(text = paste0("data$SE.", i)))),
                    !is.na(bzx+SEzx+bzy+SEzy))
       SNP <- dat[["SNP"]]
       bzx <- dat[["bzx"]]; SEzx <- dat[["SEzx"]]
       bzy <- dat[["bzy"]]; SEzy <- dat[["SEzy"]]
       nxy[k] <- paste0(X, ".", i)
       m1 <- lm(bzy~bzx-1, weights=SEzy^-2); summary.m1 <- summary(m1); coef.summary.m1 <- coef(summary.m1)
       m <- lm(bzy~bzx, weights=SEzy^-2); summary.m <- summary(m); coef.summary.m <- coef(summary.m)
       bIVW <- coef.summary.m1[1,1]
       sebIVW <- coef.summary.m1[1,2] / min(with(summary.m1, sigma), 1)
       bEGGER <- coef.summary.m[2,1]
       sebEGGER <- coef.summary.m[2,2] / min(with(summary.m, sigma), 1)
       intEGGER <- coef.summary.m[1,1]
       seintEGGER <- coef.summary.m[1,2]
       bIV <- bzy / bzx
       weights <- (SEzy / bzx)^-2
       bIVW <- sum(bzy*bzx*SEzy^-2) / sum(bzx^2 * SEzy^-2)
       bWM <- weighted.median(bIV, weights)
       sebWM <- mr.boot(bzx, SEzx, bzy, SEzy, weights)
       penalty <- pchisq(weights*(bIV-bIVW)^2, df=1, lower.tail=FALSE)
       penalty.weights <- weights*pmin(1, penalty * 20)
       bPWM <- weighted.median(bIV, penalty.weights)
       sebPWM <- mr.boot(bzx, SEzx, bzy, SEzy, penalty.weights)
       res <- metafor::rma(bIV, 1/sqrt(weights), method="FE")
       CochQ <- res$QE
       CochQp <- res$QEp
       graph_title <- paste("SNP effect on", i)
       x_title <- paste("Effect on", X)
       y_title <- paste("Effect on", i)
       if(other_plots)
       {
         metafor::funnel(res, main=paste(graph_title), xlab=y_title)
         abline(v=0, lty="dashed", col="red")
         metafor::forest(res, main=paste(graph_title), addfit=FALSE, slab=SNP, xlab=y_title)
       }
       bzxLCL <- bzx+c*SEzx; bzxUCL <- bzx-c*SEzx
       bzyLCL <- bzy-c*SEzy; bzyUCL <- bzy+c*SEzy
       ES_plot[[i]] <- ggplot2::ggplot(data.frame(bzx, SEzx, bzy, SEzy, bzxLCL, bzxUCL, bzyLCL, bzyUCL),
                                aes(x=bzx, y=bzy)) +
                                geom_point() +
                                theme_bw() +
                                geom_errorbar(aes(ymin=bzyLCL, ymax=bzyUCL)) +
                                geom_errorbarh(aes(xmin=bzxLCL, xmax=bzxUCL)) +
                                geom_abline(intercept=0, slope=0, size=1) +
                                geom_abline(intercept=0, slope=bIVW, size=1, colour="red", show.legend=TRUE) +
                                geom_abline(intercept=0, slope=bEGGER, size=1, colour="blue", show.legend=TRUE) +
                                geom_abline(intercept=0, slope=bWM, size=1, colour="green", show.legend=TRUE) +
                                geom_abline(intercept=0, slope=bPWM, size=1, colour="orange", show.legend=TRUE) +
                                geom_vline(xintercept=0, size=1) +
                                theme_half_open(12) +
                                theme(plot.margin = margin(6, 4, 6, 4)) +
                                ggtitle(graph_title) +
                                xlab(x_title) +
                                ylab(y_title)
       plot(ES_plot[[i]])
       ColOut <- parse(text = paste0(X, ".", i))
       assign(paste(ColOut),c(paste(ColOut),bIVW,sebIVW,CochQ,CochQp,bEGGER,sebEGGER,intEGGER,seintEGGER,bWM,sebWM,bPWM,sebPWM))
       k <- k+1
  }
  r <- c("IV","bIVW","sebIVW","CochQ","CochQp","bEGGER","sebEGGER","intEGGER","seintEGGER","bWM","sebWM","bPWM","sebPWM")
  for (p in nxy) r <- rbind(r, eval(parse(text = paste(p))))
  invisible(list(r=r,plots=ES_plot))
}

library(cowplot)
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
