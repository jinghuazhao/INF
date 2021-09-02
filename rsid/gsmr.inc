# Adapted from code by Felix Day 16/9/2015
# with reference from Bowden J, et al. (2015). 44(2):512-525.

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

#' Mendelian randomization analysis
#' 
#' The function initially intends to rework on GSMR outputs, but it would be appropriate for general use.
#' 
#' @md
#' @param data Data to be used.
#' @param X Exposure.
#' @param Y Outcome.
#' @param alpha type I error rate for confidence intervals.
#' @param other_plots To add funnel and forest plots.
#' @export
#' @return The result and plots.
#' @examples
#' library(cowplot)
#' library(ggplot2)
#' library(gap)
#' r <- gsmr(mr, "LIF.R", c("CAD","FEV1"))

gsmr <- function(data, X, Y, alpha=0.05, other_plots=FALSE)
{
   c <- qnorm(alpha/2,lower.tail=FALSE)
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
       x_title <- paste("Effect on", gsub("[.]","-",X))
       y_title <- paste("Effect on", i)
       if(other_plots)
       {
         metafor::funnel(res, main=paste(graph_title), xlab=y_title)
         abline(v=0, lty="dashed", col="red")
         metafor::forest(res, main=paste(graph_title), addfit=FALSE, slab=SNP, xlab=y_title)
       }
       bzxLCL <- bzx-c*SEzx; bzxUCL <- bzx+c*SEzx
       bzyLCL <- bzy-c*SEzy; bzyUCL <- bzy+c*SEzy
       intercept <- slope <- colour <- method <- NA
       group4 <- data.frame(intercept=c(0,intEGGER,0,0),
                            slope=c(bIVW,bEGGER,bWM,bPWM),
                            colour=c("red","blue","green","orange"),
                            method=c("IVW","Egger","WM","PWM"))
       ES_plot[[i]] <- ggplot2::ggplot(data.frame(bzx, SEzx, bzy, SEzy, bzxLCL, bzxUCL, bzyLCL, bzyUCL),
                       ggplot2::aes(x=bzx, y=bzy)) +
                       ggplot2::geom_point() +
                       ggplot2::theme_bw() +
                       cowplot::theme_cowplot(12) +
                       ggplot2::geom_errorbar(ggplot2::aes(ymin=bzyLCL, ymax=bzyUCL)) +
                       ggplot2::geom_errorbarh(ggplot2::aes(xmin=bzxLCL, xmax=bzxUCL)) +
                       ggplot2::geom_abline(data=group4, ggplot2::aes(intercept=intercept, slope=slope, colour=method), size=1, show.legend=TRUE) +
                       ggplot2::scale_colour_manual(values = with(group4,colour)) +
                       ggplot2::theme(legend.position="bottom", legend.direction="vertical") +
                       ggplot2::guides(colour=ggplot2::guide_legend(ncol=4)) +
                       ggplot2::geom_abline(intercept=0, slope=0, size=1) +
                       ggplot2::geom_vline(xintercept=0, size=1) +
                       ggplot2::ggtitle(graph_title) +
                       ggplot2::xlab(x_title) +
                       ggplot2::ylab(y_title)
       plot(ES_plot[[i]])
       ColOut <- parse(text = paste0(X, ".", i))
       assign(paste(ColOut),c(paste(ColOut),bIVW,sebIVW,CochQ,CochQp,bEGGER,sebEGGER,intEGGER,seintEGGER,bWM,sebWM,bPWM,sebPWM))
       k <- k+1
  }
  r <- c("IV","bIVW","sebIVW","CochQ","CochQp","bEGGER","sebEGGER","intEGGER","seintEGGER","bWM","sebWM","bPWM","sebPWM")
  for (p in nxy) r <- rbind(r, eval(parse(text = paste(p))))
  invisible(list(r=r,plots=ES_plot))
}
