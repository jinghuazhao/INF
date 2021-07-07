#' Truncated Manhattan plot
#'
#' To generate truncated Manhattan plot, e.g., of genomewide significance (P values) or a random variable that is uniformly distributed.
#'
#' The rationale of this function is to extend mhtplot() to handle extremely small p values as often seen from a protein GWAS; for R will break down when p <= 1e-324.
#'
#' @param x A data.frame.
#' @param chr Chromosome.
#' @param bp Position.
#' @param p p values, e.g., "1.23e-600".
#' @param log10p log10(p).
#' @param z z statistic, i.e., BETA/SE.
#' @param snp SNP. Pending on the setup it could either of variant or gene ID(s).
#' @param col Colours.
#' @param chrlabs Chromosome labels, 1,2,...22,23,24,25.
#' @param suggestiveline Suggestive line.
#' @param genomewideline Genomewide line.
#' @param highlight A list of SNPs to be highlighted.
#' @param annotatelog10P Threshold of -log10(P) to annotate.
#' @param annotateTop Annotate top.
#' @param cex.mtext axis label extension factor.
#' @param cex.text SNP label extension factor.
#' @param mtext.line position of the y lab.
#' @param cex.y y axis numbers.
#' @param y.ax.space interval of ticks of the y axis.
#' @param y.brk1 lower -log10(P) break point.
#' @param y.brk2 upper -log10(P) break point.
#' @param delta a value to enable column(s) of red points.
#' @param ... other options.
#' @return The plot is shown on or saved to the appropriate device.
#' @examples
#' \dontrun{
#' options(width=120)
#' require(gap.datasets)
#' mhtdata <- within(mhtdata, {z=qnorm(p/2, lower.tail=FALSE)})
#' mhtplot.trunc(mhtdata, chr = "chr", bp = "pos", z = "z", snp = "rsn",
#'               y.brk1=6, y.brk2=10, y.ax.space=1, mtext.line=2.5)
#' # https://portals.broadinstitute.org/collaboration/
#' # giant/images/0/0f/Meta-analysis_Locke_et_al+UKBiobank_2018.txt.gz
#' gz <- gzfile("work/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz")
#' BMI <- within(read.delim(gz,as.is=TRUE), {Z <- BETA/SE})
#' print(subset(BMI[c("CHR","POS","SNP","P")],CHR!=16 & P<=1e-150))
#' library(Rmpfr)
#' print(within(subset(BMI, P==0, select=c(CHR,POS,SNP,Z)),
#'              {P <- format(2*pnorm(mpfr(abs(Z),100),lower.tail=FALSE));
#'               Pvalue <- pvalue(Z); log10P <- -log10p(Z)}))
#' png("BMI.png", res=300, units="in", width=9, height=6)
#' par(oma=c(0,0,0,0), mar=c(5,6.5,1,1))
#' mhtplot.trunc(BMI, chr="CHR", bp="POS", z="Z", snp="SNP",
#'               suggestiveline=FALSE, genomewideline=-log10(1e-8),
#'               cex.mtext=1.2, cex.text=1.2,
#'               annotatelog10P=156, annotateTop = FALSE,
#'               highlight=c("rs13021737","rs17817449","rs6567160"),
#'               mtext.line=3, y.brk1=200, y.brk2=280, cex.axis=1.2, cex.y=1.2, cex=0.5,
#'               y.ax.space=20,
#'               col = c("blue4", "skyblue")
#' )
#' dev.off()
#' }
#' @author James Peters, Jing Hua Zhao.
#' @keywords hplot.
#' @seealso \code{\link[gap]{mhtplot}}.
#' @export

mhtplot.trunc <- function (x, chr = "CHR", bp = "BP", p = NULL, log10p = NULL, z = NULL, snp = "SNP",
                           col = c("gray10", "gray60"),
                           chrlabs = NULL, suggestiveline = -log10(1e-05),
                           genomewideline = -log10(5e-08), highlight = NULL,
                           annotatelog10P = NULL, annotateTop = FALSE, cex.mtext=1.5, cex.text=0.7,
                           mtext.line = 2, cex.y = 1, y.ax.space = 5, y.brk1, y.brk2, delta=0.05, ...)
{
  for (q in c("calibrate","plotrix")) {
     if (length(grep(paste("^package:", q, "$", sep=""), search())) == 0) {
        if (!requireNamespace(q, quietly = TRUE))
        warning(paste("mhtplot.trunc needs package `", q, "' to be fully functional; please install", sep=""))
     }
  }
  CHR <- BP <- BP.x <- BP.y <- SNP <- log10P <- index <- NULL
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!is.numeric(x[[chr]])) 
     stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!(bp %in% names(x)))  stop(paste("Column", bp, "not found!"))
  if (!is.null(p)) log10P <- -log10pvalue(x[[p]])
  if (is.null(p) & !is.null(log10p)) log10P <- -as.numeric(x[[log10p]])
  if (is.null(p) & is.null(log10p) & !is.null(z)) log10P <- -log10p(as.numeric(x[[z]]))
  if(is.null(p) & is.null(log10p) & is.null(z)) stop("At least one of p, log10p, or z (priority given in that order) is needed")
  if (y.brk2 <= y.brk1) stop("y.brk2 must be larger than y.brk1")
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  d <- data.frame(CHR = x[[chr]], BP = as.integer(x[[bp]]), log10P = log10P)
  d <- subset(d, !is.na(CHR) & !is.na(BP) & !is.na(log10P))
  if (!is.null(x[[snp]])) d <- transform(d, SNP = x[[snp]])
  d.order <- with(d,order(CHR, BP))
  d <- within(d[d.order,], {pos <- NA; index <- NA})
  ind <- 0
  for (i in unique(with(d,CHR))) {
    ind <- ind + 1
    d[with(d,CHR) == i, "index"] <- ind
  }
  nchr <- length(unique(with(d,CHR)))
  if (nchr == 1) {
    d["pos"] <- d["BP"]
    ticks <- floor(length(d["pos"]))/2 + 1
    xlabel <- paste("Chromosome", unique(with(d,CHR)), "position")
    labs <- ticks
  } else {
    lastbase <- 0
    ticks <- NULL
    for (i in unique(with(d,index))) {
      if (i == 1) d[with(d,index) == i, "pos"] <- d[with(d,index) == i, "BP"]
      else {
        lastbase = lastbase + tail(with(subset(d, index == i - 1),BP), 1)
        d[with(d,index) == i, "pos"] = d[with(d,index) == i, "BP"] + lastbase
      }
      ticks <- c(ticks, (min(d[with(d,index) == i, "pos"]) + max(d[with(d,index) == i, "pos"]))/2 + 1)
    }
    xlabel <- "Chromosome"
    labs <- unique(with(d,CHR))
  }
  xmax <- ceiling(max(with(d,pos)) * 1.03)
  xmin <- floor(max(with(d,pos)) * -0.03)
  max.y <- ceiling(max(with(d,log10P), na.rm=TRUE))
  if (y.brk2 > max.y ){
      message(paste("max.y is", max.y))
      stop("User error: Upper breakpoint must be lower than maximum -log10(P-value)")
  }
  offset <- y.brk2-y.brk1
  d <- within(d, {
    gapped <- log10P > y.brk1 & log10P < y.brk2
    above <- log10P > y.brk2
    log10P[gapped] <- NA
    log10P[above] <- log10P[above] - offset
  })
  def_args <- list(xaxt = "n", yaxt="n", bty = "n", xaxs = "i",
                   las = 1, pch = 20, xlim = c(xmin, xmax),
                   ylim = c(0, ceiling(max(with(d,log10P), na.rm=TRUE))),
                   xlab = "", ylab = "")
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  mtext(text = xlabel, side = 1, line = mtext.line, cex = cex.mtext, font=2)
  mtext(text = expression(-log[10](italic(p))), side=2, line = mtext.line, cex = cex.mtext, font=2)
  y.lab.tick.pos <- seq(from = 0, by = y.ax.space, to = ceiling(max.y) - offset + y.ax.space/3)
  pre.brk.labs <- seq(from = 0, by = y.ax.space, to = y.brk1)
  y.labels <- c(pre.brk.labs, seq(from=y.brk2, by=y.ax.space, length.out=length(y.lab.tick.pos)-length(pre.brk.labs)))
  axis(side=2, at=y.lab.tick.pos, labels=y.labels, cex.axis=cex.y, las=1)
  plotrix::axis.break(axis = 2, breakpos = y.brk1, style = "slash")
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
       if (length(chrlabs) == length(labs)) labs <- chrlabs
       else warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
    }
    else warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
  }
  if (nchr == 1) axis(1, ...) else axis(1, at = ticks, labels = labs)
  col <- rep(col, max(with(d,CHR)))
  if (nchr == 1) with(d, points(pos, log10P, pch = 20, col = col[1], ...))
  else {
    icol = 1
    for (i in unique(with(d,index))) {
      with(d[with(d,index) == unique(with(d,index))[i], ], points(pos, log10P, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline) abline(h = suggestiveline, col = "blue")
  if (genomewideline) abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% with(d,SNP)))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(with(d,SNP) %in% highlight), ]
    with(d.highlight, points(pos, log10P, col = "red", pch = 20, ...))
    d.column <- subset(merge(d,d.highlight[c("CHR","BP")],by=c("CHR")),BP.x>(1-delta)*BP.y & BP.x<(1+delta)*BP.y)
    print(nrow(d.column))
    with(d.column,points(pos, log10P, col = "red", pch = 20, ...))
  }
  if (!is.null(annotatelog10P)) {
    topHits = subset(d, log10P >= annotatelog10P)
    if (!annotateTop) {
      toHighlight <- subset(topHits,SNP %in% highlight)
      toLift <- c("SH2B3","TARF3","FLT3","RAD51B")
      part1 <- subset(toHighlight, ! SNP %in% toLift)
      part2 <- subset(toHighlight, SNP %in% toLift)
      with(part1, calibrate::textxy(pos, log10P, offset = 0.625, pos = 3, labs = SNP, cex = cex.text, font = 4))
      with(part2, calibrate::textxy(pos, log10P+10, offset = 0.625, pos = 3, labs = SNP, cex = cex.text, font = 4))
    }
    else {
      topHits <- topHits[order(with(topHits,log10P)), ]
      topSNPs <- NULL
      for (i in unique(with(topHits,CHR))) {
        chrSNPs <- topHits[with(topHits,CHR) == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      with(topSNPs,calibrate::textxy(pos, log10P, offset = 0.625, pos = 3, labs = SNP, cex = cex.text, font = 4))
    }
  }
}
