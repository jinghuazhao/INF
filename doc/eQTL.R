# 17-12-2019 JHZ

m <- 6
designs <- 1:6
N <- 50 * designs
n.grids <- 100

index <- 1:n.grids
grids <- index / n.grids

MAF <- c(0.5, 1.0, 2.0, 5.0, 10.0, 20.0)/100
MAF <- seq(0.005, n.grids/2, by=0.5) / n.grids

require(powerEQTL)
png("eQTL.png", res=300, height=8, width=6, units="in")
plot(MAF,grids,type="n",ylab="Power")
title(main="Power Estimation for eQTL Studies", cex.main = 0.9,
      sub="alpha=0.05, effect size=0.8, nSNP=240 (one−way unbalanced ANOVA)", cex.sub=0.6)
colors <- hcl.colors(m)
for (design in designs)
{
  power <- rep(NA,n.grids)
  for (j in index) power[j] <- powerEQTL.ANOVA2(effsize = 0.8,
                                                MAF = MAF[j], typeI = 0.05, nTests = 240,
                                                myntotal = N[design], verbose = FALSE)
  lines(MAF,power,col=colors[design])
}
legend("topleft", inset=.02, title="Sample size (N)", paste(N), col=colors, horiz=FALSE, cex=0.8, lty=designs)
dev.off()
