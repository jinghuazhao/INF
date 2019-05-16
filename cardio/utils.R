# MAF ~ N

MAF <- seq(0.01,0.1,0.01)
N <- seq(10,100,10)
PROD <- (MAF^2+MAF*(1-MAF))%o%N
colnames(PROD) <- N
rownames(PROD) <- MAF
PROD
