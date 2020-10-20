# log(p), log10(p) and p for a very large z from Normal(0,1)

When a Normally distributed association statistic z is very large, its corresponding p value is very small. A genomewide significance is declared at 0.05/1000000=5e-8 with Bonferroni correction assuming 1 million SNPs are tested. This short note describes how to get -log10(p), which can be used in a Q-Q plot and software such as DEPICT. The solution here is generic since z is also the square root of a chi-squared statistic, for instance.

## log(p) and log10(p)
First thing first, here are the answers for log(p) and log10(p) given z,
```r
# log(p) for a standard normal deviate z based on log()
logp <- function(z) log(2)+pnorm(-abs(z), lower.tail=TRUE, log.p=TRUE)

# log10(p) for a standard normal deviate z based on log()
log10p <- function(z) log(2, base=10)+pnorm(-abs(z), lower.tail=TRUE, log.p=TRUE)/log(10)
```
Note `logp()` will be used for functions such as `qnorm()` as in R/gap function `cs()` whereas `log10p()` is more appropriate for Manhattan plot and used in R/gap `sentinels()`.

## Rationale
We start with z=1.96 whose corresponding p value is approximately 0.05.
```r
2*pnorm(-1.96,lower.tail=TRUE)
```
giving an acceptable value 0.04999579, so we proceed to get log10(p)
```r
log10(2)+log10(pnorm(-abs(z),lower.tail=TRUE))
```
leading to the expression above from the fact that log10(X)=log(X)/log(10) since log(),
being the natural log function, ln() -- so log(exp(1)) = 1, in R, works far better on
the numerator of the second term. The use of -abs() just makes sure we are working on
the lower tail of the standard Normal distribution from which our p value is calculated.

## Benchmark
Now we have a stress test,
```r
z <- 20000
-log10p(z)
```
giving -log10(p) = 86858901.

## p, log(p), log10(p) and the multiple precision arithmetic

We would be curious about the p value itself as well, which is furnished with the Rmpfr package
```r
require(Rmpfr)
2*pnorm(mpfr(-abs(z),100),lower.tail=TRUE,log.p=FALSE)
mpfr(log(2),100) + pnorm(mpfr(-abs(z),100),lower.tail=TRUE,log.p=TRUE)
```
giving p = 1.660579603192917090365313727164e-86858901 and -log(p) = -200000010.1292789076808554854177,
respectively. To carry on we have -log10(p) = -log(p)/log(10)=86858901.

To make -log10(p) usable in R we obtain it directly through
```r
as.numeric(-log10(2*pnorm(mpfr(-abs(z),100),lower.tail=TRUE)))
```
which actually yields exactly the same 86858901.

If we go very far to have z=50000. then -log10(p)=542868107 but we have less luck with Rmpfr.

One may wonder the P value in this case, which is 6.6666145952e-542868108 or simply 6.67e-542868108.

The magic function for doing this is defined as follows,
```r
pvalue <- function(z,decimals=2)
{
  lp <- log10p(z)
  exponent <- ceiling(lp)
  base <- 10^(lp - exponent)
  paste0(round(base*10,decimals),"e",-1+exponent)
}
```
