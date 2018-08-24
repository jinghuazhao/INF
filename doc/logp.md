# p and log10(p) for very large z from Normal(0,1)

First thing first, here is the anwser,
```r
# -log10(p) for a standard normal deviate z based on log()
logp <- -log(2, base=10)-pnorm(-abs(z), lower.tail=TRUE, log=TRUE)/log(10)
```
We start with z=1.96 whose corresponding p value is known approximately 0.05.
```r
2*pnorm(-1.96,lower.tail=TRUE)
```
or 0.04999579. We proceed with the log version
```r
log10(2)+log10(pnorm(z,lower.tail=TRUE))
```
leading to form above from the fact that log10(X)=ln(X)/ln(10) since ln(), or 
equivalently log() in R, works far better on the numerator of the second term.

To test, now let
```r
z. <- 20000
-log10(2)+pnorm(-abs(z), log=TRUE)/log(10)
```
giving -log10(p)=86858901.

To contrast with Rmpfr package on the actual p and log10(p),
```r
require(Rmpfr)
format(2*pnorm(mpfr(-abs(z),100),lower.tail=TRUE,log.p=FALSE))
format(2*pnorm(mpfr(-abs(z),100),lower.tail=TRUE,log.p=TRUE))
```
which are 1.660579603192917090365313727164e-86858901 and -400000021.6448521764816015432890, respectively.

We can conclude that the order of magnitude is the same 86858901 even with such a big z, the base
making the slight difference. In reality, there might well not be such a huge z score!
