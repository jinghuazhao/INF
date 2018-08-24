# log10(p) and p for a very large z from Normal(0,1)

First thing first, here is the anwser for log10(p) given z,
```r
# -log10(p) for a standard normal deviate z based on log()
-log(2, base=10)-pnorm(-abs(z), lower.tail=TRUE, log=TRUE)/log(10)
```
We start with z=1.96 whose corresponding p value is known to be 0.05 approximately.
```r
2*pnorm(-1.96,lower.tail=TRUE)
```
or 0.04999579. We proceed to get -log10(p)
```r
-log10(2)-log10(pnorm(-abs(z),lower.tail=TRUE))
```
leading to form above from the fact that log10(X)=log(X)/log(10) since log() is 
also the natural log, ln(), in R, works far better on the numerator of the second term.
The use of -abs() simply makes sure we are working on the lower tail of the standard
Normal distribution from which our p value is calculated.

For a stress test, now let
```r
z <- 20000
-log10(2)-pnorm(-abs(z), lower.tail=TRUE, log=TRUE)/log(10)
```
giving -log10(p) = 86858901.

We would be curious about the p value itself as well, which is furnished together with log10(p) via Rmpfr package
```r
require(Rmpfr)
format(2*pnorm(mpfr(-abs(z),100),lower.tail=TRUE,log.p=FALSE))
format(2*pnorm(mpfr(-abs(z),100),lower.tail=TRUE,log.p=TRUE))
```
giving 1.660579603192917090365313727164e-86858901 and -400000021.6448521764816015432890, respectively.

We can conclude that the order of magnitude is the same 86858901 even with such a big z, the base
making the slight difference. In reality, there might well be no such a huge z score!
