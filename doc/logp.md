# -log10(p) and p for a very large z from Normal(0,1)

## -log10(p)
First thing first, here is the anwser for -log10(p) given z,
```r
# -log10(p) for a standard normal deviate z based on log()
-log(2, base=10)-pnorm(-abs(z), lower.tail=TRUE, log.p=TRUE)/log(10)
```

## Rationale
We start with z=1.96 whose corresponding p value is approximately 0.05.
```r
2*pnorm(-1.96,lower.tail=TRUE)
```
giving an acceptable value 0.04999579, so we proceed to get -log10(p)
```r
-log10(2)-log10(pnorm(-abs(z),lower.tail=TRUE))
```
leading to the expression above from the fact that log10(X)=log(X)/log(10) since log(),
being the natural log function, ln() -- so log(exp(1)) = 1, in R, works far better on
the numerator of the second term. The use of -abs() just makes sure we are working on
the lower tail of the standard Normal distribution from which our p value is calculated.

## Benchmark
Now we have a stress test,
```r
z <- 20000
-log10(2)-pnorm(-abs(z), lower.tail=TRUE, log=TRUE)/log(10)
```
giving -log10(p) = 86858901.

## p, -log10(p) and the multiple precision arithmetic

We would be curious about the p value itself as well, which is furnished with the Rmpfr package
```r
require(Rmpfr)
format(2*pnorm(mpfr(-abs(z),100),lower.tail=TRUE,log.p=FALSE))
format(2*pnorm(mpfr(-abs(z),100),lower.tail=TRUE,log.p=TRUE))
```
which gives p = 1.660579603192917090365313727164e-86858901 and -log10(p) = 400000021.6448521764816015432890,
respectively.

We can then conclude that the order of magnitude is the same 86858901 even with a z big as this, the
base making the slight difference. To make -log10(p) usable in R we obtain it directly through
```r
as.numeric(-log10(2*pnorm(mpfr(-abs(z),100),lower.tail=TRUE)))
```
which actually yields exactly the same 86858901.

If we go very far to have z=50000. then -log10(p)=542868107 but we have less luck with Rmpfr.
