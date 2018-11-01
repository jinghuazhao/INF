gc.lambda <- function(p) {
  p <- p[!is.na(p)]
  n <- length(p)
  obs <- qchisq(p,1,lower.tail=FALSE)
  exp <- qchisq(1:n/n,1,lower.tail=FALSE)
  lambda <- median(obs)/median(exp)
  return(lambda)
}
