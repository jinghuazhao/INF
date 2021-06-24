# 1. https://cran.r-project.org/web/packages/httr/vignettes/quickstart.html
library(httr)
r <- GET("http://httpbin.org/get")
r
status_code(r)
http_status(r)
headers(r)
headers(r)$date
content(r, "text")
cookies(r)
url <- "http://httpbin.org/post"
body <- list(a = 1, b = 2, c = 3)
r <- POST(url, body = body, encode = "form")
r <- POST(url, body = body, encode = "multipart")
r <- POST(url, body = body, encode = "json")
# 
# 2. https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/LitVar/api.html
# 3. http://gwasapi.mrcieu.ac.uk/docs/
#    https://mrcieu.github.io/ieugwasr/articles/guide.html

rsid <- "rs10076557"
IL12B <- subset(pQTLtools::inf1,gene=="IL12B")
r <- with(IL12B,paste0(chr,":",start,"-",end))

library(ieugwasr)
api_status()
batches()
gwasinfo()
tophits(2)
variants_rsid(rsid)
p <- phewas(rsid,pval=5e-8)
data.frame(p)
a <- associations(r, c(2,7))
data.frame(a)
