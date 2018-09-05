# http://jokergoo.github.io/circlize/example/genomic_link.html

library(circlize)

bed1 = generateRandomBed(nr = 100)
bed1 = bed1[sample(nrow(bed1), 1), ]
bed1 <- rbind(bed1,bed1,bed1,bed1,bed1,bed1)
bed2 = generateRandomBed(nr = 100)
bed2 = bed2[sample(nrow(bed2), 6), ]
circos.par("track.height" = 0.1, cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram()

circos.genomicLink(bed1, bed2, col = sample(1:6, 6, replace = TRUE), border = NA)
circos.clear()
