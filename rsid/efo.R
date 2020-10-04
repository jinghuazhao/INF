library(ontologyIndex)

id <- function(ontology)
{
  inflammatory <- grep(ontology$name,pattern="inflammatory")
  immune <- grep(ontology$name,pattern="immune")
  inf <- union(inflammatory,immune)
  list(id=ontology$id[inf],name=ontology$name[inf])
}
# GO
data(go)
goidname <- id(go)
# EFO
file <- "annotate/efo.obo"
get_relation_names(file)
efo <- get_ontology(file, extract_tags="everything")
length(efo) # 89
length(efo$id) # 27962
efoidname <- id(efo)
diseases <- get_descendants(efo,"EFO:0000408")
efo_0000540 <- get_descendants(efo,"EFO:0000540")
efo_0000540name <- efo$name[efo_0000540]
isd <- data.frame(efo_0000540,efo_0000540name)
save(efo,diseases,isd,efoidname,goidname, file="work/efo.rda")
write.csv(isd,file="work/efo_0000540.csv",col.names=FALSE,row.names=FALSE)
pdf("work/efo_0000540.pdf",height=15,width=15)
library(ontologyPlot)
onto_plot(efo,efo_0000540)
dev.off()

# http://www.ebi.ac.uk/efo/releases/v3.14.0/efo.owl
# http://www.ebi.ac.uk/efo/efo.obo
# Malone J, et al. (2010) Modeling sample variables with an Experimental Factor Ontology. Bioinformatics, 26(8): 1112–1118.
# Fang H, et al. (2019). A genetics-led approach defines the drug target landscape of 30 immune-related traits. Nature Genetics 51(7): 1082-1091.
