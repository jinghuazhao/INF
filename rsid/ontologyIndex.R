# https://bioconductor.org/packages/release/bioc/vignettes/gwascat/inst/doc/gwascatOnt.html

library(gwascat)
requireNamespace("graph")
makeCurrentGwascat()
data(efo.obo.g)
efo.obo.g
head(graph::nodes(efo.obo.g))
graph::nodeData(efo.obo.g, "EFO:0000540")
ue = graph::ugraph(efo.obo.g)
neighISD = graph::adj(ue, "EFO:0000540")[[1]]
sapply(graph::nodeData(graph::subGraph(neighISD, efo.obo.g)), "[[", "name")
requireNamespace("RBGL")
p = RBGL::sp.between(efo.obo.g, "EFO:0000685", "EFO:0000001")
sapply(graph::nodeData(graph::subGraph(p[[1]]$path_detail, efo.obo.g)), "[[", "name")
data(ebicat37)
names(S4Vectors::mcols(ebicat37))
adinds = grep("autoimmu|inflamm", ebicat37$MAPPED_TRAIT)
adgr = ebicat37[adinds]
adgr
write.table(as.data.frame(adgr),file="docs/gwas_autoinfl.tsv",quote=FALSE,row.names=FALSE,sep="\t")

# A list of 30 excluding Osteoarthritis based on
# http://pi.well.ox.ac.uk:3010/pidb/gateway

fang_list <- subset(read.delim("docs/fang.list.txt"),fullname!="Osteoarthritis")
nd = graph::nodeData(efo.obo.g)
alldef = sapply(nd, function(x) unlist(x[["def"]]))
allnames = sapply(nd, function(x) unlist(x[["name"]]))
alld = sapply(alldef, function(x) if(is.null(x)) return(" ") else x[1])
df = data.frame(id = names(allnames), concept=as.character(allnames), def=unlist(alld))
autoinfl = df[grep("autoimm|inflam", df$def, ignore.case=TRUE),]
requireNamespace("DT")
suppressWarnings({DT::datatable(autoinfl, rownames=FALSE, options=list(pageLength=5))})
inf1 <- data.frame()
for (str in with(fang_list,fullname)) inf1 <- rbind(inf1,df[grep(str, df$def, ignore.case=TRUE),])
options(width=200)
inf1 <- unique(inf1)
write.table(inf1,file="docs/fang.efos.txt",quote=FALSE,row.names=FALSE,sep="\t")

# https://cran.r-project.org/web/packages/ontologyIndex/vignettes/intro-to-ontologyX.html

library(ontologyIndex)
# GO
data(go)
# EFO
file <- "efo.obo"
get_relation_names(file)
efo <- get_ontology(file, extract_tags="everything")

id <- function(ontology)
{
  length(ontology)
  length(ontology$id)
  inf <- grep(ontology$name,pattern="immune|inflammatory")
  data.frame(id=ontology$id[inf],name=ontology$name[inf])
}

goidname <- id(go)
efoidname <- id(efo)
diseases <- get_descendants(efo,"EFO:0000408")
efo_0000540 <- get_descendants(efo,"EFO:0000540")
efo_0000540name <- efo$name[efo_0000540]
isd <- data.frame(efo_0000540,efo_0000540name)

fang_efo <- with(inf1,id)
fang_disease <- data.frame(id=efo$id[fang_efo],name=efo$name[fang_efo])
save(efo,diseases,fang_disease,isd,efoidname,goidname, file="work/efo.rda")
write.table(isd,file="efo_0000540.csv",col.names=FALSE,row.names=FALSE,sep=",")

efo_0000540_plot <- function()
{
  library(ontologyPlot)
  pdf("efo_0000540.pdf",height=15,width=15)
  onto_plot(efo,efo_0000540)
  dev.off()
}

# http://www.ebi.ac.uk/efo/efo.obo
# http://www.ebi.ac.uk/efo/releases/v3.14.0/efo.owl
# Malone J, et al. (2010) Modeling sample variables with an Experimental Factor Ontology. Bioinformatics, 26(8): 1112–1118.
# Fang H, et al. (2019). A genetics-led approach defines the drug target landscape of 30 immune-related traits. Nature Genetics 51(7): 1082-1091.
