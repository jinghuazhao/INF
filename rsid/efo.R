library(ontologyIndex)
file <- "annotate/efo.obo"
get_relation_names(file)
onto <- get_ontology(file, extract_tags="everything")
# 89
length(onto)
# 27962
length(onto$id)
inflammatory <- grep(onto$name,pattern="inflammatory")
immune <- grep(onto$name,pattern="immune")
inf <- union(inflammatory,immune)
ontoname <- onto$name[inf]
# 115
length(ontoname)
ontoid <- onto$id[inf]
# 115
length(ontoid)
save(ontoid, file="work/efo.rda")

# http://www.ebi.ac.uk/efo/releases/v3.14.0/efo.owl
# http://www.ebi.ac.uk/efo/efo.obo
# Malone J, et al. (2010) Modeling sample variables with an Experimental Factor Ontology. Bioinformatics, 26(8): 1112–1118
