foreach g in "IL12B" "KITLG" "TNFSF10" {
    insheet using "GREAT/`g'-all.csv", case comma clear
    di "1g'"
    list Ontology Desc BinomP HyperP Genes if BinomP < 0.0001, linesize(200) separator(0)
}
