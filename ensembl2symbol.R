# Enseble2symbols.R
# sometimes bioMart site is unaccessible (e.g. 2 May 2022). To fix this, just specify the following host argument in line 11: host = "dec2021.archive.ensembl.org"
getEnsemblGenes <- function() {
  library(biomaRt)
  getBM(
    attributes = c(
      "external_gene_name", 
      "ensembl_gene_id"
    ),
    #		mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    mart = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
  )
}
ens2sym <- function(ensids) {
  mapper <- getEnsemblGenes()
  v <- unlist(sapply(ensids, function(ensid) {
    paste(unique(mapper$external_gene_name[ensid == mapper$ensembl_gene_id]), collapse = "|")
  }))
  v[v == ""] <- NA
  v
}
