library(HIPCMatrix)
hgncAlias2Symbol <- create_gene_alias_map()
hgncAlias2Symbol_version <- Sys.Date()

usethis::use_data(hgncAlias2Symbol, hgncAlias2Symbol_version, overwrite = TRUE)
