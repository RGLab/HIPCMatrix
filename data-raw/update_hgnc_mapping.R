# Updates HGNC mapping table and gene signatures with new annotation.

library(HIPCMatrix)

# First, update hgnc mapping table
hgncAlias2Symbol <- create_gene_alias_map()
hgncAlias2Symbol_version <- Sys.Date()

usethis::use_data(
  hgncAlias2Symbol,
  hgncAlias2Symbol_version,
  overwrite = TRUE
)


# Use mapping to update gene sets
updated_btm_df <- update_geneset_symbols(
  "orig_btm_list",
  hgncAlias2Symbol
)
updated_btm_list <- plyr::dlply(updated_btm_df, 1, function(x) x$SYMBOL)

chaussabel_modules <- plyr::dlply(
  update_geneset_symbols(
    "orig_chaussabel",
    hgncAlias2Symbol
  ),
  1,
  function(x) x$SYMBOL
)

emory_blood_transcript_modules <- plyr::dlply(
  update_geneset_symbols(
    "orig_emory",
    hgncAlias2Symbol
  ),
  1,
  function(x) x$SYMBOL
)

msigdb_immunologic_signatures <- plyr::dlply(
  update_geneset_symbols(
    "orig_msigdb",
    hgncAlias2Symbol
  ),
  1,
  function(x) x$SYMBOL
)

usethis::use_data(
  updated_btm_df,
  updated_btm_list,
  chaussabel_modules,
  emory_blood_transcript_modules,
  msigdb_immunologic_signatures,
  overwrite = TRUE
)
