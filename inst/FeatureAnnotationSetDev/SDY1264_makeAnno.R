# Evan Henrich
# July 2018
# Notes: SDY1264 is a study done by Bali Pulendran at Stanford and used in the ImmuneSignatures2 project.
# Because HIPC collaborators only uploaded a GEF using GEO accessions, the annotation came from notes in GEO
# see - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL7567.
# At issue is that about 1/4 of the probes / unigene values are not mapped to official gene symbols.

library(org.Hs.eg.db)
library(data.table)

# get unigene from norm_exprs by cutting off suffix of probes
em <- data.table::fread("FeatureAnnotationSetDev/SDY1264_Trial1.tsv")
unigenes <- gsub("_at", "", em$feature_id)

# Get unigene-2-symbol map from org.Hs.eg.db
std <- suppressMessages(data.table(AnnotationDbi::select(org.Hs.eg.db,
                                                         keys = AnnotationDbi::keys(org.Hs.eg.db,
                                                                                    keytype = "SYMBOL"),
                                                         columns = c("UNIGENE", "SYMBOL"),
                                                         keytype = "SYMBOL")))
syms <- std$SYMBOL[ match(unigenes, std$UNIGENE)]

# Remove unmapped
res <- data.frame(Probe_ID = em$feature_id,
                  Gene_Symbol = syms,
                  stringsAsFactors = F)
res <- res[ !is.na(res$Gene_Symbol), ]

# Write out
write.table(res,
              file = "FeatureAnnotationSetDev/SDY1264_customAnno.tsv",
              quote = FALSE,
              sep = "\t",
              row.names = FALSE)

# Remove unmapped
res <- data.frame(Probe_ID = std$UNIGENE,
                  Gene_Symbol = std$SYMBOL,
                  stringsAsFactors = F)
res <- res[ !is.na(res$Gene_Symbol), ]

# Write out
write.table(res,
            file = "FeatureAnnotationSetDev/unigene2symbol.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

