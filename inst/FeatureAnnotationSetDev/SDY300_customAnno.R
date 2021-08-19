# SDY300
# August 2018
# Evan Henrich

# read in orig tpm counts file
df <- data.table::fread("FeatureAnnotationSetDev/SDY300_EXP13588_RNA_seq.703274.tsv")

# Remove unmapped
res <- data.frame(Probe_ID = df$GENE_SYMBOL,
                  Gene_Symbol = df$GENE_SYMBOL,
                  stringsAsFactors = F)
res <- res[ !is.na(res$Gene_Symbol), ]
res <- res[ !duplicated(res$Probe_ID), ] # there were two

# Note: there are likely misnamed probes still in res, e.g. "10-Mar" is likely a forced datetime

# Write out
write.table(res,
            file = "FeatureAnnotationSetDev/SDY300_CustomAnno.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
