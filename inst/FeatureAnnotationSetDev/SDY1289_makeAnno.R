# E.Henrich
# July 26, 2018
# FeatureAnnotationSet for SDY1289: GSE13699 in GEO
# Platform: Illumina HumanRef-8 v3.0 expression beadchip
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6883

# Downloaded the "full table" from above link
# GPL6883-11606.txt

df <- data.table::fread("FeatureAnnotationSetDev/GPL6883-11606.txt")

# Remove unmapped
res <- data.frame(Probe_ID = df$ID,
                  Gene_Symbol = df$ILMN_Gene,
                  stringsAsFactors = F)
res <- res[ !is.na(res$Gene_Symbol), ]

# Write out
write.table(res,
            file = "FeatureAnnotationSetDev/IlluminaHumanRef8_v3_anno.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
