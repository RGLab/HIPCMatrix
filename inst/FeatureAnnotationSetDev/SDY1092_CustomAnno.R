# SDY1092
# February 2019
# Evan Henrich

# Download raw counts from GEO
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97158&format=file&file=GSE97158%5Fraw%5Fcounts%2Etsv%2Egz"
fl <- file.path(getwd(), "FeatureAnnotationSetDev", "GSE97158.tsv.gz")
download.file(url = url, destfile = fl)

# Gunzip
GEOquery::gunzip(fl)

# read in orig raw counts file
df <- data.table::fread(gsub("\\.gz","",fl), header = TRUE)

# Remove unmapped if necessary
res <- data.frame(Probe_ID = df$gene_symbol,
                  Gene_Symbol = df$gene_symbol,
                  stringsAsFactors = F)
res <- res[ !is.na(res$Gene_Symbol), ]
res <- res[ !duplicated(res$Probe_ID), ]

# Write out
write.table(res,
            file = "FeatureAnnotationSetDev/SDY1092_CustomAnno.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

# Remove orig file
file.remove(gsub("\\.gz","",fl))
