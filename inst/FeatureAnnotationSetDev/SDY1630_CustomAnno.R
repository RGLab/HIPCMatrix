# SDY1630
# January 2021 (DR36)
# J. Kim

# Download raw counts from GEO
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE133383&format=file&file=GSE133383%5Ffiltered%5FNKcounts%5Ftable%2Ecsv%2Egz"
fl <- file.path(getwd(), "FeatureAnnotationSetDev", "GSE133383.tsv.gz")
download.file(url = url, destfile = fl)

# Gunzip
GEOquery::gunzip(fl)

# read in orig raw counts file
df <- data.table::fread(gsub("\\.gz", "", fl))

# Remove unmapped if necessary
res <- data.frame(
  Probe_ID = df$V1,
  Gene_Symbol = df$V1,
  stringsAsFactors = FALSE
)
res <- res[!is.na(res$Gene_Symbol), ]
res <- res[!duplicated(res$Probe_ID), ]

# Write out
write.table(res,
  file = "FeatureAnnotationSetDev/SDY1630_CustomAnno.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# Remove orig file
file.remove(gsub("\\.gz", "", fl))
