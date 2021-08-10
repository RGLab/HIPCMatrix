# Helen Miller
# May 2020
# Notes: SDY787 is an RNASeq study and annotation is based off raw
# data held in the supplemental files of the GSE
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113891

library(httr)
library(data.table)

# get the raw file on the GSE record
link <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE113891&format=file&file=GSE113891%5FPT%5Fall%5FCount%2Etsv%2Egz"
gzfl <- tempfile()
GET(link, write_disk(gzfl, overwrite = TRUE))
fl <- tempfile()
GEOquery::gunzip(filename = gzfl, destname = fl, overwrite = TRUE)
raw <- fread(fl)

# Create anno df
anno <- data.table(Probe_ID = raw$GENES,
                   Gene_Symbol = raw$GENES,
                   stringsAsFactors = FALSE)

# Write out
write.table(anno,
            file = "FeatureAnnotationSetDev/SDY787_customAnno.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

