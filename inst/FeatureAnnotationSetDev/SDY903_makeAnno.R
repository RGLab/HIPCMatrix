# Evan Henrich
# September 2019
# Notes: SDY903 is an RNASeq study and annotation is based off raw
# data held in the supplemental files of the GSE
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105884

library(httr)
library(data.table)

# get the raw file on the GSE record
link <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE105884&format=file&file=GSE105884%5FAll%5FTPM%2Etsv%2Egz"
gzfl <- tempfile()
GET(link, write_disk(gzfl, overwrite=TRUE))
fl <- tempfile()
GEOquery::gunzip(filename = gzfl, destname = fl, overwrite = TRUE)
raw <- fread(fl)

# Remove unmapped and duplicated probe ids
res <- data.frame(Probe_ID = raw$V1,
                  Gene_Symbol = raw$V1,
                  stringsAsFactors = F)
res <- res[ !is.na(res$Gene_Symbol), ]
res <- res[ !duplicated(res$Probe_ID), ]

# Write out
write.table(res,
            file = "FeatureAnnotationSetDev/SDY903_customAnno.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

