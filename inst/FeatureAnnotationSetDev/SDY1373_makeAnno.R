# Evan Henrich
# Oct 2018
# Notes: SDY1373 used ensembl version 79 according to the GSE accession:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97590

library(biomaRt)
library(data.table)

# Get ensembl id and gene id
mart <- useEnsembl(biomart = "ensembl", version = 79, dataset = "hsapiens_gene_ensembl")
ens2gene <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),
             mart = mart)

# Remove unmapped and duplicated probe ids
res <- data.frame(Probe_ID = ens2gene$ensembl_gene_id,
                  Gene_Symbol = ens2gene$hgnc_symbol,
                  stringsAsFactors = FALSE)
res <- res[ !is.na(res$Gene_Symbol), ]
res <- res[ !duplicated(res$Probe_ID), ]

# Write out
write.table(res,
            file = "FeatureAnnotationSetDev/SDY1373_customAnno.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

